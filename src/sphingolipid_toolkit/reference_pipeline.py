"""Input/output adapters; scientific decisions live in the reference engine."""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .engine import SphingolipidIdentifier, compatible_frame, config_from_workflow
from .io_utils import parse_ms2_text_to_dataframe, read_precursor_table, write_excel_atomic, require_file
from .result_writer import write_identifications


def run_text_file(config, file_index, logger=None):
    logger = logger or logging.getLogger(__name__)
    text = Path(config.input_dir) / config.text_pattern.format(i=file_index)
    require_file(text,'MS/MS txt file')
    targets = read_precursor_table(Path(config.precursor_dir or config.input_dir) / config.target_pattern.format(i=file_index))
    fragments = parse_ms2_text_to_dataframe(text, encoding=config.encoding)
    rows = []
    with SphingolipidIdentifier(config_from_workflow(config)) as identifier:
        for (scan_id, precursor), group in fragments.groupby(['scan_id','precursor_mz'],sort=False):
            errors = (targets['target'].astype(float) - precursor).abs() / precursor * 1e6
            rt = group['rt'].iloc[0]
            rt_error = (targets['RT'].astype(float)-float(rt)).abs() if pd.notna(rt) else pd.Series(0.,index=targets.index)
            eligible = targets.loc[(errors <= config.precursor_ppm) & (rt_error <= config.feature_rt_tolerance_min)]
            if eligible.empty: continue
            # Keep independently measured target identity; do not replace the
            # MS2 selected-ion mass with a theoretical mass from the old table.
            index = (errors.loc[eligible.index]/config.precursor_ppm + rt_error.loc[eligible.index]/max(config.feature_rt_tolerance_min,1e-12)).idxmin()
            feature = targets.loc[index]
            result = identifier.identify(float(precursor), list(zip(group.fragment_mz,group.fragment_intensity)),
                scan_id=scan_id,rt_minutes=float(rt) if pd.notna(rt) else float(feature.RT))
            for row in result:
                row.update(source_file=text.name, input_target_row=int(index), target_mz_measured=float(feature.target),
                    target_rt_measured=float(feature.RT), signal_intensity_measured=float(feature.Abund))
            frame = compatible_frame(result, abundance=float(feature.Abund))
            if not frame.empty: rows.append(frame)
    data = pd.concat(rows,ignore_index=True) if rows else compatible_frame([])
    path = Path(config.output_dir)/config.result_pattern.format(i=file_index)
    write_identifications(data,path)
    logger.info('SphinGOlipID reference engine: %s candidates exported to %s',len(data),path)
    return path


def annotate_targeted_spectrum(spectrum, feature, config, identifier):
    result = identifier.identify(float(spectrum.precursor_mz), list(zip(spectrum.fragment_mz,spectrum.fragment_intensity)),
        scan_id=spectrum.scan_id,rt_minutes=spectrum.ms2_rt,
        polarity=getattr(spectrum,'polarity','+'),precursor_charge=getattr(spectrum,'precursor_charge',None),
        metadata=getattr(spectrum,'metadata',None))
    data = compatible_frame(result,abundance=None)
    if data.empty:return data
    for key,value in feature.items():data[key]=value
    data['ms2_rt']=spectrum.ms2_rt
    data['ms2_precursor_mz']=spectrum.precursor_mz
    data['ms2_scan_id']=spectrum.scan_id
    data['ms2_file']=str(spectrum.mzml_file)
    data['ms2_file_name']=spectrum.mzml_file.name
    data['collision_energy']=spectrum.collision_energy
    data['ms1_rt']=float(feature['ms1_rt'])
    data['feature_mz']=float(feature['feature_mz'])
    data['ms1_library_mz_error_ppm']=data['ppm_error']
    return data


def iter_prepared_mzml(path, min_relative_intensity=0.002):
    """Use the reference MS1 refiner, retaining measured masses and charge.

    Two streaming passes avoid loading the MS2 run. Only survey arrays are
    retained for the reference algorithm's preceding/adjacent MS1 check.
    """
    import numpy as np
    from pyteomics import mzml
    from ._vendor.lipidgate.ms2.precursor_refinement import MS1Survey
    from ._vendor.lipidgate.ms2.spectrum_preparation import make_refiner, prepare_spectrum
    from .targeted_mzml_pipeline import _extract_precursor, _extract_scan_rt
    surveys = []
    with mzml.MzML(str(path), decode_binary=True) as reader:
        for index, raw in enumerate(reader, 1):
            if raw.get('ms level') != 1:continue
            rt = _extract_scan_rt(raw)
            if rt is None:continue
            masses = np.asarray(raw.get('m/z array', []), dtype=float)
            intensities = np.asarray(raw.get('intensity array', []), dtype=float)
            order = np.argsort(masses, kind='stable')
            surveys.append(MS1Survey(str(raw.get('id','')),f'scan_{index}',rt*60,
                masses[order],intensities[order],'centroid spectrum' in raw))
    refiner = make_refiner(surveys)
    with mzml.MzML(str(path), decode_binary=True) as reader:
        for index, raw in enumerate(reader, 1):
            if raw.get('ms level') != 2:continue
            precursor = _extract_precursor(raw)
            rt = _extract_scan_rt(raw)
            if precursor is None or rt is None:continue
            raw_parent = raw.get('precursorList',{}).get('precursor',[{}])[0]
            ion = raw_parent.get('selectedIonList',{}).get('selectedIon',[{}])[0]
            charge = ion.get('charge state')
            polarity = '+' if 'positive scan' in raw else '-' if 'negative scan' in raw else ''
            if polarity not in ('+','-'):
                if charge is None:continue
                polarity = '+' if int(charge)>0 else '-'
            prepared = prepare_spectrum(scan_id=f'scan_{index}',raw_mz=precursor['precursor_mz'],
                rt_seconds=rt*60,raw_peaks=list(zip(raw.get('m/z array',[]),raw.get('intensity array',[]))),
                polarity=polarity,charge=charge,parent_native=raw_parent.get('spectrumRef',''),
                refiner=refiner,source=path,min_relative_intensity=min_relative_intensity)
            if prepared is not None:
                prepared.metadata['native_scan_id'] = str(raw.get('id',''))
                prepared.metadata['collision_energy'] = precursor.get('collision_energy')
                yield prepared


def annotate_targeted_directory(features,config,logger):
    from .targeted_mzml_pipeline import _find_input_files,find_matching_feature
    rows=[]
    with SphingolipidIdentifier(config_from_workflow(config)) as identifier:
        for path in _find_input_files(config.mzml_dir,'.mzml'):
            logger.info('Reading mzML file with reference engine: %s',path)
            for spectrum in iter_prepared_mzml(path):
                # The feature must belong to this acquisition. A target from
                # another file is not evidence of local detection.
                local=features
                if 'source_file' in features:
                    local=features.loc[features['source_file'].astype(str).map(lambda value:Path(value).stem==path.stem)]
                feature=find_matching_feature(local,spectrum.precursor_mz,spectrum.rt_minutes,
                    config.feature_mz_ppm,config.feature_rt_tolerance_min)
                if feature is None:continue
                data=compatible_frame(identifier.score_spectrum(spectrum))
                if data.empty:continue
                for key,value in feature.items():data[key]=value
                data['ms2_rt']=spectrum.rt_minutes
                data['ms2_precursor_mz']=spectrum.precursor_mz
                data['ms2_scan_id']=spectrum.metadata['native_scan_id']
                data['ms2_file']=str(path)
                data['ms2_file_name']=path.name
                data['collision_energy']=spectrum.metadata['collision_energy']
                data['ms1_rt']=float(feature['ms1_rt'])
                data['feature_mz']=float(feature['feature_mz'])
                rows.append(data)
    result=pd.concat(rows,ignore_index=True) if rows else compatible_frame([])
    if not result.empty:result.insert(0,'identification_id',[f'ID{i:08d}' for i in range(1,len(result)+1)])
    return result
