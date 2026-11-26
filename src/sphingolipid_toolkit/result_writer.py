"""Preserve complete MS2 evidence when a JSON cell exceeds Excel's limit."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile, TemporaryDirectory
from collections.abc import Mapping

from openpyxl import load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
import pandas as pd

from .io_utils import identification_excel_writer, write_excel_atomic


SCAN_ID_DTYPES = {name: str for name in (
    'scan_id', 'native_scan_id', 'apex_native_scan_id', 'ms2_scan_id', 'preceding_ms1_native_id',
)}


def _xlsx_destination(path):
    path = Path(path)
    if path.suffix.lower() != '.xlsx':
        raise ValueError('结果文件扩展名须为 .xlsx。')
    return path


def _prepare_evidence_tables(tables, path):
    path=Path(path)
    # Workbooks omit the DataFrame index. Give each row a unique position so
    # duplicate labels cannot attach one candidate's evidence to another.
    prepared={name: data.reset_index(drop=True) for name,data in tables.items()}
    evidence={}
    external_indices={name: [] for name in tables}
    for name,output in prepared.items():
        if 'ms2_evidence_json' not in output:
            continue
        for index,value in output['ms2_evidence_json'].items():
            if isinstance(value,str) and len(value)>32767:
                key=hashlib.sha256(value.encode('utf-8')).hexdigest()
                evidence[key]=value
                output.loc[index,'ms2_evidence_json']=json.dumps({'external_evidence_sha256':key},separators=(',',':'))
                output.loc[index,'evidence_sha256']=key
                external_indices[name].append(index)
    created=None
    if evidence:
        payload=json.dumps(evidence,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode('utf-8')
        digest=hashlib.sha256(payload).hexdigest()
        sidecar=path.parent/f'ms2-evidence-{digest[:16]}.json.gz'
        path.parent.mkdir(parents=True,exist_ok=True)
        if sidecar.exists():
            if gzip.decompress(sidecar.read_bytes())!=payload:
                raise ValueError('An existing evidence file has different contents')
        else:
            with sidecar.open('xb') as handle:handle.write(gzip.compress(payload,mtime=0))
            created=sidecar
        for name,indices in external_indices.items():
            if indices:
                prepared[name].loc[prepared[name].index.isin(indices),'evidence_file']=sidecar.name
    return prepared,created


def write_identifications(data,path):
    path=_xlsx_destination(path)
    prepared,created=_prepare_evidence_tables({'Sheet1': data},path)
    try:
        write_excel_atomic(prepared['Sheet1'],path,index=False,literal_strings=True)
    except BaseException:
        if created is not None and created.resolve().parent==path.resolve().parent:
            created.unlink(missing_ok=True)
        raise
    return path


def write_identification_workbook(tables, path):
    """Atomically export all result tables with one complete shared evidence file."""
    path=_xlsx_destination(path)
    prepared,created=_prepare_evidence_tables(tables,path)
    temporary=None
    try:
        path.parent.mkdir(parents=True,exist_ok=True)
        with NamedTemporaryFile(dir=path.parent,prefix='.sphingolipid_result_',suffix='.xlsx',delete=False) as handle:
            temporary=Path(handle.name)
        with identification_excel_writer(temporary) as writer:
            for name,data in prepared.items():
                data.to_excel(writer,sheet_name=name,index=False)
        temporary.replace(path)
    except BaseException:
        if created is not None:
            created.unlink(missing_ok=True)
        raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path


def _file_digest(path: Path) -> bytes:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.digest()


def workbook_stamp(path: Path) -> tuple[int, int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ino


def _write_report_sheets(workbook, tables: Mapping[str, pd.DataFrame]) -> None:
    for position, (name, data) in enumerate(tables.items()):
        if name in workbook:
            del workbook[name]
        if len(data) >= 1_048_576 or len(data.columns) > 16_384:
            raise ValueError(f'工作表超出 Excel 的行列限制：{name}')
        sheet = workbook.create_sheet(name, index=position)
        sheet.append(list(data.columns))
        for cell in sheet[1]:
            cell.data_type = 's'
            cell.font = Font(bold=True)
            sheet.column_dimensions[get_column_letter(cell.column)].width = min(44, max(16, len(cell.value) + 3))
        for row_number, row in enumerate(data.itertuples(index=False, name=None), 2):
            sheet.append([None if pd.isna(value) else value for value in row])
            for cell in sheet[row_number]:
                if cell.data_type in {'f', 'e'}:
                    cell.data_type = 's'
        sheet.freeze_panes = 'B2' if name == 'identification_table' else 'A2'
        sheet.auto_filter.ref = sheet.dimensions
    workbook.active = 0


def copy_identification_bundle(source: str | Path, destination: str | Path, *,
                               report_tables: Mapping[str, pd.DataFrame] | None = None,
                               expected_stamp: tuple[int, int, int] | None = None) -> Path:
    """Save all worksheets and referenced evidence, publishing the XLSX last.

    Safe relative references keep the workbook byte-for-byte intact. Absolute
    references or references outside the destination folder are localized;
    only those reference cells are changed. Derived report worksheets can be added
    in the same transaction as the copy.
    """
    source, destination = Path(source), _xlsx_destination(destination)
    if source.resolve() == destination.resolve():
        if report_tables:
            raise ValueError('鉴定表须保存为新文件。')
        return destination
    if expected_stamp is not None and workbook_stamp(source) != expected_stamp:
        raise ValueError('来源工作簿已更新，请重新打开 ECN 复核与导出。')
    root = destination.parent.resolve()
    evidence = {}
    digests = {}
    rebases = []
    workbook = load_workbook(source, read_only=True)
    try:
        for sheet in workbook.worksheets:
            rows = sheet.iter_rows()
            header = next(rows, ())
            columns = [index for index, cell in enumerate(header) if cell.value == 'evidence_file']
            if not columns:
                continue
            for row in rows:
                for column in columns:
                    cell = row[column]
                    if cell.value is None or not str(cell.value).strip():
                        continue
                    reference = Path(str(cell.value))
                    origin = (source.parent / reference).resolve()
                    if not origin.is_file():
                        raise FileNotFoundError(f'缺少碎片证据文件：{origin}')
                    if origin not in digests:
                        digests[origin] = _file_digest(origin)
                    digest = digests[origin]
                    target = (root / reference).resolve()
                    if reference.is_absolute() or not target.is_relative_to(root):
                        target = root / f'ms2-evidence-{digest.hex()[:16]}.json.gz'
                        rebases.append((sheet.title, cell.coordinate, target.name))
                    if target == destination.resolve():
                        raise ValueError('结果文件不能覆盖其碎片证据文件。')
                    if target in evidence and evidence[target][1] != digest:
                        raise FileExistsError(f'碎片证据文件名称冲突：{target}')
                    if target.exists():
                        if target not in digests:
                            digests[target] = _file_digest(target)
                        if digests[target] != digest:
                            raise FileExistsError(f'已有不同内容的碎片证据文件，未覆盖：{target}')
                    evidence[target] = (origin, digest)
    finally:
        workbook.close()

    destination.parent.mkdir(parents=True, exist_ok=True)
    created = []
    try:
        with TemporaryDirectory(prefix='.sphingolipid_export_', dir=destination.parent) as temporary:
            staged = Path(temporary) / destination.name
            shutil.copy2(source, staged)
            if rebases or report_tables:
                workbook = load_workbook(staged)
                try:
                    for sheet, coordinate, reference in rebases:
                        workbook[sheet][coordinate] = reference
                    if report_tables:
                        _write_report_sheets(workbook, report_tables)
                    workbook.save(staged)
                finally:
                    workbook.close()
            for target, (origin, digest) in evidence.items():
                if target.exists():
                    if _file_digest(target) != digest:
                        raise FileExistsError(f'碎片证据文件在保存期间发生冲突：{target}')
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with origin.open('rb') as reader, target.open('xb') as writer:
                    created.append(target)
                    shutil.copyfileobj(reader, writer)
                shutil.copystat(origin, target)
            if expected_stamp is not None and workbook_stamp(source) != expected_stamp:
                raise ValueError('来源工作簿已更新，请重新打开 ECN 复核与导出。')
            staged.replace(destination)
    except BaseException:
        for target in created:
            target.unlink(missing_ok=True)
        raise
    return destination
