"""Preserve complete MS2 evidence when a JSON cell exceeds Excel's limit."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from openpyxl import load_workbook

from .io_utils import write_excel_atomic


def write_identifications(data,path):
    path=Path(path)
    output=data.copy()
    evidence={}
    external_indices=[]
    if 'ms2_evidence_json' in output:
        for index,value in output['ms2_evidence_json'].items():
            if isinstance(value,str) and len(value)>32767:
                key=hashlib.sha256(value.encode('utf-8')).hexdigest()
                evidence[key]=value
                output.loc[index,'ms2_evidence_json']=json.dumps({'external_evidence_sha256':key},separators=(',',':'))
                output.loc[index,'evidence_sha256']=key
                external_indices.append(index)
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
        mask=output.index.isin(external_indices)
        output.loc[mask,'evidence_file']=sidecar.name
    try:
        write_excel_atomic(output,path,index=False)
    except BaseException:
        if created is not None and created.resolve().parent==path.resolve().parent:
            created.unlink(missing_ok=True)
        raise
    return path


def _file_digest(path: Path) -> bytes:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.digest()


def copy_identification_bundle(source: str | Path, destination: str | Path) -> Path:
    """Save all worksheets and referenced evidence, publishing the XLSX last.

    Safe relative references keep the workbook byte-for-byte intact. Absolute
    references or references outside the destination folder are localized;
    only those reference cells are changed. Existing unequal evidence files
    are never overwritten, and a failed save leaves the previous XLSX intact.
    """
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        return destination
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
            if rebases:
                workbook = load_workbook(staged)
                try:
                    for sheet, coordinate, reference in rebases:
                        workbook[sheet][coordinate] = reference
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
            staged.replace(destination)
    except BaseException:
        for target in created:
            target.unlink(missing_ok=True)
        raise
    return destination
