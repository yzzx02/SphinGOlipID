"""Optional live XCMS integration; pure-Python tests never require R."""
import os
import subprocess

import pytest


@pytest.fixture(scope="session")
def xcms_runtime():
    from sphingolipid_toolkit.untargeted import find_rscript
    try:
        executable=find_rscript()
    except RuntimeError as exc:
        pytest.skip(str(exc))
    code=('p<-c("xcms","MSnbase","mzR","jsonlite","BiocParallel"); '
          'if (!all(vapply(p,requireNamespace,logical(1),quietly=TRUE))) quit(status=2); '
          'cat(as.character(packageVersion("xcms")))')
    options={'creationflags':subprocess.CREATE_NO_WINDOW} if os.name=='nt' else {}
    try:
        result=subprocess.run([str(executable),'--vanilla','-e',code],capture_output=True,
                              text=True,encoding='utf-8',errors='replace',timeout=45,**options)
    except (OSError,subprocess.TimeoutExpired) as exc:
        pytest.skip(f'Live XCMS runtime unavailable: {exc}')
    if result.returncode:
        pytest.skip('Live XCMS packages unavailable')
    return {'executable':executable,'version':result.stdout.strip()}
