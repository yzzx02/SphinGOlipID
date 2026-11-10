"""One identification backend; retired interfaces fail explicitly."""
CURRENT_BACKEND = 'lipidgate'

class RemovedBackendError(ValueError):
    """A retired identification route or parameter was requested."""

def removed_backend(*args, **kwargs):
    raise RemovedBackendError(
        '旧鉴定后端已移除；请使用当前鞘脂分池算法与 min_total_score、'
        'precursor_ppm、fragment_ppm/fragment_da、top_n 参数。')

def require_current_backend(value):
    if value != CURRENT_BACKEND:
        removed_backend()

def reject_retired_parameters(config):
    for name in ('min_fragment_intensity','min_matched_fragments','min_match_score'):
        if getattr(config,name,None) is not None:
            raise RemovedBackendError(f'{name} 是旧后端参数，已停用；请使用当前分池门控和 min_total_score。')
