"""Retired notebook backend; kept only to explain failures to old callers."""
from .backend_policy import removed_backend

def __getattr__(name):
    if name.startswith('__'):raise AttributeError(name)
    removed_backend()
