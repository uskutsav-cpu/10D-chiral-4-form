"""Compatibility entrypoint: all exact degree-12 source targets use hardened v2."""
from .degree12_exact_targets import *  # noqa: F401,F403
from .degree12_exact_targets import recover_exact_targets

recover_missing_targets = recover_exact_targets
