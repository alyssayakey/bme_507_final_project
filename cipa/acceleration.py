"""Optional JIT: retain runnable equations when native Numba loading is blocked."""
import warnings

try:
    from numba import njit
except (ImportError, OSError) as exc:
    warnings.warn(f"Numba unavailable; using slower Python equations: {exc}", RuntimeWarning, stacklevel=2)

    def njit(*args, **kwargs):
        def decorate(function):
            return function
        return decorate
