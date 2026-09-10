"""Small C ABI wrapper; optimization and differentiation execute in C++."""
from __future__ import annotations

import ctypes as ct
import math
import os
from pathlib import Path


class _Options(ct.Structure):
    _fields_ = [("max_iterations", ct.c_int),
                ("gradient_tolerance", ct.c_double),
                ("step_tolerance", ct.c_double),
                ("difference_step", ct.c_double),
                ("regularization", ct.c_double)]


class _Result(ct.Structure):
    _fields_ = [("status", ct.c_int), ("iterations", ct.c_int),
                ("forward_evaluations", ct.c_int),
                ("initial_objective", ct.c_double),
                ("final_objective", ct.c_double),
                ("projected_gradient_norm", ct.c_double)]


_DoublePtr = ct.POINTER(ct.c_double)
_Predict = ct.CFUNCTYPE(ct.c_int, ct.c_int, _DoublePtr, ct.c_int,
                      _DoublePtr, ct.c_void_p)
STATUS = {0: "converged", 1: "iteration_limit", 2: "step_tolerance",
          3: "nonfinite_step_or_gradient", 4: "user_stopped", 5: "unknown"}


def default_library() -> Path:
    if os.environ.get("UDAR_ROL_LIBRARY"):
        return Path(os.environ["UDAR_ROL_LIBRARY"]).resolve()
    filename = "udar_rol.dll" if os.name == "nt" else "libudar_rol.so"
    return Path(__file__).resolve().parents[3] / "build" / "udar-inversion" / filename


def minimize(predict, initial, lower, upper, observed, stddev, *, prior=None,
             library=None, log_path="rol.log", max_iterations=60,
             gradient_tolerance=1e-6, step_tolerance=1e-10,
             difference_step=1e-3, regularization=0.0):
    """Minimize 0.5*||((predict(m)-observed)/stddev)||^2 + beta/2*||m-prior||^2.

    All parameter vectors are natural-log conductivities. predict(m) returns
    real scalar channels in observation order. Calls are synchronous; callback
    exceptions are propagated after C++ unwinds. A completed call is not
    necessarily converged: inspect the returned status/converged fields.
    """
    initial, lower, upper = list(initial), list(lower), list(upper)
    observed, stddev = list(observed), list(stddev)
    prior = list(initial if prior is None else prior)
    n, nd = len(initial), len(observed)
    if n == 0 or nd == 0 or any(len(v) != n for v in (lower, upper, prior)):
        raise ValueError("Nonempty, matching parameter dimensions required")
    if len(stddev) != nd:
        raise ValueError("Standard deviations must match observations")
    if not isinstance(max_iterations, int) or not 0 < max_iterations < 2**31:
        raise ValueError("max_iterations must be a positive 32-bit integer")
    # Validate before converting to fixed-size buffers, and before loading DLLs.
    for v in (initial, lower, upper, prior, observed, stddev):
        if not all(math.isfinite(float(x)) for x in v):
            raise ValueError("All vectors must contain finite numbers")
    if any(not lo < hi or not lo <= x <= hi
           for x, lo, hi in zip(initial, lower, upper)):
        raise ValueError("Initial parameters must lie inside strict lower/upper bounds")
    if any(x <= 0 for x in stddev):
        raise ValueError("Standard deviations must be positive")
    for name, value in (("gradient_tolerance", gradient_tolerance),
                        ("step_tolerance", step_tolerance),
                        ("difference_step", difference_step)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(regularization) or regularization < 0:
        raise ValueError("regularization must be finite and nonnegative")

    dll_path = Path(library or default_library()).resolve()
    # Python 3.8+ uses explicit DLL search directories on Windows. Keep these
    # handles alive through the call, including any late-loaded oneMKL DLLs.
    dll_dirs = []
    if os.name == "nt":
        directories = [str(dll_path.parent)]
        if os.environ.get("MKLROOT"):
            directories.append(str(Path(os.environ["MKLROOT"]) / "bin"))
        directories += os.environ.get("UDAR_DLL_DIRS", "").split(os.pathsep)
        for directory in dict.fromkeys(directories):
            if directory and Path(directory).is_dir():
                dll_dirs.append(os.add_dll_directory(directory))
    try:
        lib = ct.CDLL(str(dll_path))
        lib.udar_abi_version.restype = ct.c_int
        if lib.udar_abi_version() != 1:
            raise RuntimeError("Incompatible UDAR C ABI version; rebuild the library")
        lib.udar_minimize.argtypes = [ct.c_int, ct.c_int] + [_DoublePtr]*6 + [
            ct.POINTER(_Options), _Predict, ct.c_void_p, ct.c_char_p,
            ct.POINTER(_Result), ct.POINTER(ct.c_char), ct.c_int]
        lib.udar_minimize.restype = ct.c_int
        model = (ct.c_double*n)(*initial)
        buffers = [(ct.c_double*len(v))(*v)
                   for v in (lower, upper, prior, observed, stddev)]
        failure = []

        @_Predict
        def callback(count, x, count_data, output, _):
            try:
                prediction = list(predict([x[j] for j in range(count)]))
                if len(prediction) != count_data:
                    raise ValueError("Forward callback returned the wrong data size")
                for i, v in enumerate(prediction):
                    if not math.isfinite(v):
                        raise ValueError("Forward callback returned non-finite data")
                    output[i] = v
                return 0
            except BaseException as exc:
                failure.append(exc)
                return 1

        options = _Options(max_iterations, gradient_tolerance, step_tolerance,
                           difference_step, regularization)
        result, error = _Result(), ct.create_string_buffer(4096)
        rc = lib.udar_minimize(n, nd, model, *buffers, ct.byref(options), callback,
                              None, str(Path(log_path).resolve()).encode("utf-8"),
                              ct.byref(result), error, len(error))
        if failure:
            raise failure[0]
        if rc != 0:
            raise RuntimeError(error.value.decode("utf-8", errors="replace"))
        record = {name: getattr(result, name) for name, _ in _Result._fields_}
        record.update(log_conductivity=list(model),
                      status=STATUS.get(result.status, "unknown"),
                      status_code=result.status, converged=result.status == 0,
                      backend="Trilinos/ROL projected quasi-Newton (L-BFGS)",
                      derivative="C++ finite differences of residuals",
                      library=str(dll_path))
        return record
    finally:
        for handle in dll_dirs:
            handle.close()
