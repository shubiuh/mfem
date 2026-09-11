"""Configuration, synthetic data generation, and auditable inversion outputs."""
from __future__ import annotations

import json
import math
import csv
from pathlib import Path

from .core import minimize
from .forward import CHANNELS, PARAMETERS, ForwardModel, TensorData


def load_config(path):
    path = Path(path).resolve()
    with path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    unknown = config.keys() - {"forward", "parameters", "channels", "noise",
                               "optimizer", "synthetic_truth", "backend"}
    if unknown:
        raise ValueError(f"Unknown configuration keys: {sorted(unknown)}")
    if config.get("backend", "rol") not in {"rol", "tao"}:
        raise ValueError("backend must be 'rol' or 'tao'")
    parameters = config.get("parameters", [])
    if not parameters:
        raise ValueError("At least one inversion parameter is required")
    names = [p["name"] for p in parameters]
    if len(set(names)) != len(names) or set(names) - PARAMETERS.keys():
        raise ValueError("Use unique sigma_background and/or sigma_target parameters")
    for p in parameters:
        if p.keys() - {"name", "initial", "lower", "upper", "prior"}:
            raise ValueError(f"Unknown parameter fields for {p['name']}")
        for field in ("initial", "lower", "upper", "prior"):
            value = p.get(field, p["initial"])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{p['name']}.{field} must be positive and finite")
        if not p["lower"] < p["upper"] or not p["lower"] <= p["initial"] <= p["upper"]:
            raise ValueError(f"Invalid conductivity bounds for {p['name']}")
    return config, path.parent


def _write_json(path, record):
    with Path(path).open("w", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _write_iteration_models(path, forward, names):
    """Write the model attached to every forward evaluation for ParaView review."""
    trace = forward.work_dir / "evaluations.jsonl"
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["evaluation", *names, "vtk_directory"])
        writer.writeheader()
        if not trace.exists():
            return
        for line in trace.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            row = {"evaluation": record["evaluation"],
                   "vtk_directory": record["vtk"]}
            row.update(record["conductivity"])
            writer.writerow(row)


def run_inversion(config_path, observed_path, output_dir, *, library=None,
                  backend=None):
    config, base = load_config(config_path)
    observations = TensorData.read(observed_path)
    channels = config.get("channels", CHANNELS)
    noise = config.get("noise", {"relative": 0.01, "absolute": 1e-8})
    if noise.keys() - {"relative", "absolute"}:
        raise ValueError("Noise fields are relative and absolute")
    observed = observations.vector(channels)
    stddev = observations.standard_deviations(channels, noise["relative"], noise["absolute"])
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    _write_json(output / "configuration.json", config)
    observations.write(output / "observed.csv")
    forward = ForwardModel(config["forward"], output / "forward", base)
    parameters = config["parameters"]
    names = [p["name"] for p in parameters]

    def predict(log_model):
        physical = dict(zip(names, map(math.exp, log_model)))
        return forward.predict(physical).vector(channels, reference=observations)

    initial = [math.log(p["initial"]) for p in parameters]
    lower = [math.log(p["lower"]) for p in parameters]
    upper = [math.log(p["upper"]) for p in parameters]
    prior = [math.log(p.get("prior", p["initial"])) for p in parameters]
    try:
        selected_backend = (backend or config.get("backend", "rol")).lower()
        result = minimize(predict, initial, lower, upper, observed, stddev, prior=prior,
                          library=library, backend=selected_backend,
                          log_path=output / f"{selected_backend}.log",
                          **config.get("optimizer", {}))
        recovered = dict(zip(names, map(math.exp, result["log_conductivity"])))
        predicted = forward.predict(recovered)
        predicted.write(output / "predicted.csv")
        residual = [(p-o)/s for p, o, s in
                    zip(predicted.vector(channels, observations), observed, stddev)]
        result.update(conductivity_s_per_m=recovered,
                      resistivity_ohm_m={name: 1/v for name, v in recovered.items()},
                      data_count=len(observed), channels=channels,
                      weighted_rms=math.sqrt(sum(r*r for r in residual)/len(residual)),
                      forward_server_starts=forward.server_starts,
                      config_source=str(Path(config_path).resolve()),
                      observation_source=str(Path(observed_path).resolve()),
                      parameter_order=names,
                      forward_executable=str(forward.executable), mesh=str(forward.mesh))
        _write_iteration_models(output / "iteration_models.csv", forward, names)
        _write_json(output / "result.json", result)
        _write_json(output / "residuals.json", {"weighted_residuals": residual,
                                               "standard_deviations": stddev,
                                               "ordering": "sorted station/frequency/spacing, then channels"})
        forward.close()
        return result
    except BaseException as exc:
        forward.close(force=True)
        _write_json(output / "failure.json", {"type": type(exc).__name__, "message": str(exc)})
        raise


def synthetic_recovery(config_path, output_dir, *, library=None, backend=None):
    """Noiseless same-discretization recovery: a software check, not physics validation."""
    config, base = load_config(config_path)
    truth = config.get("synthetic_truth")
    names = [p["name"] for p in config["parameters"]]
    if not isinstance(truth, dict) or set(truth) != set(names):
        raise ValueError("synthetic_truth must specify exactly the inverted parameters")
    for p in config["parameters"]:
        if not p["lower"] <= truth[p["name"]] <= p["upper"]:
            raise ValueError("Synthetic truth must lie within the inversion bounds")
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    forward = ForwardModel(config["forward"], output / "true_model", base)
    try:
        forward.predict(truth).write(output / "observed.csv")
    finally:
        forward.close(force=True)
    result = run_inversion(config_path, output / "observed.csv",
                           output / "inversion", library=library,
                           backend=backend)
    errors = {name: abs(result["conductivity_s_per_m"][name]/truth[name]-1) for name in names}
    result.update(synthetic_truth=truth, relative_parameter_error=errors,
                  recovery_passed=result["converged"] and max(errors.values()) < 0.01,
                  validation_scope="Noiseless recovery using the same mesh and forward solver")
    _write_json(output / "recovery.json", result)
    return result
