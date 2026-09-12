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
                               "optimizer", "synthetic_truth", "backend",
                               "incremental_update"}
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
    update = config.get("incremental_update")
    if update is not None:
        if not isinstance(update, dict):
            raise ValueError("incremental_update must be an object")
        allowed = {"previous_result", "active_stations", "freeze_parameters",
                   "initialize_from_previous"}
        if update.keys() - allowed:
            raise ValueError("Unknown incremental_update fields")
        if not isinstance(update.get("previous_result"), str):
            raise ValueError("incremental_update.previous_result must be a path")
        frozen = update.get("freeze_parameters", [])
        if not isinstance(frozen, list) or len(set(frozen)) != len(frozen) or \
                set(frozen) - set(names):
            raise ValueError("freeze_parameters must be unique configured parameters")
        if "active_stations" in update:
            if not isinstance(update["active_stations"], list) or not update["active_stations"]:
                raise ValueError("active_stations must be a nonempty list")
            for station in update["active_stations"]:
                if not isinstance(station, int) or isinstance(station, bool) or station < 0:
                    raise ValueError("active_stations must contain nonnegative integers")
        if "initialize_from_previous" in update and \
                not isinstance(update["initialize_from_previous"], bool):
            raise ValueError("initialize_from_previous must be true or false")
    return config, path.parent


def _write_json(path, record):
    with Path(path).open("w", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _write_iteration_models(path, forward, names):
    """Map each evaluated model to its one or more ParaView mesh outputs."""
    trace = forward.work_dir / "evaluations.jsonl"
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["evaluation", *names, "frequency_mesh_group",
                        "frequencies_hz", "vtk_directory"])
        writer.writeheader()
        if not trace.exists():
            return
        for line in trace.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            contexts = record.get("frequency_mesh_contexts")
            if contexts is None:
                contexts = [{"name": "group_00", "frequencies_hz": [],
                             "vtk": record["vtk"]}]
            for context in contexts:
                row = {"evaluation": record["evaluation"],
                       "frequency_mesh_group": context["name"],
                       "frequencies_hz": ",".join(
                           map(str, context.get("frequencies_hz", []))),
                       "vtk_directory": context["vtk"]}
                row.update(record["conductivity"])
                writer.writerow(row)


def _prepare_incremental_update(config, base, observations):
    """Load a prior result, freeze trusted values, and select new stations.

    The current executable exposes scalar formation parameters, so a frozen
    parameter represents a trusted part of the present parameterization.  The
    output records this explicitly rather than claiming a spatial cell mask.
    """
    update = config.get("incremental_update")
    parameters = [dict(parameter) for parameter in config["parameters"]]
    forward = dict(config["forward"])
    if update is None:
        return observations, parameters, forward, {}, {}

    previous_path = (base / update["previous_result"]).resolve()
    try:
        with previous_path.open(encoding="utf-8") as stream:
            previous = json.load(stream)
    except OSError as exc:
        raise FileNotFoundError(
            f"Cannot read incremental prior result: {previous_path}") from exc
    conductivity = previous.get("conductivity_s_per_m")
    if not isinstance(conductivity, dict):
        raise ValueError("previous_result lacks conductivity_s_per_m")

    previous_values = {}
    for parameter in parameters:
        name = parameter["name"]
        value = conductivity.get(name)
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"previous_result has no valid value for {name}")
        previous_values[name] = float(value)

    frozen_names = update.get("freeze_parameters", [])
    frozen = {name: previous_values[name] for name in frozen_names}
    forward.update(frozen)
    if update.get("initialize_from_previous", True):
        for parameter in parameters:
            if parameter["name"] not in frozen:
                parameter["initial"] = previous_values[parameter["name"]]
                if "prior" not in parameter:
                    parameter["prior"] = parameter["initial"]
                if not parameter["lower"] <= parameter["initial"] <= parameter["upper"]:
                    raise ValueError(
                        f"previous_result value for {parameter['name']} is outside its bounds")
    active = update.get("active_stations")
    if active is not None:
        observations = observations.select_stations(active)
        forward["station_indices"] = active
    parameters = [parameter for parameter in parameters
                  if parameter["name"] not in frozen]
    if not parameters:
        raise ValueError("incremental_update freezes every inversion parameter")
    metadata = {
        "previous_result": str(previous_path),
        "active_stations": active,
        "frozen_conductivity_s_per_m": frozen,
        "initialize_from_previous": update.get("initialize_from_previous", True),
    }
    return observations, parameters, forward, frozen, metadata


def run_inversion(config_path, observed_path, output_dir, *, library=None,
                  backend=None):
    config, base = load_config(config_path)
    observations = TensorData.read(observed_path)
    observations, parameters, forward_config, frozen, update_metadata = \
        _prepare_incremental_update(config, base, observations)
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
    forward = ForwardModel(forward_config, output / "forward", base)
    names = [p["name"] for p in parameters]

    def predict(log_model):
        physical = dict(frozen)
        physical.update(zip(names, map(math.exp, log_model)))
        return forward.predict(physical).select_keys(observations).vector(
            channels, reference=observations)

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
        recovered = dict(frozen)
        recovered.update(zip(names, map(math.exp, result["log_conductivity"])))
        predicted = forward.predict(recovered).select_keys(observations)
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
        result["forward_contexts"] = forward.forward_contexts
        if update_metadata:
            result["incremental_update"] = update_metadata
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
