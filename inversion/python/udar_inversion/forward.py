"""Validated persistent-process adapter for the MPI-enabled udar3p solver."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import queue
import subprocess
import threading
import time

CHANNELS = [f"H{a}{b}_{part}" for a in "xyz" for b in "xyz"
            for part in ("real", "imag")]
KEYS = ["station", "frequency_hz", "spacing_m"]
GEOMETRY = ["md_m", "tx_x", "tx_y", "tx_z", "rx_x", "rx_y", "rx_z"]
PARAMETERS = {"sigma_background": "-sig-bg", "sigma_target": "-sig-target"}
OPTIONS = {
    "frequencies": "-freq", "spacings": "-sp", "stations": "-ns",
    "survey_start": "-p0", "survey_end": "-p1", "tool_up": "-up",
    "formation_model": "-model", "sigma_background": "-sig-bg",
    "sigma_target": "-sig-target", "slab_center": "-slab-c",
    "slab_normal": "-slab-n", "slab_thickness": "-slab-t",
    "interface_width": "-iw", "box_minimum": "-box-min",
    "box_maximum": "-box-max", "order": "-o",
    "serial_refinements": "-rs", "parallel_refinements": "-rp",
    "boundary_condition": "-bc", "source_type": "-src",
    "loop_radius": "-loop-r", "wire_radius": "-wire-r",
    "relative_permittivity": "-epsr", "relative_permeability": "-mur",
    "relative_tolerance": "-rtol", "maximum_iterations": "-maxit",
    "vtk_directory": "-vtk",
}


class TensorData:
    def __init__(self, rows, fieldnames=None):
        self.rows = rows
        self.fieldnames = fieldnames or KEYS + GEOMETRY + CHANNELS
        if not rows:
            raise ValueError("Tensor data is empty")
        self.by_key = {}
        for row in rows:
            for name in KEYS + GEOMETRY + CHANNELS:
                if name not in row or not math.isfinite(float(row[name])):
                    raise ValueError(f"Missing or non-finite tensor field: {name}")
            key = tuple(float(row[name]) for name in KEYS)
            if key in self.by_key:
                raise ValueError(f"Duplicate station/frequency/spacing: {key}")
            self.by_key[key] = row

    @classmethod
    def read(cls, path):
        with Path(path).open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            return cls(list(reader), reader.fieldnames)

    def write(self, path):
        with Path(path).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=self.fieldnames)
            writer.writeheader()
            writer.writerows(self.rows)

    def vector(self, channels=CHANNELS, reference=None):
        if not channels or len(set(channels)) != len(channels) or any(
                name not in CHANNELS for name in channels):
            raise ValueError("Select unique Hxx_real ... Hzz_imag channels")
        ref = reference or self
        if self.by_key.keys() != ref.by_key.keys():
            raise ValueError("Observed/predicted station, frequency or spacing keys differ")
        result = []
        # Sort physical keys, never rely on file row order.
        for key in sorted(ref.by_key):
            row, expected = self.by_key[key], ref.by_key[key]
            for name in GEOMETRY:
                if not math.isclose(float(row[name]), float(expected[name]),
                                    rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError(f"Observed/predicted geometry differs: {key}, {name}")
            result.extend(float(row[name]) for name in channels)
        return result

    def standard_deviations(self, channels, relative, absolute):
        if not math.isfinite(relative) or relative < 0 or not math.isfinite(absolute) or absolute <= 0:
            raise ValueError("Noise requires relative >= 0 and a positive absolute floor")
        self.vector(channels)
        result = []
        for key in sorted(self.by_key):
            row = self.by_key[key]
            for name in channels:
                base = name.rsplit("_", 1)[0]
                amplitude = math.hypot(float(row[base+"_real"]), float(row[base+"_imag"]))
                result.append(math.hypot(relative*amplitude, absolute))
        return result


class ForwardModel:
    """Replaceable prediction interface with a persistent udar3p process.

    ROL is serial; MPI parallelism belongs to the forward child processes.
    Constructor paths are relative to base_dir, never to temporary run folders.
    """
    def __init__(self, configuration, work_dir, base_dir="."):
        self.config = dict(configuration)
        self.base_dir = Path(base_dir).resolve()
        self.work_dir = Path(work_dir).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=False)
        self.vtk_dir = self.work_dir / "vtk"
        self.vtk_dir.mkdir()
        unknown = self.config.keys() - OPTIONS.keys() - {
            "executable", "mesh", "solver", "launcher", "timeout_seconds",
            "persistent"}
        if unknown:
            raise ValueError(f"Unknown forward options: {sorted(unknown)}")
        self.executable = self._path(self.config["executable"])
        self.mesh = self._path(self.config["mesh"])
        if not self.executable.is_file() or not self.mesh.is_file():
            raise FileNotFoundError("Forward executable and mesh must both exist")
        self.timeout = float(self.config.get("timeout_seconds", 300))
        if not math.isfinite(self.timeout) or self.timeout <= 0:
            raise ValueError("timeout_seconds must be positive and finite")
        self.launcher = self.config.get("launcher", [])
        if not isinstance(self.launcher, list) or not all(isinstance(x, str) for x in self.launcher):
            raise ValueError("launcher must be an argument list, e.g. ['mpiexec','-n','2']")
        self.solver = self.config.get("solver", "mumps")
        if self.solver not in ("mumps", "gmres"):
            raise ValueError("solver must be mumps or gmres")
        if not isinstance(self.config.get("formation_model", 0), int) or self.config.get("formation_model", 0) not in (0, 1, 2):
            raise ValueError("formation_model must be 0, 1 or 2")
        if bool(self.config.get("survey_start")) != bool(self.config.get("survey_end")):
            raise ValueError("Supply both survey_start and survey_end")
        self.persistent = self.config.get("persistent", True)
        if not isinstance(self.persistent, bool):
            raise ValueError("persistent must be true or false")
        self.cache = {}
        self.evaluations = 0
        self.server_starts = 0
        self._process = None
        self._reader = None
        self._lines = None
        self._server_log = None

    def _path(self, path):
        return (self.base_dir / path).resolve()

    def _arguments(self, conductivities=None, server=False, output=None):
        args = list(self.launcher) + [str(self.executable), "-m", str(self.mesh),
                                     "-"+self.solver]
        if server:
            args.extend(["-server", str(self.work_dir)])
        else:
            args.extend(["-out", str(output)])
        settings = {"formation_model": 0, "vtk_directory": str(self.vtk_dir),
                    **self.config,
                    **(conductivities or {})}
        for name, flag in OPTIONS.items():
            if name in settings:
                value = settings[name]
                text = ",".join(map(str, value)) if isinstance(value, list) else str(value)
                args.extend([flag, text])
        return args

    @staticmethod
    def _read_output(process, lines, server_log):
        try:
            for line in process.stdout:
                server_log.write(line)
                server_log.flush()
                lines.put(line.rstrip("\r\n"))
        except (OSError, ValueError):
            # Closing a failed process can close its pipe while this daemon
            # thread is blocked in a read. The process status carries the error.
            pass
        finally:
            lines.put(None)

    def _wait_for(self, marker):
        deadline = time.monotonic() + self.timeout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f"Timed out waiting for {marker}; see {self._server_log.name}")
            try:
                line = self._lines.get(timeout=remaining)
            except queue.Empty as exc:
                raise TimeoutError(f"Timed out waiting for {marker}; see {self._server_log.name}") from exc
            if line is None:
                code = self._process.poll()
                raise RuntimeError(f"Persistent forward server exited with {code}; see {self._server_log.name}")
            if line == marker:
                return
            if line.startswith("UDAR_SERVER_ERROR"):
                raise RuntimeError(f"Persistent forward server rejected a request: {line}")

    def _start_server(self):
        if self._process is not None and self._process.poll() is None:
            return
        self.close(force=True)
        args = self._arguments(server=True)
        self._server_log = (self.work_dir / "server.log").open("a", encoding="utf-8")
        self._lines = queue.Queue()
        try:
            self._process = subprocess.Popen(
                args, cwd=self.base_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                bufsize=1,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        except BaseException:
            self._server_log.close()
            self._server_log = None
            raise
        self.server_starts += 1
        self._reader = threading.Thread(
            target=self._read_output,
            args=(self._process, self._lines, self._server_log), daemon=True)
        self._reader.start()
        try:
            self._wait_for("UDAR_SERVER_READY 1")
        except BaseException:
            self.close(force=True)
            raise

    def _server_prediction(self, conductivities, request):
        background = conductivities.get("sigma_background",
                                        self.config.get("sigma_background", 0.1))
        target = conductivities.get("sigma_target",
                                    self.config.get("sigma_target", background))
        if any(not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0
               for x in (background, target)):
            raise ValueError("Persistent mode currently requires scalar isotropic conductivities")
        last_error = None
        for attempt in range(2):
            try:
                self._start_server()
                self._process.stdin.write(f"E {request} {background:.17g} {target:.17g}\n")
                self._process.stdin.flush()
                self._wait_for(f"UDAR_SERVER_DONE {request}")
                return self.work_dir / f"evaluation_{request:05d}.csv"
            except (OSError, ValueError, RuntimeError, TimeoutError) as exc:
                last_error = exc
                self.close(force=True)
                if attempt == 0:
                    continue
        raise RuntimeError("Persistent forward server failed twice") from last_error

    def close(self, force=False):
        process, self._process = self._process, None
        if process is not None and process.poll() is None:
            if not force:
                try:
                    process.stdin.write("Q\n")
                    process.stdin.flush()
                    process.wait(timeout=5)
                except (BrokenPipeError, OSError, ValueError,
                        subprocess.TimeoutExpired):
                    force = True
            if force and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
        reader, self._reader = self._reader, None
        if reader is not None:
            reader.join(timeout=1)
        for stream in (() if process is None else (process.stdin, process.stdout)):
            if stream:
                try:
                    stream.close()
                except (OSError, ValueError):
                    pass
        if reader is not None and reader.is_alive():
            reader.join(timeout=1)
        if self._server_log is not None:
            self._server_log.close()
            self._server_log = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def __del__(self):
        try:
            self.close(force=True)
        except Exception:
            pass

    def predict(self, conductivities):
        if not conductivities or conductivities.keys() - PARAMETERS.keys():
            raise ValueError("Supported parameters: sigma_background, sigma_target")
        if "sigma_target" in conductivities and self.config.get("formation_model", 0) == 0:
            raise ValueError("A uniform formation has no target conductivity parameter")
        if any(not math.isfinite(v) or v <= 0 for v in conductivities.values()):
            raise ValueError("Conductivities must be positive finite values in S/m")
        key = tuple(sorted(conductivities.items()))
        if key in self.cache:
            return self.cache[key]
        self.evaluations += 1
        stem = self.work_dir / f"evaluation_{self.evaluations:05d}"
        output = stem.with_suffix(".csv")
        start = time.monotonic()
        if self.persistent:
            output = self._server_prediction(conductivities, self.evaluations)
            args = ["persistent-server", f"request={self.evaluations}"]
        else:
            args = self._arguments(conductivities, output=output)
            # Communicator and library runtime are initialized inside the executable.
            # shell=False keeps paths, signs, and arguments literal on both platforms.
            with stem.with_suffix(".log").open("w", encoding="utf-8") as log:
                try:
                    proc = subprocess.run(args, cwd=self.base_dir, stdout=log,
                                          stderr=subprocess.STDOUT, timeout=self.timeout,
                                          creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
                except subprocess.TimeoutExpired as exc:
                    raise RuntimeError(f"Forward solve timed out; see {log.name}") from exc
            if proc.returncode != 0:
                raise RuntimeError(f"Forward solve exited with {proc.returncode}; see {stem.with_suffix('.log')}")
        data = TensorData.read(output)
        record = {"evaluation": self.evaluations, "conductivity": conductivities,
                  "persistent": self.persistent,
                  "seconds": time.monotonic()-start, "command": args,
                  "output": str(output),
                  "vtk": str(self.vtk_dir / f"evaluation_{self.evaluations:05d}")}
        with (self.work_dir / "evaluations.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record, allow_nan=False)+"\n")
        if len(self.cache) >= 256:
            self.cache.clear()
        self.cache[key] = data
        return data
