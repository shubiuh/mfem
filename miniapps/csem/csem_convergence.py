#!/usr/bin/env python3
"""Run paper Scenario 2 h- and p-convergence regressions for CSEM."""

import argparse
import csv
import math
import re
import subprocess
import time
from pathlib import Path


PAPER_AXIS = [
    -20000.0, -10000.0, -5000.0, -2500.0, -1000.0, -500.0, -100.0,
    0.0, 100.0, 500.0, 1000.0, 2500.0, 5000.0, 10000.0, 20000.0,
]
PAPER_Z = [
    -20000.0, -10000.0, -5000.0, -2500.0, -1000.0, -500.0, -100.0,
    0.0, 250.0, 500.0, 2000.0 / 3.0, 2500.0 / 3.0, 1000.0,
    5000.0, 20000.0,
]
QUICK_AXIS = [
    -20000.0, -5000.0, -1000.0, -100.0, 0.0,
    100.0, 1000.0, 5000.0, 20000.0,
]
QUICK_Z = [
    -20000.0, -5000.0, -1000.0, -100.0, 0.0,
    500.0, 1000.0, 5000.0, 20000.0,
]
COMPONENTS = {
    "Ex": ("E", 0.0, 0.0),
    "Hy": ("H", 90.0, 0.0),
    "Hz": ("H", 0.0, 90.0),
}
# Table 9 prints 70,466 for p=2. A conforming tensor-product 14^3 mesh has
# 3*(14*p)*(14*p+1)^2 = 70,644 unknowns; the other printed counts follow that
# formula exactly. Accept the printed value as well in case the supplied mesh
# is replaced by the authors' original mesh.
PAPER_DOFS = {
    1: (9450,),
    2: (70466, 70644),
    3: (232974,),
    4: (545832,),
    5: (1058610,),
}


def refined_axis(axis, levels):
    result = list(axis)
    for _ in range(levels):
        result = [
            value
            for left, right in zip(result[:-1], result[1:])
            for value in (left, 0.5 * (left + right))
        ] + [result[-1]]
    return result


def vertex_id(i, j, k, nx, ny):
    return (k * ny + j) * nx + i


def write_mesh(path, x, y, z):
    nx, ny, nz = len(x), len(y), len(z)
    elements = []
    for k in range(nz - 1):
        if z[k + 1] <= 0.0:
            attribute = 1
        elif z[k] < 500.0:
            attribute = 2
        elif z[k] < 1000.0:
            attribute = 3
        else:
            attribute = 4
        for j in range(ny - 1):
            for i in range(nx - 1):
                vertices = [
                    vertex_id(i, j, k, nx, ny),
                    vertex_id(i + 1, j, k, nx, ny),
                    vertex_id(i + 1, j + 1, k, nx, ny),
                    vertex_id(i, j + 1, k, nx, ny),
                    vertex_id(i, j, k + 1, nx, ny),
                    vertex_id(i + 1, j, k + 1, nx, ny),
                    vertex_id(i + 1, j + 1, k + 1, nx, ny),
                    vertex_id(i, j + 1, k + 1, nx, ny),
                ]
                elements.append((attribute, vertices))

    boundary = []
    for k in range(nz - 1):
        for j in range(ny - 1):
            boundary.append((1, [vertex_id(0, j + 1, k, nx, ny),
                                 vertex_id(0, j, k, nx, ny),
                                 vertex_id(0, j, k + 1, nx, ny),
                                 vertex_id(0, j + 1, k + 1, nx, ny)]))
            boundary.append((2, [vertex_id(nx - 1, j, k, nx, ny),
                                 vertex_id(nx - 1, j + 1, k, nx, ny),
                                 vertex_id(nx - 1, j + 1, k + 1, nx, ny),
                                 vertex_id(nx - 1, j, k + 1, nx, ny)]))
    for k in range(nz - 1):
        for i in range(nx - 1):
            boundary.append((3, [vertex_id(i, 0, k, nx, ny),
                                 vertex_id(i + 1, 0, k, nx, ny),
                                 vertex_id(i + 1, 0, k + 1, nx, ny),
                                 vertex_id(i, 0, k + 1, nx, ny)]))
            boundary.append((4, [vertex_id(i + 1, ny - 1, k, nx, ny),
                                 vertex_id(i, ny - 1, k, nx, ny),
                                 vertex_id(i, ny - 1, k + 1, nx, ny),
                                 vertex_id(i + 1, ny - 1, k + 1, nx, ny)]))
    for j in range(ny - 1):
        for i in range(nx - 1):
            boundary.append((5, [vertex_id(i, j + 1, 0, nx, ny),
                                 vertex_id(i + 1, j + 1, 0, nx, ny),
                                 vertex_id(i + 1, j, 0, nx, ny),
                                 vertex_id(i, j, 0, nx, ny)]))
            boundary.append((6, [vertex_id(i, j, nz - 1, nx, ny),
                                 vertex_id(i + 1, j, nz - 1, nx, ny),
                                 vertex_id(i + 1, j + 1, nz - 1, nx, ny),
                                 vertex_id(i, j + 1, nz - 1, nx, ny)]))

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("MFEM mesh v1.0\n\ndimension\n3\n\nelements\n")
        stream.write(f"{len(elements)}\n")
        for attribute, vertices in elements:
            stream.write(f"{attribute} 5 {' '.join(map(str, vertices))}\n")
        stream.write(f"\nboundary\n{len(boundary)}\n")
        for attribute, vertices in boundary:
            stream.write(f"{attribute} 3 {' '.join(map(str, vertices))}\n")
        stream.write(f"\nvertices\n{nx * ny * nz}\n3\n")
        for zk in z:
            for yj in y:
                for xi in x:
                    stream.write(f"{xi:.16g} {yj:.16g} {zk:.16g}\n")
    return len(elements)


def read_reference(path):
    with path.open(encoding="utf-8") as stream:
        rows = [line for line in stream if not line.startswith("#")]
    return list(csv.DictReader(rows))


def write_inputs(work, component, reference):
    materials = work / "scenario2_materials.csv"
    materials.write_text(
        "attribute,Rh_ohm_m,Rv_ohm_m,dip_degrees,relative_permittivity,invertible\n"
        "1,100000000,100000000,0,1,0\n"
        "2,10000,10000,0,1,1\n"
        "3,100,100,0,1,1\n"
        "4,10000,10000,0,1,1\n",
        encoding="utf-8",
    )
    _, azimuth, dip = COMPONENTS[component]
    receivers = work / f"scenario2_receivers_{component.lower()}.csv"
    with receivers.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("x_m,y_m,z_m,azimuth_degrees,dip_degrees\n")
        for row in reference:
            stream.write(f"0,{row['offset_m']},0,{azimuth},{dip}\n")
    source = work / "scenario2_source.csv"
    source.write_text("-100,0,0\n0,0,0\n100,0,0\n", encoding="utf-8")
    return materials, receivers, source


def parse_response(path):
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def error_metrics(response, reference, component):
    if len(response) != len(reference):
        raise RuntimeError("Computed and reference receiver counts differ")
    amplitude_errors = []
    real_errors = []
    imaginary_errors = []
    for computed, expected in zip(response, reference):
        value = complex(float(computed["real"]), float(computed["imag"]))
        target = complex(float(expected[f"{component}_real"]),
                         float(expected[f"{component}_imag"]))
        amplitude_errors.append(abs(abs(value) - abs(target)) / abs(target))
        real_errors.append(abs(value.real - target.real) / abs(target.real))
        imaginary_errors.append(
            abs(value.imag - target.imag) / abs(target.imag))

    def summarize(errors):
        mean = sum(errors) / len(errors)
        rms = math.sqrt(sum(error * error for error in errors) / len(errors))
        return 100.0 * mean, 100.0 * rms, 100.0 * max(errors)

    return summarize(amplitude_errors), summarize(real_errors), summarize(
        imaginary_errors)


def run_case(args, work, reference, study, level, order, axes, component):
    mesh = work / f"scenario2_{study}_h{level}.mesh"
    elements = write_mesh(mesh, *axes)
    materials, receivers, source = write_inputs(work, component, reference)
    prefix = work / f"scenario2_{study}_h{level}_p{order}_{component.lower()}"
    field, _, _ = COMPONENTS[component]
    command = [
        str(args.executable), "-m", str(mesh), "-mat", str(materials),
        "-mp", "region", "-o", str(order), "-f", "1000",
        "-s", args.solver, "-bc", "pec", "-rf", str(receivers),
        "-rm", "point", "-rc", field, "-sm", "line", "-sqo", "1",
        "-sf", str(source), "-si", "0.5", "-sr", "1", "-no-vis",
        "-p", str(prefix),
    ]
    start = time.perf_counter()
    completed = subprocess.run(command, cwd=args.executable.parent.parent.parent,
                               text=True, capture_output=True)
    elapsed = time.perf_counter() - start
    log_file = Path(f"{prefix}.log")
    log_file.write_text(completed.stdout + completed.stderr, encoding="utf-8")
    if completed.returncode:
        print(completed.stdout)
        print(completed.stderr)
        raise RuntimeError(
            f"CSEM failed with return code {completed.returncode} for "
            f"{study} h={level} p={order} {component}; see {log_file}")
    match = re.search(r"CSEM mesh: .*?, ([0-9]+) true H\(curl\) DOFs",
                      completed.stdout)
    if not match:
        raise RuntimeError("Could not read the CSEM degree-of-freedom count")
    response_file = Path(f"{prefix}_response.csv")
    response = parse_response(response_file)
    amplitude, real_part, imaginary_part = error_metrics(
        response, reference, component)
    mean, rms, maximum = amplitude
    print(f"{study:2s} h={level} p={order} {component}: "
          f"amplitude mean={mean:.6g}% rms={rms:.6g}% max={maximum:.6g}%")
    return {
        "study": study,
        "h_level": level,
        "order": order,
        "elements": elements,
        "dofs": int(match.group(1)),
        "component": component,
        "mean_amplitude_error_percent": mean,
        "rms_amplitude_error_percent": rms,
        "max_amplitude_error_percent": maximum,
        "mean_real_error_percent": real_part[0],
        "rms_real_error_percent": real_part[1],
        "max_real_error_percent": real_part[2],
        "mean_imaginary_error_percent": imaginary_part[0],
        "rms_imaginary_error_percent": imaginary_part[1],
        "max_imaginary_error_percent": imaginary_part[2],
        "elapsed_seconds": elapsed,
        "response_file": str(response_file),
        "log_file": str(log_file),
    }


def parse_ints(text):
    return [int(value) for value in text.replace(",", " ").split()]


def main():
    source = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--reference", type=Path,
                        default=source / "paper_scenario2_reference.csv")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("csem-convergence-results"))
    parser.add_argument("--suite", choices=("quick", "paper"),
                        default="quick")
    parser.add_argument("--orders", default="")
    parser.add_argument("--h-levels", default="")
    parser.add_argument("--h-order", type=int, default=2)
    parser.add_argument("--components", default="Ex,Hy,Hz")
    parser.add_argument("--solver", default="mumps")
    parser.add_argument("--p-tolerance-percent", type=float, default=0.0,
                        help="maximum final p-study mean error; 0 disables")
    parser.add_argument("--h-tolerance-percent", type=float, default=0.0,
                        help="maximum final h-study mean error; 0 disables")
    parser.add_argument("--minimum-reduction-fraction", type=float,
                        default=0.01)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    args.executable = args.executable.resolve()
    args.reference = args.reference.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not args.executable.is_file():
        raise FileNotFoundError(args.executable)
    if not 0.0 <= args.minimum_reduction_fraction < 1.0:
        raise ValueError("Minimum reduction fraction must lie in [0, 1)")

    orders = parse_ints(args.orders) if args.orders else (
        [1, 2, 3] if args.suite == "quick" else [1, 2, 3, 4, 5])
    h_levels = parse_ints(args.h_levels) if args.h_levels else [0, 1]
    components = [item.strip() for item in args.components.split(",")]
    if any(item not in COMPONENTS for item in components):
        raise ValueError("Components must be selected from Ex, Hy, and Hz")
    base_axis = QUICK_AXIS if args.suite == "quick" else PAPER_AXIS
    base_z = QUICK_Z if args.suite == "quick" else PAPER_Z
    reference = read_reference(args.reference)

    rows = []
    for component in components:
        axes = (base_axis, base_axis, base_z)
        for order in orders:
            rows.append(run_case(args, args.output_dir, reference, "p", 0,
                                 order, axes, component))
        for level in h_levels:
            if level == 0 and args.h_order in orders:
                continue
            axes = tuple(refined_axis(axis, level)
                         for axis in (base_axis, base_axis, base_z))
            rows.append(run_case(args, args.output_dir, reference, "h", level,
                                 args.h_order, axes, component))

    if args.suite == "paper":
        for row in rows:
            if row["study"] != "p":
                continue
            expected_dofs = PAPER_DOFS.get(row["order"])
            if expected_dofs is not None and row["dofs"] not in expected_dofs:
                raise RuntimeError(
                    f"Paper p={row['order']} DOF mismatch: "
                    f"expected one of {expected_dofs}, got {row['dofs']}")

    summary = args.output_dir / "scenario2_convergence.csv"
    with summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    if not args.report_only:
        for component in components:
            p_rows = sorted((row for row in rows if row["study"] == "p" and
                             row["component"] == component),
                            key=lambda row: row["order"])
            h_rows = sorted([row for row in rows if row["study"] == "h" and
                             row["component"] == component] +
                            [row for row in p_rows
                             if row["order"] == args.h_order],
                            key=lambda row: row["h_level"])
            if (args.p_tolerance_percent > 0.0 and len(p_rows) > 1 and
                    p_rows[-1]["mean_amplitude_error_percent"] >
                    args.p_tolerance_percent):
                raise RuntimeError(
                    f"{component} final p-error exceeds tolerance")
            if (args.h_tolerance_percent > 0.0 and len(h_rows) > 1 and
                    h_rows[-1]["mean_amplitude_error_percent"] >
                    args.h_tolerance_percent):
                raise RuntimeError(
                    f"{component} final h-error exceeds tolerance")
            reduction = 1.0 - args.minimum_reduction_fraction
            best_p = min(row["mean_amplitude_error_percent"] for row in p_rows)
            if len(p_rows) > 1 and best_p > reduction * p_rows[0]["mean_amplitude_error_percent"]:
                raise RuntimeError(
                    f"{component} p-refinement did not reduce error enough")
            if (len(h_rows) > 1 and
                    h_rows[-1]["mean_amplitude_error_percent"] >
                    reduction * h_rows[0]["mean_amplitude_error_percent"]):
                raise RuntimeError(
                    f"{component} h-refinement did not reduce error enough")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
