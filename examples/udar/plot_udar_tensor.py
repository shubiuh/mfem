#!/usr/bin/env python3
"""Plot GeoSphere-style UDAR tensor amplitude and phase curve panels."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


COMPONENTS = [f"H{receiver}{transmitter}" for receiver in "xyz" for transmitter in "xyz"]


def unwrap_degrees(values):
    if not values:
        return values
    result = [values[0]]
    for value in values[1:]:
        while value - result[-1] > 180.0:
            value -= 360.0
        while value - result[-1] < -180.0:
            value += 360.0
        result.append(value)
    return result


def load_rows(filename):
    with open(filename, newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{filename} contains no response rows")
    return rows


def curve_groups(rows):
    groups = defaultdict(list)
    for row in rows:
        key = (float(row["frequency_hz"]), float(row["spacing_m"]))
        groups[key].append(row)
    for rows_in_group in groups.values():
        rows_in_group.sort(key=lambda row: float(row["md_m"]))
    return dict(sorted(groups.items()))


def plot_quantity(groups, quantity, output):
    figure, axes = plt.subplots(3, 3, figsize=(14, 11), sharex=True)
    for component_index, component in enumerate(COMPONENTS):
        axis = axes.flat[component_index]
        for (frequency, spacing), rows in groups.items():
            md = [float(row["md_m"]) for row in rows]
            values = []
            for row in rows:
                real = float(row[f"{component}_real"])
                imag = float(row[f"{component}_imag"])
                if quantity == "amplitude":
                    values.append(20.0 * math.log10(max(math.hypot(real, imag), 1e-300)))
                else:
                    values.append(math.degrees(math.atan2(imag, real)))
            if quantity == "phase":
                values = unwrap_degrees(values)
            label = f"{frequency:g} Hz / {spacing:g} m"
            axis.plot(md, values, linewidth=1.5, label=label)
        axis.set_title(component)
        axis.grid(True, alpha=0.3)
        if component_index // 3 == 2:
            axis.set_xlabel("Measured depth [m]")
        if component_index % 3 == 0:
            ylabel = "Amplitude [dB re 1 m$^{-3}$]" if quantity == "amplitude" else "Phase [deg]"
            axis.set_ylabel(ylabel)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", ncol=min(4, len(labels)))
    title = "UDAR full 3x3 magnetic transfer tensor — " + quantity
    figure.suptitle(title, y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(output, dpi=180, bbox_inches="tight")
    print(f"Wrote {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="?", default="udar_tensor.csv")
    parser.add_argument("--prefix", default="udar_tensor")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    groups = curve_groups(load_rows(args.csv))
    prefix = Path(args.prefix)
    plot_quantity(groups, "amplitude", prefix.with_name(prefix.name + "_amplitude.png"))
    plot_quantity(groups, "phase", prefix.with_name(prefix.name + "_phase.png"))
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
