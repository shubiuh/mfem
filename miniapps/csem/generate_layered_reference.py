#!/usr/bin/env python3
"""Regenerate the paper Scenario 2 layered-earth reference response.

This optional maintainer utility requires empymod 2.6.0 and NumPy. Routine
regression runs consume the committed CSV and do not require either package.
"""

import argparse
import csv
from pathlib import Path

import empymod
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("paper_scenario2_reference.csv"),
    )
    args = parser.parse_args()
    if empymod.__version__ != "2.6.0":
        raise RuntimeError("Reference generation is pinned to empymod 2.6.0")

    offsets = np.arange(1000.0, 5000.1, 250.0)
    source = [-100.0, 100.0, 0.0, 0.0, 0.0, 0.0]
    common = {
        "src": source,
        "depth": [0.0, 500.0, 1000.0],
        "res": [1.0e8, 1.0e4, 100.0, 1.0e4],
        "freqtime": 1000.0,
        "epermH": [1.0] * 4,
        "epermV": [1.0] * 4,
        "mpermH": [1.0] * 4,
        "mpermV": [1.0] * 4,
        "srcpts": 31,
        "strength": 0.5,
        "ht": "qwe",
        "htarg": {"rtol": 1.0e-10, "atol": 1.0e-30,
                  "nquad": 81, "maxint": 80},
        "verb": 0,
    }
    receivers = {
        "Ex": ([np.zeros(offsets.size), offsets, 0.0, 0.0, 0.0], False),
        "Hy": ([np.zeros(offsets.size), offsets, 0.0, 90.0, 0.0], True),
        "Hz": ([np.zeros(offsets.size), offsets, 0.0, 0.0, 90.0], True),
    }
    fields = {
        name: empymod.bipole(rec=receiver, mrec=magnetic, **common)
        for name, (receiver, magnetic) in receivers.items()
    }

    comments = [
        "# Weiss, Kalscheuer & Ren (2023), Scenario 2 at 1000 Hz.",
        "# Independent 1-D layered-earth solution from empymod 2.6.0 bipole/QWE.",
        "# exp(+i*omega*t), x-directed 200 m cable, current 0.5 A (moment 100 A m).",
        "# Air: 1e-8 S/m; earth: 1e-4 S/m; 500-1000 m layer: 0.01 S/m.",
        "# Source and receivers are at z=0; receivers are cross-line at x=0.",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        stream.write("\n".join(comments) + "\n")
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["offset_m", "Ex_real", "Ex_imag", "Hy_real",
                         "Hy_imag", "Hz_real", "Hz_imag"])
        for index, offset in enumerate(offsets):
            writer.writerow([
                f"{offset:.0f}",
                f"{fields['Ex'][index].real:.16e}",
                f"{fields['Ex'][index].imag:.16e}",
                f"{fields['Hy'][index].real:.16e}",
                f"{fields['Hy'][index].imag:.16e}",
                f"{fields['Hz'][index].real:.16e}",
                f"{fields['Hz'][index].imag:.16e}",
            ])


if __name__ == "__main__":
    main()
