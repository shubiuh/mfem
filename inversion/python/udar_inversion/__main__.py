import argparse
import json

from .workflow import run_inversion, synthetic_recovery


def main():
    parser = argparse.ArgumentParser(
        description="MFEM UDAR inversion using a C++ ROL or PETSc/TAO backend")
    parser.add_argument("mode", choices=["invert", "synthetic"])
    parser.add_argument("configuration", help="JSON configuration file")
    parser.add_argument("--observed", help="udar3p tensor CSV (required for invert)")
    parser.add_argument("--output", required=True, help="New directory for this run")
    parser.add_argument("--backend", choices=["rol", "tao"],
                        help="Override configuration backend (default: rol)")
    parser.add_argument("--library", help="Path to the selected backend library")
    args = parser.parse_args()
    if args.mode == "invert" and not args.observed:
        parser.error("invert requires --observed")
    if args.mode == "synthetic":
        result = synthetic_recovery(args.configuration, args.output,
                                    library=args.library, backend=args.backend)
    else:
        result = run_inversion(args.configuration, args.observed, args.output,
                               library=args.library, backend=args.backend)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result.get("recovery_passed", result["converged"]) else 2


if __name__ == "__main__":
    raise SystemExit(main())
