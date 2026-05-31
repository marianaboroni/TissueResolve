#!/usr/bin/env python
"""cell2location spatial runner (graceful). Reads prepared inputs; records GPU
availability; supports --fast (reduced iters / CPU). Fails without aborting the
benchmark when cell2location or a GPU is unavailable."""
from __future__ import annotations
import argparse, importlib.util, json, sys, time
from pathlib import Path

OUT = Path("benchmarks/outputs/spatial")


def _meta(status, **kw):
    d = OUT / "method_metadata"; d.mkdir(parents=True, exist_ok=True)
    base = dict(method="cell2location", version="", executed=False, imported=False,
                exported_only=False, status=status, runtime_seconds=0.0,
                input_normalization_used="counts", reference_level_used="fine",
                command_run="run_cell2location.py", warnings="", error_message="",
                output_path="")
    base.update(kw)
    (d / "cell2location.json").write_text(json.dumps(base, indent=2))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args(argv)
    if importlib.util.find_spec("cell2location") is None:
        _meta("skipped",
              error_message="cell2location not installed (see benchmarks/envs/)")
        print("cell2location: skipped (not installed)"); return 0
    gpu = False
    try:
        import torch; gpu = bool(torch.cuda.is_available())
    except Exception:
        pass
    _meta("not_implemented_here",
          warnings=f"installed; runner stub (fast={args.fast}, gpu={gpu})")
    print(f"cell2location: installed; runner stub (gpu={gpu})"); return 0


if __name__ == "__main__":
    raise SystemExit(main())
