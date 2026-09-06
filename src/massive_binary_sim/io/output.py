"""HDF5 output with full run metadata + provenance."""
from __future__ import annotations

import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from .. import __version__


def git_hash():
    try:
        h = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                    cwd=Path(__file__).parent, stderr=subprocess.DEVNULL)
        return h.decode().strip()
    except Exception:
        return "UNKNOWN(no-git)"


def run_metadata(cfg, wallclock_s, validation_results=None):
    return dict(
        code_version=__version__,
        git_hash=git_hash(),
        python=sys.version,
        platform=platform.platform(),
        timestamp_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        seed=cfg.numerics.seed,
        wallclock_s=wallclock_s,
        config=cfg.to_dict(),
        unverified_free_params=cfg.unverified_params(),
        validation=validation_results or {},
    )


def write_hdf5(path, result, cfg, validation_results=None):
    import h5py
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = run_metadata(cfg, result.wallclock_s, validation_results)
    with h5py.File(path, "w") as f:
        f.attrs["metadata_json"] = json.dumps(meta, default=str)
        f.attrs["fate"] = result.fate
        g = f.create_group("history")
        for k, v in result.history.as_arrays().items():
            arr = v
            if arr.dtype.kind in ("U", "O"):
                arr = np.asarray(arr, dtype="S32")
            g.create_dataset(k, data=arr)
        f.create_dataset("events_json", data=np.bytes_(json.dumps(result.events, default=str)))
        f.create_dataset("final_json", data=np.bytes_(json.dumps(result.final, default=str)))
        f.create_dataset("warnings_json", data=np.bytes_(json.dumps(result.warnings)))
    return path
