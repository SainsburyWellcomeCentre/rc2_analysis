"""Generate the regression-test baseline CSV from the current 4-probe run.

Pulls the 12 (probe_id, cluster_id) pairs listed in
``python/tests/test_pipeline_regression.REGRESSION_CLUSTERS`` from the
aggregated ``glm/current/glm_model_comparison.csv`` and writes
``python/tests/data/regression_baseline.csv``.

Run this:
- once after the first 4-probe production run with the new defaults
  (history ON, onset OFF) to bootstrap the test;
- whenever the model intentionally changes (and we have explicit
  buy-in to bump the baseline).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python" / "tests"))
from test_pipeline_regression import REGRESSION_CLUSTERS  # noqa: E402

CURRENT_CSV = Path(
    "/Users/lauraporta/local_data/motion_clouds/figures/glm/current"
    "/glm_model_comparison.csv"
)
OUT_CSV = REPO / "python" / "tests" / "data" / "regression_baseline.csv"


def main() -> int:
    if not CURRENT_CSV.is_file():
        print(f"ERROR: {CURRENT_CSV} does not exist. Run rc2-glm first.",
              file=sys.stderr)
        return 1

    df = pd.read_csv(CURRENT_CSV)
    keep_cols = ["probe_id", "cluster_id", "time_Selected_cv_bps",
                 "time_selected_vars", "time_n_spikes"]
    df = df[keep_cols]

    keys = pd.DataFrame(
        list(REGRESSION_CLUSTERS), columns=["probe_id", "cluster_id"]
    )
    baseline = keys.merge(df, on=["probe_id", "cluster_id"], how="left")

    missing = baseline[baseline["time_Selected_cv_bps"].isna()]
    if len(missing):
        print(f"ERROR: {len(missing)} regression clusters missing from "
              f"{CURRENT_CSV.name}:\n{missing[['probe_id', 'cluster_id']]}",
              file=sys.stderr)
        return 2

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    baseline.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(baseline)} rows)")
    print()
    print("=== baseline contents ===")
    print(baseline.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
