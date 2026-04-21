# rc2_python

Python side of `rc2_analysis`. One distribution (`rc2_python`) shipping
two packages:

- **`rc2_formatted_data_reader`** — lazy reader for MATLAB v7.3
  formatted `.mat` files (spikes, trials, per-sample traces, stimulus
  lookup).
- **`rc2_glm`** — replication of the MATLAB Poisson GLM pipeline in
  `scripts/glm_single_cluster_analysis.m` (v10 Hierarchical Forward
  Selection). Imports the reader above.

## Install

Requires a conda env (the project uses `rc2_analysis`).

```bash
conda activate rc2_analysis
pip install -e python
```

One install covers both packages. The CLI entry point `rc2-glm` is
registered on the PATH of the active env.

## Configure local paths (.env)

Copy `python/.env.example` → `python/.env` and fill in your paths:

```bash
cp python/.env.example python/.env
```

Recognised keys:

| Key | Meaning |
|---|---|
| `RC2_FORMATTED_DATA_DIR` | Directory of formatted `.mat` files. Bare-filename CLI args are resolved here; with no CLI arg the first `*.mat` is used. |
| `RC2_FORMATTED_DATA_PATH` | Pin a single `.mat` file (takes precedence over `DIR`). |
| `RC2_MC_SEQUENCE_PATH` | Motion-cloud sequence `.mat` (V/VT stimulus params). |
| `RC2_MC_FOLDERS_PATH` | `image_folders.mat` (paired with the sequence). |
| `RC2_GLM_OUTPUT_DIR` | Default output directory for CSVs. |

`python/.env` is gitignored. CLI flags always override anything in `.env`.

## Run the GLM pipeline — CLI

With `python/.env` configured, just run:

```bash
rc2-glm                 # picks first *.mat in RC2_FORMATTED_DATA_DIR
rc2-glm probe01.mat     # bare filename is resolved against the same dir
rc2-glm /abs/probe.mat /abs/out
```

Flags:
- `--backend irls|nemos` (default `irls`; `nemos` absorbs the constant
  `log(bin_width)` into the intercept)
- `--all-clusters` — include non-VISp clusters (default is VISp only)
- `--prefilter` — run the stationary-vs-motion Wilcoxon test and keep
  only clusters the MATLAB decision tree would run a GLM on
- `--mc-sequence` / `--mc-folders` — override stimulus lookup paths
- `--no-plots` — skip matplotlib figures (CSVs + log still written)
- `--plot-format pdf|png|both` (default `pdf`)
- `--plot-clusters N` — cap per-cluster plots to the first N clusters

## Run the GLM pipeline — Python

```python
from rc2_glm import GLMConfig
from rc2_glm.pipeline import run_pipeline
from rc2_formatted_data_reader import StimulusLookup

lookup = StimulusLookup(
    sequence_mat="~/local_data/motion_clouds/motion_cloud_sequence_250414.mat",
    folders_mat="~/local_data/motion_clouds/image_folders.mat",
)
run_pipeline(
    mat_path="~/local_data/motion_clouds/formatted_data/<probe>.mat",
    config=GLMConfig(),
    output_dir="./glm_out",
    stimulus_lookup=lookup,
    backend="irls",
)
```

Outputs in `out_dir`:

- `glm_model_comparison.csv` — one row per cluster (null / selected /
  additive / full-interaction CV bits-per-spike, selected variables,
  tuning flags, interaction flags)
- `glm_selection_history.csv` — one row per forward-selection round
- `glm_coefficients.csv` — one row per fitted coefficient
- `prefilter_decision_tree.csv` — only when `--prefilter` is used
- `run.log` — stdout mirror with per-section banners and per-cluster
  timing (always written when `out_dir` is set)
- `figs/*.pdf` — matplotlib figures (unless `--no-plots`). Add
  `--plot-format png` (or `both`) to also emit PNG rasters.

## Plots

`figs/` contains:

| File | MATLAB ref | Contents |
|---|---|---|
| `basis_functions.pdf` | Fig 1 | Speed / TF / SF-OR / onset bases |
| `forward_selection_summary.pdf` | Fig 2 | Cross-cluster null-vs-selected, #vars, variable counts |
| `cluster_<id>_overview.pdf` | Fig 4 (script 2727-3290) | 4×6 per-cluster grid: β swarm (col 1) + predicted-vs-observed scatter by {condition, speed, TF, SF, OR}. Selected row outlined in red with a ★ prefix; R² in each scatter title. Observed FR is per-trial boxcar-smoothed (100 ms) and binned trial-first (≥2 trials per cell). |
| `cluster_<id>_tuning.pdf` | Fig 5 (script 3298-3980) | 5×4 reconstructed tuning curves (Observed + Null / Selected / Additive / FullInteraction × Speed / TF / SF / OR); per-condition traces (T_Vstatic green, V orange, VT blue; black stationary dot at speed=0). Uses the stored training β — no refit. |
| `cluster_<id>_selection.pdf` | slice of Fig 4 | CV bps across forward-selection rounds |
| `cluster_<id>_coefficients.pdf` | slice of Fig 4 | Fitted β ± SE for Selected and Additive |

Not yet ported: Fig 3 (speed-profile CV panels), Fig 6 (tuning-curve
correlation box plots), Fig 7 (selected-vs-additive paired comparison),
Fig 8 (trial-level predictions). The Python pipeline doesn't yet
export the per-fold / per-trial traces those plots need.

`plot_tuning_curves` and `plot_cluster_model_overview` both reuse the
**trained** β vectors and column names stored on `ClusterFit`
(`model_betas`, `model_col_names`, `model_predictions`). No refit
happens at plot time. The tuning-curve sweeps build a prediction design
with ``is_prediction=True`` and align its columns by name against the
stored training design — missing columns are filled with zeros, extras
dropped — so ``X_aligned @ β`` stays well-defined regardless of which
columns training dropped under ``_drop_zero_variance``. The onset kernel
is evaluated at `steady_state_time = 1.5 s` for motion predictions and at
zeros for the stationary dot, matching MATLAB Section 8b.

## Read formatted data directly

```python
from rc2_formatted_data_reader import FormattedDataReader

with FormattedDataReader("<probe>.mat") as r:
    print(r.probe_id, r.fs, r.n_clusters, r.n_trials)
    spikes = r.spike_times(cluster_idx=0)      # lazy
    s, e = r.trial_bounds(trial_idx=0)
    v = r.velocity(s, e)                        # sliced, no full load
```

## Layout

```
python/
├── pyproject.toml          # single distribution
├── README.md
├── src/
│   ├── rc2_formatted_data_reader/
│   └── rc2_glm/
└── tests/
    ├── conftest.py
    ├── test_reader.py
    └── test_glm_*.py
```

## rc2_glm module map

| Module | Mirrors MATLAB |
|---|---|
| `config` | parameter block at top of the script |
| `io` | `FormattedData` + `Trial.treadmill_motion_mask` |
| `time_binning` | Section 3b (lines ~830-1107) |
| `basis` | `make_raised_cosine_basis`, `make_onset_kernel_basis` |
| `design_matrix` | `assemble_design_matrix(_selected)` |
| `fitting` | `fit_poisson_glm` (IRLS is a direct port) |
| `cross_validation` | `cross_validate_glm` |
| `forward_selection` | `forward_select_model` |
| `prefilter` | `is_stationary_vs_motion_significant` + decision tree |
| `pipeline` | Section 9 CSV exports + logging (+ CLI `main`) |
| `plots` | Plotly ports of Fig 1 / Fig 2 / Fig 5 and Fig 4 slices |

## Tests

From `python/`:

```bash
conda activate rc2_analysis
pytest
```

Reader tests skip cleanly without a `.mat` fixture. Set
`RC2_FORMATTED_DATA_PATH` or `RC2_FORMATTED_DATA_DIR` to point at a
real file, or drop one under `~/local_data/motion_clouds/formatted_data/`.

## Known divergences from MATLAB

- NeMoS backend absorbs a constant `log(bin_width)` offset into the
  intercept. For uniform 100 ms bins this is exact; for non-uniform
  bins use the IRLS backend.
- `AlignedTrial` offset and the analysis-window solenoid logic from
  `FormattedData.trial_by_id` are not yet replicated — trial bounds
  come straight from `FormattedData.probe_t` slicing.
- Plot coverage: Figs 1, 2, 4, 5 + waterfall/coefficient summaries are
  ported; Figs 3, 6, 7, 8 are not.
- Fig 4 colormaps: MATLAB uses `parula` for the VT speed scatter;
  Python uses `viridis` (perceptually uniform stand-in). SF and OR use
  the same fixed RGB palette as MATLAB
  (`sf_unique_all = [0.003, 0.006, 0.012]`,
  `or_unique_all = [-π/4, 0, π/4, π/2]`).
