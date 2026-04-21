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

## Run the GLM pipeline — CLI

```bash
rc2-glm <formatted.mat> <out_dir> \
    --backend irls \
    --mc-sequence motion_cloud_sequence_250414.mat \
    --mc-folders image_folders.mat
```

Flags:
- `--backend irls|nemos` (default `irls`; `nemos` absorbs the constant
  `log(bin_width)` into the intercept)
- `--all-clusters` — include non-VISp clusters (default is VISp only)
- `--prefilter` — run the stationary-vs-motion Wilcoxon test and keep
  only clusters the MATLAB decision tree would run a GLM on
- `--mc-sequence` / `--mc-folders` — stimulus lookup for V/VT params

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
| `pipeline` | Section 9 CSV exports (+ CLI `main`) |

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
