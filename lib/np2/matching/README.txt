Cross-session unit tracking for Neuropixels 2.0 (4-shank) data
================================================================

Track the same units across recordings of one animal, using DeepUnitMatch
(recommended) or classic UnitMatchPy (fallback).

Scripts in this directory:

  extract_raw_waveforms.py   -- prep: write RawWaveforms/ for (Deep)UnitMatch
  run_deep_unit_match.py     -- cross-session tracking with DeepUnitMatch (recommended)
  run_unit_match.py          -- cross-session tracking with classic UnitMatchPy (fallback)
  unitmatch_plots.py         -- plotting helpers for matching outputs

For environment setup and installation, see the main repo README:
  ../../../README.md

For the sorting stage that must run first, see:
  ../sorting/README.txt


======================================================================
 ENVIRONMENT -- matching-specific tools
======================================================================

Runs in the same conda environment as sorting (`spikeinterface`, Python 3.10).

  TOOL              USED FOR                         INSTALLED AS
  ----------------- -------------------------------- ----------------------------
  UnitMatchPy       cross-session matching + I/O      pip install UnitMatchPy
                    (also pulls mtscomp, mat73)
  DeepUnitMatch     deep-NN matching + trained model  git clone of UnitMatch repo
                    (model is NP2.0 4-shank only)
  PyTorch (CUDA)    GPU for DeepUnitMatch              pip install torch (CUDA build,
                                                        same as used for KS4)

Notes:
  * UnitMatchPy + DeepUnitMatch are imported from the cloned UnitMatch repo (so we
    get the DeepUnitMatch package + pretrained model), but `pip install UnitMatchPy`
    is still the easy way to pull the runtime deps (mtscomp, mat73, scikit-learn...).
  * extract_raw_waveforms.py reads from the CatGT bin (band-pass only) -- not the
    destriped recording used for sorting -- so the cross-channel spatial footprint
    that UnitMatch matches on is preserved. Called automatically by the matchers.

INSTALL (once you already have the `spikeinterface` env from the sorting stage):
   pip install UnitMatchPy            (matching + I/O; pulls mtscomp, mat73, ...)
   pip install torch                  (only if not already installed for KS4)
   git clone https://github.com/EnnyvanBeest/UnitMatch    (e.g. into C:\Users\Lab\SWC\UnitMatch)

   The clone provides the DeepUnitMatch package and the pretrained model
   (UnitMatchPy/DeepUnitMatch/utils/model, trained for Npix 2.0 4-shank).
   The default clone location is set in the scripts as DEFAULT_UNITMATCH_REPO;
   override per run with --unitmatch-repo.

   If you did NOT `pip install UnitMatchPy`, install its deps directly:
     pip install mtscomp mat73 scikit-learn joblib h5py

VERIFY (DeepUnitMatch; set the path to your clone):
   python -c "import sys; sys.path.insert(0, r'C:\Users\Lab\SWC\UnitMatch\UnitMatchPy'); \
              from DeepUnitMatch.testing import test; test.load_trained_model(device='cpu'); print('DeepUnitMatch OK')"


======================================================================
 RUNNING
======================================================================

Do this only after ALL recordings of one animal are sorted (see ../sorting).
This step is OPTIONAL and INDEPENDENT -- it never runs during sorting.

You do NOT edit any paths. Point the script at ONE folder that contains all
the recordings you want to match; it finds every imec*_ks4 session under it:

    conda activate spikeinterface
    cd C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\np2\matching
    python run_deep_unit_match.py  D:\path\to\recordings_root

What happens:
  1. discovers every sorted session (imec*_ks4) under the folder and PRINTS the
     list in match order (sorted by path = chronological if folders are named by
     date/sequence -- always check the printed list);
  2. extracts each session's RawWaveforms/ from its CatGT bin + KS4 sorting
     (automatic; skipped if already present);
  3. runs DeepUnitMatch (model -> similarity -> drift-correct + Bayes -> matches);
  4. saves results to <recordings_root>\unit_match_deep\.

Options:
    --save-dir DIR        output folder (default <recordings_root>\unit_match_deep)
    --threshold 0.5       match-probability threshold
    --dist-thresh 50      max drift-corrected centroid distance (um)
    --include-mua         match mua units too (default: good units only)
    --re-extract          force-rebuild RawWaveforms/
    --unitmatch-repo DIR  clone of EnnyvanBeest/UnitMatch (default already set)

Classic UnitMatchPy instead of the deep model (same interface):
    python run_unit_match.py  D:\path\to\recordings_root

Rebuild RawWaveforms only (rarely needed; the matchers do it automatically):
    python extract_raw_waveforms.py  D:\path\to\recordings_root

Outputs (in the save dir):
    MatchTable.csv          unit pairs with match probability
    UniqueIDConversion.*    cluster IDs with cross-session unique IDs
    MatchingOverview.png    similarity / probability / final-match matrices

Requirements recap: each session must already have KS4 output and Bombcell
labels (cluster_group.tsv); RawWaveforms are generated for you. The trigger /
sorting itself is NOT re-run -- matching uses the spike times and clusters
from the sorting stage plus the CatGT voltage.
