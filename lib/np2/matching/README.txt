Cross-session unit tracking for Neuropixels 2.0 data
================================================================

Track the same units across recordings of one animal, using UnitMatchPy's
official SpikeInterface integration.

  match_sessions.py   -- build UnitMatch inputs from each session's
                          SortingAnalyzer, then run UnitMatchPy matching

For environment setup and installation, see the main repo README:
  ../../../README.md

For the sorting stage that must run first, see:
  ../sorting/README.txt


======================================================================
 ENVIRONMENT -- matching-specific tools
======================================================================

Runs in the same conda environment as sorting (`spikeinterface`).

  pip install -U UnitMatchPy
  git clone https://github.com/EnnyvanBeest/UnitMatch    (e.g. into C:\Users\Lab\SWC\UnitMatch)

The clone is only needed for its UnitMatchPy source (imported directly, not
via the pip package, so a `git pull` there picks up upstream fixes without a
reinstall). The default clone location is set in match_sessions.py as
DEFAULT_UNITMATCH_REPO; override per run with --unitmatch-repo.

DeepUnitMatch (the neural-network matcher) is NOT used here: UnitMatchPy's
official SpikeInterface integration (make_UnitMatch_folder_from_sorting_
analyzers) only exists for classic UnitMatchPy, not for DeepUnitMatch, which
still expects hand-built RawWaveforms/ files. Revisit if that changes upstream.

Each session's SortingAnalyzer for matching is built on the bandpass_only
recording (see ../sorting/README.txt) rather than the destriped one used for
sorting: destriping deliberately removes the cross-channel spatial structure
UnitMatch matches on. It is built sparse (radius_um=150), matching
UnitMatchPy's own internal channel radius -- no accuracy gain from a denser
footprint, only extra computation.

"Good" units for matching are read from each session's selected_clusters.txt
-- the pipeline's actual curation result (Bombcell + ClusterSelectionTable,
plus any manual keep/discard edit in clusters_to_check.csv), not Bombcell's
own default thresholds recomputed on the spot. See
apply_selected_clusters_labels() in match_sessions.py.


======================================================================
 RUNNING
======================================================================

From MATLAB (recommended -- sorts any missing session automatically):

    ctl = RC2Preprocess();
    ctl.match_sessions({probe_id_1, probe_id_2, ...})   % chronological order, N >= 2

See 'help RC2Preprocess.match_sessions'.

Standalone (all sessions must already be sorted):

    conda activate spikeinterface
    cd C:\Users\Lab\SWC\rc2_analysis_si_ks4+bc_um\lib\np2\matching
    python match_sessions.py <imec0_ks4_dir_1> <imec0_ks4_dir_2> [...]

What happens:
  1. checks every session has a completed sort (sorting_analyzer/,
     cluster_groups.csv, selected_clusters.txt, bandpass_only.ap/);
  2. builds a sparse SortingAnalyzer per session on its bandpass_only
     recording, and a UnitMatch input folder from it
     (make_UnitMatch_folder_from_sorting_analyzers);
  3. overwrites each session's bombcell_labels.tsv so "good" means "in this
     session's selected_clusters.txt", not Bombcell's own default thresholds;
  4. runs UnitMatchPy (waveform properties -> metric scores -> Naive Bayes ->
     matches -> unique IDs across sessions);
  5. saves results.

Options:
    --save-dir DIR         output folder (default <parent of session 1>\unit_match)
    --threshold 0.75       match-probability threshold
    --unitmatch-repo DIR   clone of EnnyvanBeest/UnitMatch (default already set)

Outputs (in the save dir):
    MatchTable.csv          unit pairs with match probability
    MatchingOverview.png    total-score / probability / final-match matrices
    ClusInfo.pickle, UMparam.pickle, MatchProb.npy, Matches.npy,
    UM Scores.npz, WaveformInfo.npz   intermediate data

More than 2 sessions: pass all of them in one call. UnitMatchPy's own
assign_unique_id does the transitive closure across the whole group (so a
10-day chronic experiment is one call with 10 session paths, not 9 pairwise
runs) -- match order follows the order the paths are given in, so always
list sessions chronologically.
