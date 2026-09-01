# rc2_analysis

Preprocessing and analysis pipeline for electrophysiological data acquired with the
rollercoaster (RC2) setup on **Neuropixels 2.0 (4-shank)** probes recorded with **SpikeGLX**.

| Stage | Tool |
|---|---|
| AP band-pass + phase shift, destriping, orchestration, waveforms, quality metrics, Phy export | [SpikeInterface](https://github.com/SpikeInterface/spikeinterface) |
| Spike sorting | [Kilosort 4](https://github.com/MouseLand/Kilosort) |
| Automated curation (good / mua / noise / non-soma labels) | [Bombcell](https://github.com/Julie-Fabre/bombcell) *(via SpikeInterface's `bombcell_label_units`)* |
| Manual cluster inspection (optional) | [Phy](https://github.com/cortex-lab/phy) |
| Cross-session unit tracking (optional) | [UnitMatch](https://github.com/EnnyvanBeest/UnitMatch) *(via SpikeInterface's `SortingAnalyzer` integration)* |

All spike sorting happens in **Python** — MATLAB drives the pipeline and does the downstream
formatting, quality control and analysis.

## Documentation

Installation, configuration and usage are in the **[wiki](../../wiki)**.

## Citing

If you use this pipeline, please cite:

**Velez-Fort, Cossell, Porta, Clopath, Margrie (2025), *Motor and vestibular signals in the 
visual cortex permit the separation of self versus externally generated visual motion*, Cell**

Note that the published analyses were produced with an earlier version of this pipeline 
(Kilosort 2 + `ecephys_spike_sorting`), preserved in the git history.

Please also cite the tools it relies on:

- **SpikeInterface** — Buccino et al. (2020), *eLife* 9:e61834.
- **Kilosort 4** — Pachitariu et al. (2024), *Nature Methods*
  ([github](https://github.com/MouseLand/Kilosort))
- **Bombcell** — Fabre, van Beest, Peters, Carandini & Harris (2023), *Bombcell: automated curation
  and cell classification of spike-sorted electrophysiology data*, Zenodo
  ([10.5281/zenodo.8172821](https://doi.org/10.5281/zenodo.8172821))
- **UnitMatch** — van Beest, Bimbard et al. (2024), *Nature Methods*
  ([10.1038/s41592-024-02440-1](https://www.nature.com/articles/s41592-024-02440-1)) — if you use
  cross-session tracking
- **Phy** ([github](https://github.com/cortex-lab/phy)) — if you used it for manual inspection


