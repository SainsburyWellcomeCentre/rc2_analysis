This directory contains copies of python scripts which were modified to use the fork of ecephys_spike_sorting (https://github.com/jenniferColonell/ecephys_spike_sorting, commit 60c40251ed568fb036b4364e615a261d3afc4800).

We used this to process our Neuropixels 2.0, 4 shank data.

The file 'spikeGLX_pipeline_np2.py' is used by our analysis scripts. It is copied and written to the ecephys_spike_sorting_np2\ecephys_spike_soriting\scripts directory before ecephys_spike_sorting is run.

The following files are here just for reference, and were modified after cloning the ecephys_spike_sorting repo, to run on our computer.

create_input_json:   from ecephys_spike_sorting_np2\ecephys_spike_sorting\scripts\create_input_json.py
SpikeGLX_utils.py:   from ecephys_spike_sorting_np2\ecephys_spike_sorting\scripts\helpers\SpikeGLX_utils.py
