"""Lazy reader for MATLAB v7.3 / HDF5 formatted data files."""

from rc2_formatted_data_reader.reader import FormattedDataReader
from rc2_formatted_data_reader.trial_conditions import (
    CONDITION_MAP,
    StimulusLookup,
    parse_cloud_name,
)
from rc2_formatted_data_reader.masks import (
    motion_mask,
    stationary_mask,
    treadmill_motion_mask,
)

__all__ = [
    "FormattedDataReader",
    "CONDITION_MAP",
    "StimulusLookup",
    "parse_cloud_name",
    "motion_mask",
    "stationary_mask",
    "treadmill_motion_mask",
]
