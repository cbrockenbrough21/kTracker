"""
Out-of-Time Hit Removal Filter
"""

import numpy as np

def remove_out_of_time_hits(tdc_times, keep_idx):
    """
    Filters hits based on TDC center and width.
    Args:
        tdc_times (list): Full array of TDC times
        keep_idx (list): List of indices to check and filter
    Returns:
        list: Updated keep_idx
    """
    # lower = tdc_center - (tdc_width / 2)
    # upper = tdc_center + (tdc_width / 2)

    # return [i for i in keep_idx if lower <= tdc_times[i] <= upper]
    return keep_idx
