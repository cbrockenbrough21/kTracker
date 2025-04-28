"""
Out-of-Time Hit Removal Filter
"""

import numpy as np

def remove_out_of_time_hits(tdc_times, tdc_center, tdc_width, keep_idx):
    """
    Filters hits based on TDC center and width.
    Args:
        tdc_times (np.array): Full array of TDC times
        tdc_center (float): Center of the allowed window
        tdc_width (float): Full width
        keep_idx (list): List of indices to check and filter
    Returns:
        list: Updated keep_idx
    """
    lower = tdc_center - (tdc_width / 2)
    upper = tdc_center + (tdc_width / 2)

    return [i for i in keep_idx if lower <= tdc_times[i] <= upper]
