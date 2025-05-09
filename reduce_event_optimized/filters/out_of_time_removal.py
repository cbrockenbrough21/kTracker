"""
Out-of-Time Hit Removal Filter
"""

import numpy as np



def remove_out_of_time_hits(tdc_times, keep_idx):
    """
    Filters hits based on TDC center and width.
    Args:
        tdc_times (np.array): Full array of TDC times
        keep_idx (list): List of indices to check and filter
    Returns:
        list: Updated keep_idx
    """
    # tdc_center=950.0
    # tdc_width=7.4
    
    # lower = tdc_center - (tdc_width / 2)
    # upper = tdc_center + (tdc_width / 2)

    # return [i for i in keep_idx if lower <= tdc_times[i] <= upper]
    
    return keep_idx

# def remove_out_of_time_hits(detectorid, elementid, tdctime, keep_idx):
#     """
#     Filters hits based on detector-specific TDC timing windows.
#     Only keeps indices from keep_idx that are in time.

#     Args:
#         detectorid (array-like): detector IDs for all hits
#         elementid (array-like): element IDs (not used here, but passed for compatibility)
#         tdctime (array-like): TDC times for all hits
#         keep_idx (list[int]): current list of indices to be filtered

#     Returns:
#         list[int]: filtered list of indices with only in-time hits
#     """
#     result = []
#     for i in keep_idx:
#         det = int(detectorid[i])
#         tdc = tdctime[i]

#         if (det < 7 and 1700 < tdc < 1820) \
#            or (12 < det < 19 and 1450 < tdc < 1710) \
#            or (18 < det < 25 and 1360 < tdc < 1580) \
#            or (24 < det < 31 and 1490 < tdc < 1700) \
#            or (46 < det < 55 and 560 < tdc < 1200) \
#            or (30 < det < 47):  # hodoscopes
#             result.append(i)

#     return result

