import numpy as np

def deduplicate_hits(detectorIDs, elementIDs, keep_idx):
    """
    Remove duplicate hits based on (detectorID, elementID) pairs.
    Keeps the first occurrence.

    Args:
        detectorIDs (np.ndarray)
        elementIDs (np.ndarray)
        keep_idx (list[int])

    Returns:
        list[int]: Deduplicated list of indices
    """
    seen = set()
    result = []
    for i in keep_idx:
        key = (detectorIDs[i], elementIDs[i])
        if key not in seen:
            seen.add(key)
            result.append(i)
    return result

# import numpy as np

# def deduplicate_hits(detectorIDs, elementIDs, keep_idx):
#     """
#     Vectorized deduplication based on (detectorID, elementID) pairs.
#     Keeps the first occurrence.

#     Args:
#         detectorIDs (np.ndarray)
#         elementIDs (np.ndarray)
#         keep_idx (np.ndarray)

#     Returns:
#         np.ndarray: Deduplicated indices
#     """
#     # Extract only the relevant subset
#     det_subset = detectorIDs[keep_idx]
#     elem_subset = elementIDs[keep_idx]

#     # Combine into structured array of (detectorID, elementID)
#     combined = np.core.records.fromarrays([det_subset, elem_subset], names='det,elem')

#     # Find unique pairs and their first occurrence
#     _, unique_idx = np.unique(combined, return_index=True)

#     # Map back to the original keep_idx values
#     return keep_idx[unique_idx]
