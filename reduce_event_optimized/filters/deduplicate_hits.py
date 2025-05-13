import numpy as np

def deduplicate_hits(detectorIDs, elementIDs, keep_idx):
    """
    Vectorized deduplication based on (detectorID, elementID) pairs.
    Keeps the first occurrence of each pair in keep_idx.
    """

    # Extract the relevant detector and element IDs
    det = detectorIDs[keep_idx]
    ele = elementIDs[keep_idx]

    # Create a compound key array
    combined = det.astype(np.int64) * 1000 + ele.astype(np.int64)  # Assumes <1000 elementIDs

    # Get unique keys and first indices
    _, unique_idx = np.unique(combined, return_index=True)

    # Return the original indices corresponding to the first unique hits
    return keep_idx[unique_idx]