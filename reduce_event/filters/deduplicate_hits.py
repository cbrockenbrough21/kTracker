def deduplicate_hits(detectorIDs, elementIDs, keep_idx):
    """
    Remove duplicate hits based on (detectorID, elementID) pairs.
    Keeps the first occurrence.

    Args:
        detectorIDs (list[int])
        elementIDs (list[int])
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
