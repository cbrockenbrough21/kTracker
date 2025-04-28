"""
Decluster Adjacent Hits
"""

def decluster_hits(detector_ids, element_ids, drift_distances, indices):
    """
    Retain the best hit within clusters of adjacent hits.

    Args:
        detector_ids (list[int])
        element_ids (list[int])
        drift_distances (list[float])
        indices (list[int]): Indices of hits to consider

    Returns:
        list[int]: Filtered list of indices to keep
    """
    result = []
    cluster = []

    for i in indices:
        if not cluster:
            cluster.append(i)
            continue

        last = cluster[-1]
        same_detector = detector_ids[i] == detector_ids[last]
        close = abs(element_ids[i] - element_ids[last]) <= 1

        if same_detector and close:
            cluster.append(i)
        else:
            # End of cluster
            result.append(min(cluster, key=lambda j: drift_distances[j]))
            cluster = [i]

    if cluster:
        result.append(min(cluster, key=lambda j: drift_distances[j]))

    return result
