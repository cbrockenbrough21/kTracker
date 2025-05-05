import numpy as np

def decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx):
    """
    Group adjacent hits and apply declustering logic to remove noise-like clusters.
    
    Args:
        detectorIDs, elementIDs, driftDistances, tdcTimes: Full hit info arrays
        keep_idx (list[int]): Indices of hits that passed previous filters (e.g. out-of-time)

    Returns:
        list[int]: Indices of hits to keep after declustering
    """

    # Step 1: Sort the hits (so clusters of adjacent elements in the same detector are grouped)
    sorted_hits = sorted(keep_idx, key=lambda i: (detectorIDs[i], elementIDs[i]))

    result = []    # Final filtered indices
    cluster = []   # Temporary cluster of adjacent hits

    for i in sorted_hits:
        if not cluster:
            # First hit starts a new cluster
            cluster.append(i)
            continue

        last = cluster[-1]
        # Check if current hit is adjacent to last hit in the same detector
        if detectorIDs[i] == detectorIDs[last] and abs(elementIDs[i] - elementIDs[last]) <= 1:
            cluster.append(i)
        else:
            # Cluster ended — process it and start a new one
            result.extend(process_cluster(cluster, detectorIDs, driftDistances, tdcTimes))
            cluster = [i]

    # Process the final cluster at end of loop
    if cluster:
        result.extend(process_cluster(cluster, detectorIDs, driftDistances, tdcTimes))

    return result


def process_cluster(cluster, detectorIDs, driftDistances, tdcTimes):
    """
    Apply station-specific rules to determine which hits in a cluster to keep.

    Rules:
    - If 1 hit → keep it.
    - If 2 hits:
        - If D3p (detectors 19–24) and |tdc diff| < 8 → remove both (likely noise)
        - Else, use drift distances to keep one or both.
    - If >=3 hits:
        - If avg tdc diff < 10 → remove all (electronic noise)
        - Else, keep first and last hit only
    """

    # Half-cell widths by detector group, used for drift-based filtering
    half_cell_widths = {
        0: 0.635 / 2,  # D0
        1: 2.083 / 2,  # D2X
        2: 2.021 / 2,  # D2U
        3: 2.021 / 2,  # D2V
        4: 2.000 / 2,  # D3p
        5: 2.000 / 2   # D3m
    }

    n = len(cluster)
    if n == 1:
        return [cluster[0]]  # Nothing to decluster

    det_id = detectorIDs[cluster[0]]
    same_detector = all(detectorIDs[i] == det_id for i in cluster)

    if not same_detector:
        # Safety check: not a real cluster — return all
        return cluster

    hw = half_cell_widths.get(det_id // 5, 0.635 / 2)  # Fall back to D0 width

    if n == 2:
        i0, i1 = cluster
        drift0, drift1 = driftDistances[i0], driftDistances[i1]
        tdc0, tdc1 = tdcTimes[i0], tdcTimes[i1]

        # D3p (detectorIDs 19–24): remove both hits if TDCs are too close
        if 19 <= det_id <= 24 and abs(tdc0 - tdc1) < 8:
            return []

        # Drift-based logic: one large + one medium drift → remove weaker one
        w_max = 0.9 * hw
        w_min = 0.4 * hw

        if drift0 > w_max and drift1 > w_min:
            return [i1]  # Keep smaller drift
        elif drift1 > w_max and drift0 > w_min:
            return [i0]
        else:
            return [i0, i1]

    # n ≥ 3: check average TDC time difference
    elif n >= 3:
        tdcs = [tdcTimes[i] for i in cluster]
        dt_mean = np.mean(np.abs(np.diff(tdcs)))

        if dt_mean < 10:
            return []  # Likely noise — remove all
        else:
            return [cluster[0], cluster[-1]]  # Keep edges only
