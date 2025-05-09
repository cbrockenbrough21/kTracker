import numpy as np

def decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx):
    """
    Groups adjacent hits by (detectorID, elementID) and applies declustering logic
    to remove potential electronic noise or overlapping spurious hits.

    Args:
        detectorIDs, elementIDs, driftDistances, tdcTimes: np.arrays of full hit info.
        keep_idx (list[int]): Indices of hits that passed previous filters (e.g., out-of-time)

    Returns:
        list[int]: Indices of hits to keep after declustering
    """

    # Step 1: Sort hits by detectorID, then elementID to group adjacent hits together
    sorted_hits = sorted(keep_idx, key=lambda i: (detectorIDs[i], elementIDs[i]))

    result = []    # Final filtered hit indices to return
    cluster = []   # Temporarily holds adjacent hits for declustering

    for i in sorted_hits:
        if not cluster:
            cluster.append(i)
            continue

        last = cluster[-1]
        # If the hit is in the same detector and element ID is adjacent, keep clustering
        if detectorIDs[i] == detectorIDs[last] and abs(elementIDs[i] - elementIDs[last]) <= 1:
            cluster.append(i)
        else:
            # Cluster ended — process and start a new one
            result.extend(process_cluster(cluster, detectorIDs, driftDistances, tdcTimes))
            cluster = [i]

    # Process the final cluster (if any)
    if cluster:
        result.extend(process_cluster(cluster, detectorIDs, driftDistances, tdcTimes))

    return result

def process_cluster(cluster, detectorIDs, driftDistances, tdcTimes):
    """
    Applies declustering rules to a group of hits that are adjacent in elementID
    and share the same detectorID.

    Args:
        cluster (list[int]): Indices of hits in the same detector and adjacent elements
        detectorIDs, driftDistances, tdcTimes: full hit arrays

    Returns:
        list[int]: Subset of input cluster to keep
    """

    # Detector-specific cell widths (used for drift distance comparisons)
    half_cell_widths = {
        0: 0.635 / 2,  # Station 0 (D0)
        1: 2.083 / 2,  # D2X
        2: 2.021 / 2,  # D2U
        3: 2.021 / 2,  # D2V
        4: 2.000 / 2,  # D3p (important for special timing cut)
        5: 2.000 / 2   # D3m
    }

    n = len(cluster)
    if n == 1:
        return [cluster[0]]  # A single hit isn't declustered

    det_id = detectorIDs[cluster[0]]

    # Sanity check: ensure all hits in cluster share the same detector
    if not all(detectorIDs[i] == det_id for i in cluster):
        return cluster

    # Default to D0 width if not found
    hw = half_cell_widths.get(det_id // 5, 0.635 / 2)

    # ---- Two-hit cluster logic ----
    if n == 2:
        i0, i1 = cluster
        drift0, drift1 = driftDistances[i0], driftDistances[i1]
        tdc0, tdc1 = tdcTimes[i0], tdcTimes[i1]

        # D3p timing noise removal: discard both if TDCs are very close
        if 19 <= det_id <= 24 and abs(tdc0 - tdc1) < 8:
            return []

        # Drift-based preference: remove larger drift if it's noisy-looking
        w_max = 0.9 * hw
        w_min = 0.4 * hw

        if drift0 > w_max and drift1 > w_min:
            return [i1]
        elif drift1 > w_max and drift0 > w_min:
            return [i0]
        else:
            return [i0, i1]

    # ---- Larger clusters (n >= 3) ----
    tdcs = [tdcTimes[i] for i in cluster]
    dt_mean = np.mean(np.abs(np.diff(np.sort(tdcs))))

    if dt_mean < 10:
        return []  # Very likely electronic noise: discard entire cluster
    else:
        return [cluster[0], cluster[-1]]  # Retain edges
