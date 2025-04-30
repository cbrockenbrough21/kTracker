"""
Decluster adjacent hits based on drift distances and TDC timing.
"""

import numpy as np

def decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx):
    """
    Apply declustering to hits.
    
    Args:
        detectorIDs (list[int]): detector IDs
        elementIDs (list[int]): element IDs
        driftDistances (list[float]): drift distances
        tdcTimes (list[float]): tdc times
        keep_idx (list[int]): indices of hits to consider

    Returns:
        list[int]: Updated list of indices to keep
    """

    # Detector half-cell widths (in cm) for windowing
    half_cell_widths = {
        0: 0.635 / 2,  # D0
        1: 2.083 / 2,  # D2X
        2: 2.021 / 2,  # D2U
        3: 2.021 / 2,  # D2V
        4: 2.000 / 2,  # D3p
        5: 2.000 / 2,  # D3m
    }

    # Sort hits by detectorID and elementID
    sorted_hits = sorted(keep_idx, key=lambda i: (detectorIDs[i], elementIDs[i]))

    result = []  # Final hits to keep
    cluster = []

    def flush_cluster(cluster):
        if not cluster:
            return

        # If only 1 hit, keep it
        if len(cluster) == 1:
            result.append(cluster[0])
            return

        # Get cluster properties
        det_ids = [detectorIDs[i] for i in cluster]
        drifts = [driftDistances[i] for i in cluster]
        tdcs = [tdcTimes[i] for i in cluster]

        same_detector = all(det_ids[0] == det for det in det_ids)

        if not same_detector:
            for i in cluster:
                result.append(i)
            return

        if len(cluster) == 2:
            det_id = det_ids[0]
            
            # Special handling for D3p (detectorIDs 19-24)
            if 19 <= det_id <= 24 and abs(tdcs[0] - tdcs[1]) < 8:
                # Electronic noise: remove both hits
                return
            else:
                # Drift-based removal
                hw = half_cell_widths.get(det_id // 5, 0.635 / 2)
                w_max = 0.9 * hw
                w_min = (w_max / 9) * 4

                if drifts[0] > w_max and drifts[1] > w_min:
                    result.append(cluster[0])
                elif drifts[1] > w_max and drifts[0] > w_min:
                    result.append(cluster[1])
                else:
                    result.extend(cluster)

        elif len(cluster) >= 3:
            dt_sum = sum(abs(tdcs[i] - tdcs[i-1]) for i in range(1, len(tdcs)))
            dt_mean = dt_sum / (len(tdcs) - 1)

            if dt_mean < 10:
                # Electronic noise cluster: remove all
                return
            else:
                # Keep only first and last hit
                result.append(cluster[0])
                result.append(cluster[-1])

    # Walk through sorted hits and group clusters
    for i in sorted_hits:
        if not cluster:
            cluster.append(i)
            continue

        last = cluster[-1]
        if detectorIDs[i] == detectorIDs[last] and abs(elementIDs[i] - elementIDs[last]) <= 1:
            cluster.append(i)
        else:
            flush_cluster(cluster)
            cluster = [i]

    # Final cluster
    flush_cluster(cluster)

    return result
