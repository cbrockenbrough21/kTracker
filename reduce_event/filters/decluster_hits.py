import numpy as np
from reco_constants import N_CHAMBER_PLANES

# Set this to True only if you want to "fix" the C++ quirk.
_FLUSH_FINAL_CLUSTER = False

def decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, geom, keep_idx):
    """
    Python port closely following the provided C++:
      - Assumes hits are already sorted by (detectorID, elementID) in keep_idx order.
      - Builds clusters in that order.
      - Breaks when detectorID > nChamberPlanes (no flush of pending cluster).
      - 2-hit rule: DRIFT rule first, then D3p timing rule.
      - >=3-hit rule: mean adjacent ΔTDC; drop all if <10, else keep ends.
    Returns original indices (relative to the passed arrays).
    """
    n_chamber_planes = N_CHAMBER_PLANES

    result = []
    cluster = []

    def process_cluster(cluster_local):
        """Return kept indices from this cluster (C++ processCluster)."""
        m = len(cluster_local)
        if m == 0:
            return []
        if m == 1:
            return [cluster_local[0]]

        # All hits in cluster are from the same detector
        i0 = cluster_local[0]
        det_id = int(detectorIDs[i0])

        if m == 2:
            i0, i1 = cluster_local
            # Pair span (C++: 0.9 * 0.5 * (pos_back - pos_front)), no abs, no swap
            pos0 = geom.detectors[int(detectorIDs[i0])].get_wire_position(int(elementIDs[i0]))
            pos1 = geom.detectors[int(detectorIDs[i1])].get_wire_position(int(elementIDs[i1]))
            w_max = 0.9 * 0.5 * (pos1 - pos0)
            w_min = (w_max / 9.0) * 4.0

            drift0, drift1 = float(driftDistances[i0]), float(driftDistances[i1])
            tdc0,   tdc1   = float(tdcTimes[i0]),       float(tdcTimes[i1])

            # 1) DRIFT rule first — keep the smaller drift
            if (drift0 > w_max and drift1 > w_min) or (drift1 > w_max and drift0 > w_min):
                # C++ erases the larger; ties keep the first (front)
                return [i0] if drift0 <= drift1 else [i1]

            # 2) D3p timing rule — det 19..24 and |ΔTDC| < 8 → drop both
            if 19 <= det_id <= 24 and abs(tdc0 - tdc1) < 8.0:
                return []

            # Otherwise keep both
            return [i0, i1]

        # m >= 3: mean of adjacent TDC differences (using incoming order)
        tdcs = [float(tdcTimes[i]) for i in cluster_local]
        dt_mean = float(np.mean(np.abs(np.diff(tdcs)))) if m >= 2 else 0.0

        if dt_mean < 10.0:
            # Electric noise — discard all
            return []
        else:
            # Keep first and last; drop middle
            return [cluster_local[0], cluster_local[-1]]

    # --- Build clusters in the order provided by keep_idx (already sorted upstream) ---
    it = iter(keep_idx)
    for i in it:
        # HODO BREAK: mirror C++ behavior
        if int(detectorIDs[i]) > n_chamber_planes:
            # If we are not flushing, we must KEEP the pending cluster as-is (C++ keeps it)
            if _FLUSH_FINAL_CLUSTER:
                result.extend(process_cluster(cluster))
            else:
                result.extend(cluster)           # <-- keep pending cluster unchanged

            # And we must KEEP the rest of the hits unchanged (C++ never touched them)
            result.append(i)                     # include current i
            result.extend(list(it))              # include all remaining indices
            cluster = []
            break

        if not cluster:
            cluster = [i]
            continue

        last = cluster[-1]
        same_det = int(detectorIDs[i]) == int(detectorIDs[last])
        adj_elem = (int(elementIDs[i]) - int(elementIDs[last])) <= 1

        if not same_det or not adj_elem:
            result.extend(process_cluster(cluster))
            cluster = [i]
        else:
            cluster.append(i)

    # END-OF-LOOP: mirror C++ “no final flush” behavior
    if cluster:
        if _FLUSH_FINAL_CLUSTER:
            result.extend(process_cluster(cluster))
        else:
            result.extend(cluster)               # <-- keep pending final cluster unchanged

    return result

 
