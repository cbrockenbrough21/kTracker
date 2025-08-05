from reco_constants import (
    Z_TARGET, Z_DUMP,
    SAGITTA_TARGET_CENTER, SAGITTA_DUMP_CENTER,
    SAGITTA_TARGET_WIDTH, SAGITTA_DUMP_WIDTH,
    TX_MAX
)

def sagitta_reducer(hitlist, geom, keep_idx):
    """
    Python translation of EventReducer::sagittaReducer (C++).
    Only applies sagitta to chamber detectors (ID 1–30).

    Args:
        hitlist: list of (detectorID, elementID) pairs from raw input
        geom: GeometryService instance (provides plane and wire info)
        keep_idx: list of indices into hitlist to consider
    
    Returns:
        list of indices (subset of keep_idx) that passed sagitta filtering
    """
    # Split chambers vs non-chambers explicitly
    chambers = [i for i in keep_idx if hitlist[i][0] <= 30]
    nonchambers = [i for i in keep_idx if hitlist[i][0] > 30]

    # Build working hit list with position - only add chamber hits
    working_hits = []
    for idx, i in enumerate(chambers):
        detID, elemID = hitlist[i]
        if detID not in geom.detectors:
            #print(f"WARNING: detID {detID} missing from geometry, skipping")
            continue
        pos = geom.detectors[detID].get_wire_position(elemID)
        working_hits.append((len(working_hits), detID, pos))

    # Partition hits by detectorID (D1/D2/D3)
    detectorID_st1_max = 12
    detectorID_st2_max = 18
    d1_hits, d2_hits, d3_hits = [], [], []

    for i, detID, pos in working_hits:
        if detID <= detectorID_st1_max:
            d1_hits.append((i, detID, pos))
        elif detID <= detectorID_st2_max:
            d2_hits.append((i, detID, pos))
        else:
            d3_hits.append((i, detID, pos))

    # Init keep flags
    flag = [-1] * len(working_hits)

    # Sagitta triplet logic
    for i_d3, detID3, pos3 in d3_hits:
        z3 = geom.get_plane_position(detID3)
        slope_target = pos3 / (z3 - Z_TARGET)
        slope_dump   = pos3 / (z3 - Z_DUMP)

        for i_d2, detID2, pos2 in d2_hits:
            if geom.get_plane_type(detID3) != geom.get_plane_type(detID2):
                continue

            z2 = geom.get_plane_position(detID2)
            if abs((pos3 - pos2) / (z2 - z3)) > TX_MAX:
                continue

            s2_target = pos2 - slope_target * (z2 - Z_TARGET)
            s2_dump   = pos2 - slope_dump   * (z2 - Z_DUMP)

            for i_d1, detID1, pos1 in d1_hits:
                if geom.get_plane_type(detID3) != geom.get_plane_type(detID1):
                    continue
                if flag[i_d3] > 0 and flag[i_d2] > 0 and flag[i_d1] > 0:
                    continue

                z1 = geom.get_plane_position(detID1)

                pos_exp_target = SAGITTA_TARGET_CENTER * s2_target + slope_target * (z1 - Z_TARGET)
                pos_exp_dump   = SAGITTA_DUMP_CENTER   * s2_dump   + slope_dump   * (z1 - Z_DUMP)

                win_target = abs(s2_target * SAGITTA_TARGET_WIDTH)
                win_dump   = abs(s2_dump   * SAGITTA_DUMP_WIDTH)

                p_min = min(pos_exp_target - win_target, pos_exp_dump - win_dump)
                p_max = max(pos_exp_target + win_target, pos_exp_dump + win_dump)

                if p_min < pos1 < p_max:
                    flag[i_d3] = flag[i_d2] = flag[i_d1] = 1

    # flagged chambers
    idx_D3 = len(d1_hits) + len(d2_hits) + len(d3_hits)
    kept_chambers = [chambers[i] for i, f in enumerate(flag) if f >= 0 and i < idx_D3]

    # untouched non-chambers
    return kept_chambers + nonchambers