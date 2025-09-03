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

    # Build working hit list with position - only add chamber hits
    working_hits = []
    new_keep_idx = []
    for i in keep_idx:
        if hitlist[i][0] > 30:
            new_keep_idx.append(i)
        else:
            detID, elemID = hitlist[i]
            if detID not in geom.detectors:
                #print(f"WARNING: detID {detID} missing from geometry, skipping")
                continue
            pos = geom.detectors[detID].get_wire_position(elemID)
            working_hits.append((i, detID, pos))

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
    
    seen = set()

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
                if i_d3 in seen and i_d2 in seen and i_d1 in seen:
                    continue
                z1 = geom.get_plane_position(detID1)

                pos_exp_target = SAGITTA_TARGET_CENTER * s2_target + slope_target * (z1 - Z_TARGET)
                pos_exp_dump   = SAGITTA_DUMP_CENTER   * s2_dump   + slope_dump   * (z1 - Z_DUMP)

                win_target = abs(s2_target * SAGITTA_TARGET_WIDTH)
                win_dump   = abs(s2_dump   * SAGITTA_DUMP_WIDTH)

                p_min = min(pos_exp_target - win_target, pos_exp_dump - win_dump)
                p_max = max(pos_exp_target + win_target, pos_exp_dump + win_dump)

                if p_min < pos1 < p_max:
                    for idx in (i_d3, i_d2, i_d1):
                        if idx not in seen:
                            seen.add(idx)
                            new_keep_idx.append(idx)
                          
    return new_keep_idx