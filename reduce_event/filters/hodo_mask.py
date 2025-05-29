# filters/hodo_mask.py
from geom.geom_service import GeometryService

def extract_hodo_hits(detectorIDs, elementIDs, hodo_ids, keep_idx):
    """
    Extracts unique hodo hit UIDs (detectorID * 1000 + elementID) from hits that passed previous filters.
    """
    return {detectorIDs[i] * 1000 + elementIDs[i] for i in keep_idx if detectorIDs[i] in hodo_ids}


def apply_hodo_mask(detectorIDs, elementIDs, hodo_uids, c2h, keep_idx):
    """
    Applies hodo masking to chamber hits using a chamber-to-hodo LUT.

    - If a chamber UID has no masking info in c2h, the hit is kept (flexible mode).
    - Otherwise, the chamber hit is kept if it overlaps any detected hodo UID.
    """
    new_keep_idx = []
    for i in keep_idx:
        
        if detectorIDs[i] > 30:
            #print(f"[KEEP NON-CHAMBER] detID={detectorIDs[i]}, elemID={elementIDs[i]}")
            new_keep_idx.append(i)  # Always keep non-chamber hits
            continue
    
        uid = detectorIDs[i] * 1000 + elementIDs[i]
        
        if uid not in c2h:
            continue

        overlap = set(c2h[uid]) & hodo_uids
        if overlap:
            #print(f"[MATCH] Chamber UID {uid} matches hodoscope UIDs {overlap} — keeping")
            new_keep_idx.append(i)
    return new_keep_idx


def hodo_mask(detectorIDs, elementIDs, geom, hodo_ids, keep_idx):
    """
    Main entry point for hodoscope masking filter.

    Args:
        detectorIDs (list[int])
        elementIDs (list[int])
        geom (GeometryService): geometry service instance with c2h map
        hodo_ids (set[int]): detector IDs for hodoscopes
        keep_idx (list[int]): indices from previous filters to consider

    Returns:
        list[int]: indices to keep after hodo masking
    """
    hodo_uids = extract_hodo_hits(detectorIDs, elementIDs, hodo_ids, keep_idx)
    return apply_hodo_mask(detectorIDs, elementIDs, hodo_uids, geom.c2h, keep_idx)
