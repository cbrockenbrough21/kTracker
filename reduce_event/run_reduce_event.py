"""
Run full hit reduction pipeline.
"""

import ROOT
import numpy as np
import time
from filters.out_of_time_removal import remove_out_of_time_hits
from filters.decluster_hits import decluster_hits
from filters.deduplicate_hits import deduplicate_hits
from reduce_event.filters.hodo_mask import hodo_mask
from reduce_event.filters.sagitta import sagitta_reducer
from utils.io_helpers import write_reduced
from geom.geom_service import GeometryService


BRANCHES_TO_FILTER = ["detectorID", "elementID", "driftDistance", "tdcTime"] #"hitID", "hit_trackID", "processID"

def reduce_event(detectorIDs, driftDistances, tdcTimes, elementIDs, **kwargs):
    geom = kwargs.get("geom", None)

    det   = np.asarray(detectorIDs)
    elem  = np.asarray(elementIDs)
    drift = np.asarray(driftDistances, dtype=float)
    tdc   = np.asarray(tdcTimes, dtype=float)

    # Sort by detectorID (primary) then elementID (secondary)
    perm = np.lexsort((elem, det))          # sorted_index -> original_index

    det_s, elem_s = det[perm], elem[perm]
    drift_s, tdc_s = drift[perm], tdc[perm]

    keep_idx = list(range(perm.size))  # indices in SORTED space [0..N)

    if kwargs.get('dedup', False):
        keep_idx = deduplicate_hits(det_s.tolist(), elem_s.tolist(), keep_idx)
    if kwargs.get('outoftime', False):
        keep_idx = remove_out_of_time_hits(tdc_s.tolist(), keep_idx)
    if kwargs.get('decluster', False):
        keep_idx = decluster_hits(det_s.tolist(), elem_s.tolist(),
                                  drift_s.tolist(), tdc_s.tolist(),
                                  geom, keep_idx)
    if kwargs.get('hodomask', False):
        hodo_ids = kwargs.get('hodo_ids', set())
        keep_idx = hodo_mask(det_s.tolist(), elem_s.tolist(), geom, hodo_ids, keep_idx)
    if kwargs.get('sagitta', False):
        keep_idx = sagitta_reducer(det_s.tolist(), elem_s.tolist(), geom, keep_idx)

    # Map from sorted indices -> original indices
    keep_idx_orig = perm[np.asarray(keep_idx, dtype=int)].tolist()
    return keep_idx_orig


def run_reduction(input_file, output_file, tsv_path, **kwargs):
    """
    Read ROOT file, apply reduce_event, and write new ROOT file.
    """
    total_start = time.perf_counter()
    f = ROOT.TFile.Open(input_file, "READ")
    tree = f.Get("tree")
    if not tree:
        raise RuntimeError(f"Could not find 'tree' in {input_file}")

    index_data = []
    
    geom = None
    HODO_IDS = set()
    if kwargs.get('hodomask', False) or kwargs.get('sagitta', False) or kwargs.get('decluster', False):
        geom = GeometryService(tsv_path=tsv_path)
        geom.load_geometry_from_tsv()
        geom.dump_geometry_summary()
        HODO_IDS = {31, 32, 37, 38, 39, 40}

    read_filter_start = time.perf_counter()
    
    for i, event in enumerate(tree):
        detectorIDs = list(event.detectorID)
        driftDistances = list(event.driftDistance)
        tdcTimes = list(event.tdcTime)
        elementIDs = list(event.elementID)

        keep_idx = reduce_event(
            detectorIDs, driftDistances, tdcTimes, elementIDs,
            geom=geom, hodo_ids=HODO_IDS, **kwargs
        )
        
        index_data.append({
            "entry": i,
            "keep_idx": keep_idx
        })

    f.Close()
    read_filter_end = time.perf_counter()

    write_start = time.perf_counter()
    write_reduced(input_file, output_file, index_data)
    write_end = time.perf_counter()

    total_end = time.perf_counter()

    print("\n--- Timing Summary ---")
    print(f"Reduction time: {read_filter_end - read_filter_start:.2f} s")
    print(f"Write time:     {write_end - write_start:.2f} s")
    print(f"Total runtime:  {total_end - total_start:.2f} s")

if __name__ == "__main__":
    # input_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_negMuon_Dump_Feb21_10000_noisy.root" 
    # output_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/MC_negMuon_Dump_Feb21_10000_cleaned.root" 

    input_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_clusters.root"
    output_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/MC_JPsi_Pythia8_Target_April17_10000_noisy_clusters_cleaned.root" 

    run_reduction(
        input_file=input_file,
        output_file=output_file,
        tsv_path = "/project/ptgroup/Catherine/kTracker/reduce_event/geom/data/param.tsv",
        outoftime=False,
        dedup=False,
        decluster=True,
        hodomask=False,
        sagitta=False
    )
