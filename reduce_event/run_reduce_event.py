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
from utils.io_helpers import write_reduced
from geom.geom_service import GeometryService


BRANCHES_TO_FILTER = ["detectorID", "elementID", "driftDistance", "tdcTime"] #"hitID", "hit_trackID", "processID"

def reduce_event(detectorIDs, driftDistances, tdcTimes, elementIDs, **kwargs):
    """
    Apply filters to reduce events.

    Args:
        detectorIDs (list[int])
        driftDistances (list[float])
        tdcTimes (list[float])
        elementIDs (list[int])

    Returns:
        list[int]: Indices to keep
    """
    keep_idx = list(range(len(detectorIDs)))
            
    if kwargs.get('dedup', False):
        keep_idx = deduplicate_hits(detectorIDs, elementIDs, keep_idx)

    if kwargs.get('outoftime', False):
        keep_idx = remove_out_of_time_hits(tdcTimes, keep_idx)

    if kwargs.get('decluster', False):
        keep_idx = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
        
    if kwargs.get('hodomask', False):
        geom = kwargs['geom']
        hodo_ids = kwargs.get('hodo_ids', set())
        keep_idx = hodo_mask(detectorIDs, elementIDs, geom, hodo_ids, keep_idx)

        
    return keep_idx

def run_reduction(input_file, output_file, **kwargs):
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
    if kwargs.get('hodomask', False):
        tsv_path = "/project/ptgroup/Catherine/kTracker/reduce_event/geom/data/param.tsv"
        geom = GeometryService(tsv_path=tsv_path)
        geom.load_geometry_from_tsv()
        HODO_IDS = {41} 

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
    input_file = "/project/ptgroup/Catherine/kTracker/data/noisy/small_combined_varied_hits.root" 
    output_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/small_combined_cleaned.root" 

    run_reduction(
        input_file=input_file,
        output_file=output_file,
        outoftime=False,
        dedup=False,
        decluster=False,
        hodomask=True
    )
