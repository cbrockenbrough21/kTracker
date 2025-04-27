import numpy as np
import random
import numpy as np
import ROOT
from ROOT import std

def reduce_event(detectorIDs, driftDistances, tdcTimes,
                 elementIDs=None, inTimes=None,
                 outoftime=False, realization=False,
                 chamber_eff=0.94, chamber_resol=0.04,
                 dedup=False, decluster=False,
                 nChamberPlanes=30):
    """
    Filters hits based on:
    - out-of-time removal
    - realization (efficiency drop + smearing)
    - deduplication (by detectorID + elementID)
    - declustering (keep smallest drift distance among adjacent hits)
    
    Returns list of indices to keep.
    """
    keep_idx = []
    for i, det_id in enumerate(detectorIDs):
         # Out-of-time hit removal
        if outoftime and inTimes and not inTimes[i]:
            continue
        
        # Realization: simulate random inefficiency and drift smearing
        if realization and det_id <= nChamberPlanes:
            if random.random() > chamber_eff:
                continue
            driftDistances[i] += np.random.normal(0, chamber_resol)

        # Keep the hit
        keep_idx.append(i)

    # Sort for consistent cluster/dedup behavior
    if elementIDs:
        keep_idx.sort(key=lambda i: (detectorIDs[i], elementIDs[i]))

        # Deduplication: keep only one hit per (detectorID, elementID)
        if dedup:
            seen = set()
            keep_idx = [i for i in keep_idx if (key := (detectorIDs[i], elementIDs[i])) not in seen and not seen.add(key)]

        # Decluster: keep hit with smallest driftDistance among adjacent hits
        if decluster:
            keep_idx = _decluster_hits(detectorIDs, elementIDs, driftDistances, keep_idx)

    return keep_idx

def _decluster_hits(detectorIDs, elementIDs, driftDistances, indices):
    """
    Within clusters of adjacent hits (same detectorID and close elementIDs),
    retain only the hit with the smallest drift distance.
    """

    result = []       # Final indices of hits to keep
    cluster = []      # Temporary storage for a cluster of adjacent hits

    for i in indices:
        if not cluster:
            # Start new cluster
            cluster.append(i)
            continue

        last = cluster[-1]
        same_detector = detectorIDs[i] == detectorIDs[last]
        close = abs(elementIDs[i] - elementIDs[last]) <= 1  # "adjacent" in terms of strip/paddle ID

        if same_detector and close:
            # Current hit is part of the ongoing cluster
            cluster.append(i)
        else:
            # End of cluster: choose best (smallest drift) and reset cluster
            result.append(min(cluster, key=lambda j: driftDistances[j]))
            cluster = [i]

    # Final cluster (if it exists)
    if cluster:
        result.append(min(cluster, key=lambda j: driftDistances[j]))

    return result

def run_reduce_event_on_file(root_filename, **kwargs):
    """
    Processes a ROOT file and applies reduce_event() to each event.
    Returns a list of dicts: { "eventID": ..., "keep_idx": [...], "entry": index }
    """
    f = ROOT.TFile.Open(root_filename, "READ")
    tree = f.Get("tree")
    if not tree:
        raise RuntimeError(f"Could not find 'tree' in {root_filename}")

    max_events = kwargs.pop("max_events", None)
    index_data = []

    for i, event in enumerate(tree):
        if max_events and i >= max_events:
            break

        detectorIDs = list(event.detectorID)
        driftDistances = list(event.driftDistance)
        tdcTimes = list(event.tdcTime)
        elementIDs = list(event.elementID)
        inTimes = [True] * len(detectorIDs)  # Placeholder

        keep_idx = reduce_event(
            detectorIDs, driftDistances, tdcTimes,
            elementIDs, inTimes, **kwargs
        )

        index_data.append({
            "entry": i,
            "eventID": int(event.eventID),
            "keep_idx": keep_idx
        })

        print(f"Event {i} (ID {event.eventID}): {len(detectorIDs)} → {len(keep_idx)} hits kept")

    f.Close()
    return index_data

def write_reduced_to_root_all_branches(input_filename, output_filename, index_data, branches_to_filter):
    """
    Writes a new ROOT file where all branches are preserved, but specific hit-level branches are filtered
    using provided keep_idx lists per event.
    """
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("tree")
    if not tree:
        raise RuntimeError("Could not find 'tree' in input ROOT file.")

    output_file = ROOT.TFile(output_filename, "RECREATE")
    new_tree = tree.CloneTree(0)  # Clone structure, no entries

    # Setup references to new vectors
    branch_buffers = {}
    tree.GetEntry(0)
    for name in branches_to_filter:
        print(f"Inspecting branch: {name}")
        sample_value = getattr(tree, name)[0]
        if isinstance(sample_value, int):
            buffer = std.vector('int')()
        elif isinstance(sample_value, float):
            buffer = std.vector('double')()
        else:
            raise TypeError(f"Unsupported type for branch {name}")
        branch_buffers[name] = buffer
        new_tree.SetBranchAddress(name, buffer)

     # Loop over each event and fill the tree with filtered vectors
    for entry in index_data:
        i = entry["entry"]
        keep_idx = entry["keep_idx"]

        tree.GetEntry(i)

        # Replace contents of specified vector branches
        for name in branches_to_filter:
            source = getattr(tree, name)
            buffer = branch_buffers[name]
            buffer.clear()
            for j in keep_idx:
                buffer.push_back(source[j])

        new_tree.Fill()

    output_file.cd()
    new_tree.Write("", ROOT.TObject.kOverwrite)
    output_file.Close()
    input_file.Close()

    print(f"Wrote reduced ROOT file: {output_filename}")


if __name__ == "__main__":
    input_file = "MC_negMuon_Dump_Feb21.root"
    output_file = "reduced_output_all_branches.root"

    index_data = run_reduce_event_on_file(
        root_filename=input_file,
        outoftime=True,
        realization=True,
        dedup=True,
        decluster=True,
        max_events=1000
    )

    # Hit-level branches that should be filtered
    branches_to_filter = [
        "hitID",
        "hit_trackID",
        "processID",
        "detectorID",
        "elementID",
        "driftDistance",
        "tdcTime"
    ]

    write_reduced_to_root_all_branches(
        input_filename=input_file,
        output_filename=output_file,
        index_data=index_data,
        branches_to_filter=branches_to_filter
    )
