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
    Filters hits based on out-of-time, realization efficiency, deduplication, and declustering.
    Returns list of indices to keep.
    """
    keep_idx = []
    for i, det_id in enumerate(detectorIDs):
        if outoftime and inTimes and not inTimes[i]:
            continue
        if realization and det_id <= nChamberPlanes:
            if random.random() > chamber_eff:
                continue
            driftDistances[i] += np.random.normal(0, chamber_resol)
        keep_idx.append(i)

    if elementIDs:
        keep_idx.sort(key=lambda i: (detectorIDs[i], elementIDs[i]))

        if dedup:
            seen = set()
            keep_idx = [i for i in keep_idx if (key := (detectorIDs[i], elementIDs[i])) not in seen and not seen.add(key)]

        if decluster:
            keep_idx = _decluster_hits(detectorIDs, elementIDs, driftDistances, keep_idx)

    return keep_idx

def _decluster_hits(detectorIDs, elementIDs, driftDistances, indices):
    """
    Within clusters of adjacent hits (same detectorID and close elementIDs),
    retain only the hit with the smallest drift distance.
    """
    result = []
    cluster = []

    for i in indices:
        if not cluster:
            cluster.append(i)
            continue

        last = cluster[-1]
        same_detector = detectorIDs[i] == detectorIDs[last]
        close = abs(elementIDs[i] - elementIDs[last]) <= 1

        if same_detector and close:
            cluster.append(i)
        else:
            result.append(min(cluster, key=lambda j: driftDistances[j]))
            cluster = [i]

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


def write_reduced_to_root_from_idx(input_filename, output_filename, index_data):
    """
    Writes a new ROOT file using filtered hit data based on provided keep_idx per event.
    """
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("tree")
    if not tree:
        raise RuntimeError("Could not find 'tree' in input ROOT file.")

    output_file = ROOT.TFile(output_filename, "RECREATE")
    new_tree = ROOT.TTree("tree", "Reduced Event Data")

    eventID = np.zeros(1, dtype=np.int32)
    detectorID = std.vector('int')()
    driftDistance = std.vector('double')()
    tdcTime = std.vector('double')()
    elementID = std.vector('int')()

    new_tree.Branch("eventID", eventID, "eventID/I")
    new_tree.Branch("detectorID", detectorID)
    new_tree.Branch("driftDistance", driftDistance)
    new_tree.Branch("tdcTime", tdcTime)
    new_tree.Branch("elementID", elementID)

    for entry in index_data:
        i = entry["entry"]
        tree.GetEntry(i)

        eventID[0] = entry["eventID"]
        detectorID.clear()
        driftDistance.clear()
        tdcTime.clear()
        elementID.clear()

        for j in entry["keep_idx"]:
            detectorID.push_back(tree.detectorID[j])
            driftDistance.push_back(tree.driftDistance[j])
            tdcTime.push_back(tree.tdcTime[j])
            elementID.push_back(tree.elementID[j])

        new_tree.Fill()

    output_file.cd()
    new_tree.Write("", ROOT.TObject.kOverwrite)
    output_file.Close()
    input_file.Close()

    print(f"Wrote reduced ROOT file: {output_filename}")


if __name__ == "__main__":
    input_file = "MC_negMuon_Dump_Feb21.root"
    output_file = "reduced_output_from_idx.root"

    index_data = run_reduce_event_on_file(
        root_filename=input_file,
        outoftime=True,
        realization=True,
        dedup=True,
        decluster=True,
        max_events=1000
    )

    write_reduced_to_root_from_idx(input_file, output_file, index_data)
