import ROOT
from collections import Counter

# ==============================
# Core: Event-level occupancy cut
# ==============================
def accept_event(detector_ids, max_hits):
    """
    Applies occupancy cuts using detectorID list.

    Parameters:
    - detector_ids: list of int (one detector ID per hit in the event)
    - max_hits: dictionary of max allowed hits per region (D0-D3m)

    Returns:
    - True if the event passes all occupancy cuts
    - False otherwise
    """

    counts = Counter(detector_ids)

    nHitsD0  = sum(counts[i] for i in range(1, 7))     # planes 1–6
    nHitsD1  = sum(counts[i] for i in range(7, 13))    # 7–12
    nHitsD2  = sum(counts[i] for i in range(13, 19))   # 13–18
    nHitsD3p = sum(counts[i] for i in range(19, 25))   # 19–24
    nHitsD3m = sum(counts[i] for i in range(25, 31))   # 25–30

    if nHitsD0  > max_hits["D0"]:  return False
    if nHitsD1  > max_hits["D1"]:  return False
    if nHitsD2  > max_hits["D2"]:  return False
    if nHitsD3p > max_hits["D3p"]: return False
    if nHitsD3m > max_hits["D3m"]: return False

    return True

# ==============================
# Wrapper: Process full ROOT file
# ==============================
def run_accept_event_on_file(root_filename, max_hits):
    """
    Opens a ROOT file and applies accept_event() to each event.

    Parameters:
    - root_filename: path to the input ROOT file
    - max_hits: dictionary with thresholds for D0, D1, D2, D3p, D3m

    Returns:
    - List of accepted event indices
    """
    file = ROOT.TFile.Open(root_filename, "READ")
    tree = file.Get("tree")
    if not tree:
        raise RuntimeError(f"Could not find TTree 'tree' in {root_filename}")

    accepted_indices = []

    for i, event in enumerate(tree):
        detector_ids = list(event.detectorID)
        if accept_event(detector_ids, max_hits):
            accepted_indices.append(i)

    file.Close()
    return accepted_indices

# ==============================
# Optional main for testing
# ==============================
if __name__ == "__main__":
    # Define your max occupancy thresholds
    max_hits = {
        "D0": 40,
        "D1": 40,
        "D2": 40,
        "D3p": 40,
        "D3m": 40
    }

    filename = "MC_negMuon_Dump_Feb21.root"
    accepted = run_accept_event_on_file(filename, max_hits)

    print(f"Accepted {len(accepted)} events")
    print("Accepted event indices (sample):", accepted[:10])
