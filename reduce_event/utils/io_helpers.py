import ROOT
from ROOT import std

def write_reduced_to_root_all_branches(input_filename, output_filename, index_data):
    """
    Writes a new ROOT file with all branches preserved, but specific hit-level branches
    ('detectorID', 'elementID', 'driftDistance', 'tdcTime') filtered using keep_idx.
    """
    # Open input file
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree_in = input_file.Get("tree")
    if not tree_in:
        raise RuntimeError("Could not find 'tree' in input ROOT file.")

    # Prepare output file
    output_file = ROOT.TFile.Open(output_filename, "RECREATE", "", 1)
    output_file.SetCompressionLevel(5)
    output_file.cd()

    # Clone tree structure only (no entries yet)
    tree_out = tree_in.CloneTree(0)

    # Input vectors
    detectorID = std.vector("int")()
    elementID = std.vector("int")()
    driftDistance = std.vector("double")()
    tdcTime = std.vector("double")()

    tree_in.SetBranchAddress("detectorID", detectorID)
    tree_in.SetBranchAddress("elementID", elementID)
    tree_in.SetBranchAddress("driftDistance", driftDistance)
    tree_in.SetBranchAddress("tdcTime", tdcTime)

    # Output vectors
    det_out = std.vector("int")()
    ele_out = std.vector("int")()
    drift_out = std.vector("double")()
    tdc_out = std.vector("double")()

    tree_out.SetBranchAddress("detectorID", det_out)
    tree_out.SetBranchAddress("elementID", ele_out)
    tree_out.SetBranchAddress("driftDistance", drift_out)
    tree_out.SetBranchAddress("tdcTime", tdc_out)

    for entry in index_data:
        i = entry["entry"]
        keep_idx = entry["keep_idx"]

        tree_in.GetEntry(i)

        # Fill output vectors
        det_out.clear()
        ele_out.clear()
        drift_out.clear()
        tdc_out.clear()

        for j in keep_idx:
            det_out.push_back(detectorID[j])
            ele_out.push_back(elementID[j])
            drift_out.push_back(driftDistance[j])
            tdc_out.push_back(tdcTime[j])

        tree_out.Fill()

    output_file.cd()
    tree_out.Write("", ROOT.TObject.kOverwrite)
    output_file.Close()
    input_file.Close()

    print(f"Wrote reduced ROOT file to '{output_filename}'")
