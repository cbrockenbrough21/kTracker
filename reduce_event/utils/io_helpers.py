"""
Helpers for reading and writing ROOT files with filtered hits.
"""

import ROOT
from ROOT import std

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