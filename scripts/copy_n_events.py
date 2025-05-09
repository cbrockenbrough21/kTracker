import ROOT
import os

def copy_first_n_events(input_file, output_file, n_events):
    f_in = ROOT.TFile.Open(input_file, "READ")
    tree_in = f_in.Get("tree")
    if not tree_in:
        raise RuntimeError("Tree 'tree' not found in input file.")
    print(f"Total events: {tree_in.GetEntries()}")


    f_out = ROOT.TFile.Open(output_file, "RECREATE")
    tree_out = tree_in.CloneTree(0)  # Empty tree with same structure

    for i in range(min(n_events, tree_in.GetEntries())):
        tree_in.GetEntry(i)
        tree_out.Fill()

    f_out.Write()
    f_out.Close()
    f_in.Close()
    
    # Show file size
    size_bytes = os.path.getsize(output_file)
    size_mb = size_bytes / (1024 ** 2)
    print(f"Saved first {n_events} events to {output_file}")
    print(f"Output file size: {size_mb:.2f} MB")


# Example usage:
if __name__ == "__main__":
    copy_first_n_events("/project/ptgroup/Catherine/kTracker/noisy_data_gen/noisy_output.root", "data/small_noisy_output.root", 10000)
