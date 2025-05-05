import ROOT
from ROOT import TFile

def compare_tree_entries(file1, file2, entry_index, branches_to_check):
    f1 = TFile.Open(file1, "READ")
    f2 = TFile.Open(file2, "READ")
    t1 = f1.Get("tree")
    t2 = f2.Get("tree")

    print(f"\n🔍 Comparing entry {entry_index} between:")
    print(f"  File 1: {file1}")
    print(f"  File 2: {file2}")

    t1.GetEntry(entry_index)
    t2.GetEntry(entry_index)

    for branch in branches_to_check:
        try:
            v1 = list(getattr(t1, branch))
            v2 = list(getattr(t2, branch))
        except Exception as e:
            print(f"  ❌ Error accessing branch '{branch}': {e}")
            continue

        match = (v1 == v2)
        print(f"  🟢 Branch '{branch}': {'MATCH' if match else 'DIFFER'}")
        if not match:
            print(f"    File1: {v1}")
            print(f"    File2: {v2}")

    f1.Close()
    f2.Close()

# Example usage
file1 = "./../noisy_data_gen/noisy_MC_negMuon_Dump_Feb21.root"
file2 = "./cleaned_output.root"
entry_index = 1
branches = ["eventID", "detectorID", "elementID", "tdcTime", "driftDistance"]

compare_tree_entries(file1, file2, entry_index, branches)
