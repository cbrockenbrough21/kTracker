import uproot
import numpy as np


def analyze_reduction(original_file, noisy_file, reduced_file):
    tree_orig = uproot.open(original_file)["tree"]
    tree_noisy = uproot.open(noisy_file)["tree"]
    tree_reduced = uproot.open(reduced_file)["tree"]

    n_events = min(len(tree_orig), len(tree_noisy), len(tree_reduced))

    total_real_hits = 0
    total_noise_hits = 0
    total_final_hits = 0
    noise_removed = 0
    real_lost = 0

    for i in range(n_events):
        orig_det = tree_orig["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        orig_elem = tree_orig["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        noisy_det = tree_noisy["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        noisy_elem = tree_noisy["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        red_det = tree_reduced["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        red_elem = tree_reduced["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        # Sets of hits
        real_hits = set(zip(orig_det, orig_elem))
        all_noisy_hits = set(zip(noisy_det, noisy_elem))
        final_hits = set(zip(red_det, red_elem))

        noise_hits = all_noisy_hits - real_hits

        total_real_hits += len(real_hits)
        total_noise_hits += len(noise_hits)
        total_final_hits += len(final_hits)

        noise_removed += len(noise_hits - final_hits)
        real_lost += len(real_hits - final_hits)

    print(f"Events analyzed: {n_events}")
    print(f"Total real hits in original: {total_real_hits}")
    print(f"Total noise hits injected: {total_noise_hits}")
    print(f"Total hits in reduced file: {total_final_hits}")
    print(f"Noise hits removed: {noise_removed}")
    print(f"Real hits lost: {real_lost}")
    
    if total_noise_hits > 0:
        print(f"Noise removal efficiency: {noise_removed / total_noise_hits:.2%}")
    else:
        print("Noise removal efficiency: N/A (no noise hits)")
        
    if total_real_hits > 0:
        print(f"Real hit preservation rate: {(total_real_hits - real_lost) / total_real_hits:.2%}")
    else:
        print("Real hit preservation rate: N/A")


if __name__ == "__main__":
    original_file = "/project/ptgroup/Catherine/kTracker/data/small_raw/small_combined.root"
    noisy_file = "/project/ptgroup/Catherine/kTracker/data/noisy/small_combined_varied_hits.root" 
    reduced_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/small_combined_cleaned.root" 
    analyze_reduction(original_file, noisy_file, reduced_file)
  