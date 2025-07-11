import uproot
import numpy as np
from collections import defaultdict
from contextlib import redirect_stdout

def analyze_by_station(original_file, noisy_file, reduced_file, max_events=100, output_file="/project/ptgroup/Catherine/kTracker/station_analysis_summary.txt"):
    tree_orig = uproot.open(original_file)["tree"]
    tree_noisy = uproot.open(noisy_file)["tree"]
    tree_reduced = uproot.open(reduced_file)["tree"]

    n_events = min(tree_orig.num_entries, tree_noisy.num_entries, tree_reduced.num_entries, max_events)

    real_hits = defaultdict(int)
    noise_hits = defaultdict(int)
    final_hits = defaultdict(int)
    real_lost = defaultdict(int)
    noise_removed = defaultdict(int)

    for i in range(n_events):
        orig_det = tree_orig["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        orig_elem = tree_orig["elementID"].array(entry_start=i, entry_stop=i+1)[0]
        noisy_det = tree_noisy["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        noisy_elem = tree_noisy["elementID"].array(entry_start=i, entry_stop=i+1)[0]
        red_det = tree_reduced["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        red_elem = tree_reduced["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        real = {(d, e) for d, e in zip(orig_det, orig_elem)}
        noisy = {(d, e) for d, e in zip(noisy_det, noisy_elem)}
        final = {(d, e) for d, e in zip(red_det, red_elem)}
        noise = noisy - real

        for d, _ in real:
            real_hits[d] += 1
        for d, _ in noise:
            noise_hits[d] += 1
        for d, _ in final:
            final_hits[d] += 1
        for d, _ in (real - final):
            real_lost[d] += 1
        for d, _ in (noise - final):
            noise_removed[d] += 1

    def summarize(det_range, label):
        r = sum(real_hits[d] for d in det_range)
        n = sum(noise_hits[d] for d in det_range)
        f = sum(final_hits[d] for d in det_range)
        lost = sum(real_lost[d] for d in det_range)
        removed = sum(noise_removed[d] for d in det_range)
        preserv = (r - lost) / r if r else 1.0
        noise_eff = removed / n if n else 0.0

        print(f"=== {label} Summary ===")
        print(f"Detectors: {list(det_range)}")
        print(f"Total real hits: {r}")
        print(f"Total noise hits: {n}")
        print(f"Total final hits: {f}")
        print(f"Total real hits lost: {lost}")
        print(f"Total noise hits removed: {removed}")
        print(f"Preservation rate: {preserv:.2%}")
        print(f"Noise removal efficiency: {noise_eff:.2%}")
        print("")

    with open(output_file, "w") as f:
        with redirect_stdout(f):
            print(f"Analyzing up to {n_events} events...\n")
            summarize(range(1, 7), "Station 1 (D1)")
            summarize(range(7, 13), "Station 2 (D2)")
            summarize(range(13, 19), "Station 3 (D3)")
            summarize(range(1, 31), "All Chambers")
            nonchamber_ids = [d for d in real_hits if d > 30]
            summarize(nonchamber_ids, "Non-Chamber Detectors")
            print("=== Per-Detector Statistics (All Detectors Present) ===")
            print(f"{'Detector':>9} | {'Real':>6} | {'Noise':>7} | {'Final':>6} | {'Lost':>5} | {'Removed':>7} | {'Preserv%':>9} | {'NoiseRem%':>10}")
            print("-" * 80)

            all_present_detectors = sorted(set(real_hits) | set(noise_hits) | set(final_hits))

            for det in all_present_detectors:
                r, n, f = real_hits[det], noise_hits[det], final_hits[det]
                lost, removed = real_lost[det], noise_removed[det]
                preserv = 100 * (r - lost) / r if r > 0 else 100.0
                noise_eff = 100 * removed / n if n > 0 else 0.0
                print(f"{det:9} | {r:6} | {n:7} | {f:6} | {lost:5} | {removed:7} | {preserv:8.2f}% | {noise_eff:9.2f}%")


# Update these file paths as needed
if __name__ == "__main__":
    original_file = "/project/ptgroup/Catherine/kTracker/data/small_raw/MC_JPsi_Pythia8_Target_April17_10000.root"
    noisy_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root"
    reduced_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/MC_JPsi_Pythia8_Target_April17_10000_onlyElectronic_cleaned.root"
    
    analyze_by_station(original_file, noisy_file, reduced_file, max_events=100)
