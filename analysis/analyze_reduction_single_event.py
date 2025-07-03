import uproot
import numpy as np
from collections import defaultdict

def analyze_event_by_eventid(original_file, noisy_file, reduced_file, target_eventid):
    # open
    tree_orig = uproot.open(original_file)["tree"]
    tree_noisy = uproot.open(noisy_file)["tree"]
    tree_reduced = uproot.open(reduced_file)["tree"]

    # grab the eventID array
    event_ids = tree_orig["eventID"].array()

    # find entry number matching the requested EventID
    matching_entries = np.where(event_ids == target_eventid)[0]
    if len(matching_entries) == 0:
        print(f"EventID {target_eventid} not found in original file.")
        return
    entry = matching_entries[0]
    print(f"EventID {target_eventid} found at TTree entry {entry}")

    # bounds check
    n_total_events = min(tree_orig.num_entries, tree_noisy.num_entries, tree_reduced.num_entries)
    if entry >= n_total_events:
        print(f"Entry {entry} out of range (max available {n_total_events-1})")
        return

    # prepare counters
    real_hits = defaultdict(int)
    noise_hits = defaultdict(int)
    final_hits = defaultdict(int)
    real_lost = defaultdict(int)
    noise_removed = defaultdict(int)

    # read only that event
    orig_det = tree_orig["detectorID"].array(entry_start=entry, entry_stop=entry+1)[0]
    orig_elem = tree_orig["elementID"].array(entry_start=entry, entry_stop=entry+1)[0]
    noisy_det = tree_noisy["detectorID"].array(entry_start=entry, entry_stop=entry+1)[0]
    noisy_elem = tree_noisy["elementID"].array(entry_start=entry, entry_stop=entry+1)[0]
    red_det = tree_reduced["detectorID"].array(entry_start=entry, entry_stop=entry+1)[0]
    red_elem = tree_reduced["elementID"].array(entry_start=entry, entry_stop=entry+1)[0]

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

    # Chamber IDs
    chamber_ids = range(1, 31)
    nonchamber_ids = [d for d in real_hits if d > 30]

    print(f"\n=== Per-Detector Statistics for EventID {target_eventid} (entry {entry}) ===")
    print(f"{'Detector':>9} | {'Real':>6} | {'Noise':>7} | {'Final':>6} | {'Lost':>5} | {'Removed':>7} | {'Preserv%':>9} | {'NoiseRem%':>10}")
    print("-" * 80)

    for det in chamber_ids:
        r, n, f = real_hits[det], noise_hits[det], final_hits[det]
        lost, removed = real_lost[det], noise_removed[det]
        preserv = 100 * (r - lost) / r if r > 0 else 100.0
        noise_eff = 100 * removed / n if n > 0 else 0.0
        print(f"{det:9} | {r:6} | {n:7} | {f:6} | {lost:5} | {removed:7} | {preserv:8.2f}% | {noise_eff:9.2f}%")

    # chamber totals
    total_real = sum(real_hits[d] for d in chamber_ids)
    total_noise = sum(noise_hits[d] for d in chamber_ids)
    total_final = sum(final_hits[d] for d in chamber_ids)
    total_lost = sum(real_lost[d] for d in chamber_ids)
    total_removed = sum(noise_removed[d] for d in chamber_ids)

    print("\n=== Chamber Summary ===")
    print(f"Total real hits: {total_real}")
    print(f"Total noise hits: {total_noise}")
    print(f"Total final hits: {total_final}")
    print(f"Total real hits lost: {total_lost}")
    print(f"Total noise hits removed: {total_removed}")
    print(f"Preservation rate: {(total_real - total_lost) / total_real:.2%}" if total_real else "N/A")
    print(f"Noise removal efficiency: {total_removed / total_noise:.2%}" if total_noise else "N/A")

    # non-chamber summary
    nonch_real = sum(real_hits[d] for d in nonchamber_ids)
    nonch_noise = sum(noise_hits[d] for d in nonchamber_ids)
    nonch_final = sum(final_hits[d] for d in nonchamber_ids)
    nonch_lost = sum(real_lost[d] for d in nonchamber_ids)
    nonch_removed = sum(noise_removed[d] for d in nonchamber_ids)

    print("\n=== Non-Chamber Hit Summary (detectorID > 30) ===")
    print(f"Total real: {nonch_real}, lost: {nonch_lost}")
    print(f"Total noise: {nonch_noise}, removed: {nonch_removed}")
    print(f"Total final hits: {nonch_final}")
    print(f"Real preservation: {(nonch_real - nonch_lost) / nonch_real:.2%}" if nonch_real else "N/A")
    print(f"Noise removal: {nonch_removed / nonch_noise:.2%}" if nonch_noise else "N/A")


if __name__ == "__main__":
    target_eventid = 1270  # change to whatever EventID you want
    original_file = "/project/ptgroup/Catherine/kTracker/data/small_raw/MC_JPsi_Pythia8_Target_April17_10000.root"
    noisy_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root"
    reduced_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/MC_JPsi_Pythia8_Target_April17_10000_onlyElectronic_cleaned.root"
    analyze_event_by_eventid(original_file, noisy_file, reduced_file, target_eventid)
