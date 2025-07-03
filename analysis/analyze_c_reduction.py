import uproot
import numpy as np
import re
from collections import defaultdict

def parse_filtered_output(filtered_txt):
    """
    Parses filtered_output.txt
    Returns a dict:
        event_id -> set of (detectorID, elementID)
    """
    events = {}
    current_event_id = None
    pairs = set()
    with open(filtered_txt) as f:
        for line in f:
            line = line.strip()
            if match := re.match(r"RunID:\s*\d+,\s*EventID:\s*(\d+)", line):
                if current_event_id is not None:
                    events[current_event_id] = pairs
                current_event_id = int(match.group(1))
                pairs = set()
            elif re.match(r"^\d+ :", line):
                parts = line.split(":")
                if len(parts) >= 3:
                    det = int(parts[1].strip())
                    elem = int(parts[2].strip())
                    pairs.add((det, elem))
        # final
        if current_event_id is not None:
            events[current_event_id] = pairs
    return events

def analyze_aggregated(filtered_events, original_file, noisy_file):
    tree_orig = uproot.open(original_file)["tree"]
    tree_noisy = uproot.open(noisy_file)["tree"]

    orig_eventids = tree_orig["eventID"].array()
    noisy_eventids = tree_noisy["eventID"].array()

    # accumulate
    real_hits = defaultdict(int)
    noise_hits = defaultdict(int)
    final_hits = defaultdict(int)
    real_lost = defaultdict(int)
    noise_removed = defaultdict(int)

    for event_id, reduced_hits in filtered_events.items():
        # find entry in raw
        entries = np.where(orig_eventids == event_id)[0]
        if len(entries) == 0:
            print(f"EventID {event_id} not found in raw!")
            continue
        entry = entries[0]
        orig_det = tree_orig["detectorID"].array(entry_start=entry, entry_stop=entry+1)[0]
        orig_elem = tree_orig["elementID"].array(entry_start=entry, entry_stop=entry+1)[0]
        noisy_det = tree_noisy["detectorID"].array(entry_start=entry, entry_stop=entry+1)[0]
        noisy_elem = tree_noisy["elementID"].array(entry_start=entry, entry_stop=entry+1)[0]

        real = {(d, e) for d, e in zip(orig_det, orig_elem)}
        noisy = {(d, e) for d, e in zip(noisy_det, noisy_elem)}
        final = reduced_hits
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

    # chamber
    chamber_ids = range(1, 31)

    print("\n=== Aggregated Per-Detector Statistics (Detectors 1–30) ===")
    print(f"{'Detector':>9} | {'Real':>6} | {'Noise':>7} | {'Final':>6} | {'Lost':>5} | {'Removed':>7} | {'Preserv%':>9} | {'NoiseRem%':>10}")
    print("-" * 80)

    for det in chamber_ids:
        r, n, f = real_hits[det], noise_hits[det], final_hits[det]
        lost, removed = real_lost[det], noise_removed[det]
        preserv = 100 * (r - lost) / r if r > 0 else 100.0
        noise_eff = 100 * removed / n if n > 0 else 0.0
        print(f"{det:9} | {r:6} | {n:7} | {f:6} | {lost:5} | {removed:7} | {preserv:8.2f}% | {noise_eff:9.2f}%")

    total_real = sum(real_hits[d] for d in chamber_ids)
    total_noise = sum(noise_hits[d] for d in chamber_ids)
    total_final = sum(final_hits[d] for d in chamber_ids)
    total_lost = sum(real_lost[d] for d in chamber_ids)
    total_removed = sum(noise_removed[d] for d in chamber_ids)

    print("\n=== Aggregated Chamber Summary ===")
    print(f"Total real hits: {total_real}")
    print(f"Total noise hits: {total_noise}")
    print(f"Total final hits: {total_final}")
    print(f"Total real hits lost: {total_lost}")
    print(f"Total noise hits removed: {total_removed}")
    print(f"Preservation rate: {(total_real - total_lost) / total_real:.2%}" if total_real else "N/A")
    print(f"Noise removal efficiency: {total_removed / total_noise:.2%}" if total_noise else "N/A")

    # non-chamber
    nonchamber_ids = [d for d in real_hits if d > 30]
    nonch_real = sum(real_hits[d] for d in nonchamber_ids)
    nonch_noise = sum(noise_hits[d] for d in nonchamber_ids)
    nonch_final = sum(final_hits[d] for d in nonchamber_ids)
    nonch_lost = sum(real_lost[d] for d in nonchamber_ids)
    nonch_removed = sum(noise_removed[d] for d in nonchamber_ids)

    print("\n=== Aggregated Non-Chamber Hit Summary (detectorID > 30) ===")
    print(f"Total real: {nonch_real}, lost: {nonch_lost}")
    print(f"Total noise: {nonch_noise}, removed: {nonch_removed}")
    print(f"Total final hits: {nonch_final}")
    print(f"Real preservation: {(nonch_real - nonch_lost) / nonch_real:.2%}" if nonch_real else "N/A")
    print(f"Noise removal: {nonch_removed / nonch_noise:.2%}" if nonch_noise else "N/A")

if __name__ == "__main__":
    filtered_file = "/project/ptgroup/Catherine/kTracker/run_C_module/filtered_output.txt"
    reduced_events = parse_filtered_output(filtered_file)

    original_file = "/project/ptgroup/Catherine/kTracker/data/small_raw/MC_JPsi_Pythia8_Target_April17_10000.root"
    noisy_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root"

    analyze_aggregated(reduced_events, original_file, noisy_file)
