import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

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
            match = re.match(r"RunID:\s*\d+,\s*EventID:\s*(\d+)", line)
            if match:
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
        if current_event_id is not None:
            events[current_event_id] = pairs
    return events

def accumulate_stats(events, original_file, noisy_file):
    import uproot
    from collections import defaultdict

    tree_orig = uproot.open(original_file)["tree"]
    tree_noisy = uproot.open(noisy_file)["tree"]

    real_hits = defaultdict(int)
    noise_hits = defaultdict(int)
    final_hits = defaultdict(int)
    real_lost = defaultdict(int)
    noise_removed = defaultdict(int)

    # loop over requested events
    chunk_size = 1000
    for event_id, reduced_hits in events.items():
        found_entry = None
        # scan in chunks
        for chunk_idx, data in enumerate(tree_orig.iterate(["eventID"], step_size=chunk_size, library="np")):
            idx = np.where(data["eventID"] == event_id)[0]
            if len(idx) > 0:
                found_entry = chunk_idx * chunk_size + idx[0]
                break
        if found_entry is None:
            print(f"EventID {event_id} not found in raw!")
            continue

        # load this single entry
        orig_det = tree_orig["detectorID"].array(entry_start=found_entry, entry_stop=found_entry+1)[0]
        orig_elem = tree_orig["elementID"].array(entry_start=found_entry, entry_stop=found_entry+1)[0]
        noisy_det = tree_noisy["detectorID"].array(entry_start=found_entry, entry_stop=found_entry+1)[0]
        noisy_elem = tree_noisy["elementID"].array(entry_start=found_entry, entry_stop=found_entry+1)[0]

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

    return real_hits, noise_hits, final_hits, real_lost, noise_removed

def hits_from_reduced_tree(event_ids, reduced_tree):
    """
    For each eventID, collect its reduced (det,elem)
    """
    eventids = reduced_tree["eventID"].array()
    result = {}
    for event_id in event_ids:
        entries = np.where(eventids == event_id)[0]
        if len(entries)==0:
            continue
        entry = entries[0]
        d = reduced_tree["detectorID"].array(entry_start=entry, entry_stop=entry+1)[0]
        e = reduced_tree["elementID"].array(entry_start=entry, entry_stop=entry+1)[0]
        result[event_id] = set(zip(d,e))
    return result

def print_stats(stats, header, fileout):
    real_hits, noise_hits, final_hits, real_lost, noise_removed = stats
    chamber_ids = range(1, 31)
    nonchamber_ids = [d for d in real_hits if d > 30]

    fileout.write(f"\n=== {header} Per-Detector Stats (1–30) ===\n")
    fileout.write(f"{'Detector':>9} | {'Real':>6} | {'Noise':>7} | {'Final':>6} | {'Lost':>5} | {'Removed':>7} | {'Preserv%':>9} | {'NoiseRem%':>10}\n")
    fileout.write("-"*90 + "\n")

    for det in chamber_ids:
        r = real_hits[det]
        n = noise_hits[det]
        f = final_hits[det]
        lost = real_lost[det]
        rem = noise_removed[det]
        preserv = 100*(r - lost)/r if r else 100
        noiseeff = 100*rem/n if n else 0
        fileout.write(f"{det:9} | {r:6} | {n:7} | {f:6} | {lost:5} | {rem:7} | {preserv:8.2f}% | {noiseeff:9.2f}%\n")

    total_r = sum(real_hits[d] for d in chamber_ids)
    total_n = sum(noise_hits[d] for d in chamber_ids)
    total_f = sum(final_hits[d] for d in chamber_ids)
    total_lost = sum(real_lost[d] for d in chamber_ids)
    total_rem = sum(noise_removed[d] for d in chamber_ids)
    fileout.write(f"\n=== Chamber Totals ===\n")
    fileout.write(f"Total real: {total_r}, noise: {total_n}, final: {total_f}\n")
    fileout.write(f"Total lost: {total_lost}, removed: {total_rem}\n")
    fileout.write(
        f"Preservation: {( (total_r - total_lost) / total_r * 100):.2f}%\n"
        if total_r
        else "Preservation: N/A\n"
    )
    fileout.write(
        f"Noise removal: {(total_rem / total_n * 100):.2f}%\n"
        if total_n
        else "Noise removal: N/A\n"
    )

    total_rnch = sum(real_hits[d] for d in nonchamber_ids)
    total_nnch = sum(noise_hits[d] for d in nonchamber_ids)
    total_fnch = sum(final_hits[d] for d in nonchamber_ids)
    total_lostnch = sum(real_lost[d] for d in nonchamber_ids)
    total_remnch = sum(noise_removed[d] for d in nonchamber_ids)
    fileout.write(f"\n=== Non-Chamber Totals (>30) ===\n")
    fileout.write(f"Total real: {total_rnch}, noise: {total_nnch}, final: {total_fnch}\n")
    fileout.write(f"Total lost: {total_lostnch}, removed: {total_remnch}\n")
    fileout.write(
        f"Preservation: {( (total_rnch - total_lostnch) / total_rnch * 100):.2f}%\n"
        if total_rnch
        else "Preservation: N/A\n"
    )
    fileout.write(
        f"Noise removal: {(total_remnch / total_nnch * 100):.2f}%\n"
        if total_nnch
        else "Noise removal: N/A\n"
    )

def print_comparison(reduced, filtered, fileout):
    r_real, r_noise, r_final, r_lost, r_rem = reduced
    f_real, f_noise, f_final, f_lost, f_rem = filtered
    chamber_ids = range(1, 31)
    fileout.write(f"\n=== Side by Side Comparison (1–30) ===\n")
    fileout.write(f"{'Detector':>9} | {'Py Pres%':>10} | {'C++ Pres%':>10} | {'Py NoiseRem%':>13} | {'C++ NoiseRem%':>13}\n")
    fileout.write("-"*65 + "\n")
    for det in chamber_ids:
        r_r = r_real[det]
        n_r = r_noise[det]
        l_r = r_lost[det]
        rem_r = r_rem[det]
        pres_r = 100*(r_r - l_r)/r_r if r_r else 100
        nr_r = 100*rem_r/n_r if n_r else 0

        r_f = f_real[det]
        n_f = f_noise[det]
        l_f = f_lost[det]
        rem_f = f_rem[det]
        pres_f = 100*(r_f - l_f)/r_f if r_f else 100
        nr_f = 100*rem_f/n_f if n_f else 0

        fileout.write(f"{det:9} | {pres_r:10.2f} | {pres_f:10.2f} | {nr_r:13.2f} | {nr_f:13.2f}\n")

if __name__ == "__main__":
    filtered_file = "/project/ptgroup/Catherine/kTracker/run_C_module/filtered_hit_output.txt"
    filtered_events = parse_filtered_output(filtered_file)

    reduced_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/MC_JPsi_Pythia8_Target_April17_10000_onlyElectronic_cleaned.root" 
    orig_file = "/project/ptgroup/Catherine/kTracker/data/small_raw/MC_JPsi_Pythia8_Target_April17_10000.root"
    noisy_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root"

    reduced_tree = uproot.open(reduced_file)["tree"]

    # get event IDs
    event_ids = list(filtered_events.keys())

    # get reduced hits from python ROOT file
    reduced_events = hits_from_reduced_tree(event_ids, reduced_tree)

    # accumulate stats
    reduced_stats = accumulate_stats(reduced_events, orig_file, noisy_file)
    filtered_stats = accumulate_stats(filtered_events, orig_file, noisy_file)

    with open("comparison_results.txt","w") as f:
        # print each separately
        print_stats(reduced_stats,"Python Reduced",f)
        print_stats(filtered_stats,"C++ Reduced",f)
        # print comparison
        print_comparison(reduced_stats,filtered_stats,f)

    # also print to terminal
    with open("comparison_results.txt") as f:
        print(f.read())
