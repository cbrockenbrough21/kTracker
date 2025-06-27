import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from collections import Counter

def get_event_hits(file_path, event_number):
    """Return a set of (detectorID, elementID) tuples for a specific event."""
    with uproot.open(file_path) as file:
        tree = file["tree"]
        det_ids = tree["detectorID"].array(entry_start=event_number, entry_stop=event_number + 1)
        elem_ids = tree["elementID"].array(entry_start=event_number, entry_stop=event_number + 1)

        if len(det_ids) == 0 or len(elem_ids) == 0:
            return set()

        return set(zip(det_ids[0].tolist(), elem_ids[0].tolist()))

def plot_event_difference(original_file, noisy_file, reduced_file, event_number):
    # Load hits
    original_hits = get_event_hits(original_file, event_number)
    noisy_hits = get_event_hits(noisy_file, event_number)
    reduced_hits = get_event_hits(reduced_file, event_number)

    # Classify
    preserved_real = original_hits & reduced_hits
    removed_real = original_hits - reduced_hits
    injected_noise = noisy_hits - original_hits
    preserved_noise = injected_noise & reduced_hits
    removed_noise = injected_noise - reduced_hits

    # Matrix for plotting: 0 = none, 1-4 = classes
    matrix = np.zeros((200, 60))

    for det, elem in preserved_real:
        if 0 <= det < 60 and 0 <= elem < 200:
            matrix[elem - 1, det - 1] = 1
    for det, elem in removed_real:
        if 0 <= det < 60 and 0 <= elem < 200:
            matrix[elem - 1, det - 1] = 2
    for det, elem in preserved_noise:
        if 0 <= det < 60 and 0 <= elem < 200:
            matrix[elem - 1, det - 1] = 3
    for det, elem in removed_noise:
        if 0 <= det < 60 and 0 <= elem < 200:
            matrix[elem - 1, det - 1] = 4

    # Colors: white, black, red, yellow, green
    cmap = ListedColormap([
        "white",     # 0: no hit
        "#000000",   # 1: preserved real hit
        "#d62728",   # 2: removed real hit
        "#ffdf00",   # 3: preserved noise
        "#2ca02c",   # 4: removed noise
    ])

    legend = [
        mpatches.Patch(color="#000000", label="Preserved real hit"),
        mpatches.Patch(color="#d62728", label="Removed real hit"),
        mpatches.Patch(color="#ffdf00", label="Preserved noise"),
        mpatches.Patch(color="#2ca02c", label="Removed noise"),
    ]

    plt.figure(figsize=(8, 8))
    plt.imshow(matrix, aspect="auto", origin="lower", cmap=cmap)
    plt.title(f"Event {event_number}: Hit Matrix")
    plt.xlabel("Detector ID")
    plt.ylabel("Element ID")

    plt.legend(handles=legend, loc="upper right")
    plt.tight_layout()
    plt.savefig(f"event_{event_number}_difference.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    plot_event_difference(
        original_file = "/project/ptgroup/Catherine/kTracker/data/small_raw/MC_JPsi_Pythia8_Target_April17_10000.root",
        noisy_file = "/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root",
        reduced_file = "/project/ptgroup/Catherine/kTracker/data/cleaned/MC_JPsi_Pythia8_Target_April17_10000_onlyElectronic_cleaned.root", 
        event_number=20
    )
