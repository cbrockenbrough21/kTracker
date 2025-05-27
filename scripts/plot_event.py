import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_comparison(file_paths, labels, event_number):
    num_files = len(file_paths)
    fig, axes = plt.subplots(1, num_files, figsize=(5 * num_files, 8), sharey=True)
    
    for idx, (file_name, label) in enumerate(zip(file_paths, labels)):
        with uproot.open(file_name) as file:
            tree = file["tree"]
            detector_id = tree["detectorID"].array(entry_start=event_number, entry_stop=event_number + 1)
            element_id = tree["elementID"].array(entry_start=event_number, entry_stop=event_number + 1)

            if len(detector_id) == 0 or len(element_id) == 0:
                print(f"No data found in {label} for event {event_number}")
                continue

            detector_id = detector_id[0].tolist()
            element_id = element_id[0].tolist()

            matrix = np.zeros((200, 60))
            for det_id, elem_id in zip(detector_id, element_id):
                det_index = det_id - 1
                elem_index = elem_id - 1
                if 0 <= det_index < 60 and 0 <= elem_index < 200:
                    matrix[elem_index, det_index] = 1

            ax = axes[idx] if num_files > 1 else axes
            ax.imshow(matrix, aspect="auto", cmap="Blues", origin="lower")
            ax.set_title(f"{label}\nEvent {event_number}")
            ax.set_xlabel("Detector ID")
            if idx == 0:
                ax.set_ylabel("Element ID")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_comparison(
        file_paths=[
            "/project/ptgroup/Catherine/kTracker/data/small_raw/small_combined.root",
            "/project/ptgroup/Catherine/kTracker/data/noisy/small_combined_varied_hits.root",
            "/project/ptgroup/Catherine/kTracker/data/cleaned/small_combined_cleaned.root"
        ],
        labels=[
            "Original",
            "Noisy",
            "Reduced"
        ],
        event_number=5  # Change this as needed
    )
