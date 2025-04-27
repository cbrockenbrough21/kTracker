import uproot
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file
file = uproot.open("run_002281_spill_000000000_spin_vector.root")
tree = file["tree"]

# Read arrays
detector_ids = tree["detectorIDs"].array(library="np")
tdc_times = tree["tdcTimes"].array(library="np")

# Flatten
detector_ids_flat = np.concatenate(detector_ids)
tdc_times_flat = np.concatenate(tdc_times)

# Focus on a specific detector
detector_id_focus = 6
mask = (detector_ids_flat == detector_id_focus)
tdc_filtered = tdc_times_flat[mask]

unique_detector_ids = np.unique(detector_ids_flat)
print("Unique detector IDs in file:", unique_detector_ids)

# Get min and max of TDCs
tdc_min = tdc_filtered.min()
tdc_max = tdc_filtered.max()

print(f"TDC range for detector {detector_id_focus}: {tdc_min:.2f} to {tdc_max:.2f} ns")

# Now plot, automatically using (min, max)
plt.hist(tdc_filtered, bins=100, range=(tdc_min, tdc_max))
plt.xlabel("TDC Time (ns)")
plt.ylabel("Counts")
plt.title(f"TDC Time Distribution for Detector ID {detector_id_focus}")
plt.grid(True)
plt.savefig(f"tdc_distribution_detector_{detector_id_focus}.png", dpi=300)
plt.close()
