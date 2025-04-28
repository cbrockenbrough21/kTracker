import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Gaussian function to fit
def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

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
detector_id_focus = 32
mask = (detector_ids_flat == detector_id_focus)
tdc_filtered = tdc_times_flat[mask]

# Plot histogram and fit
counts, bin_edges = np.histogram(tdc_filtered, bins=100, range=(tdc_filtered.min(), tdc_filtered.max()))
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Initial guess for fit: (amplitude, mean, sigma)
p0 = [counts.max(), bin_centers[np.argmax(counts)], 10.0]

# Fit
popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0)

# Extract fitted parameters
amp, mean, sigma = popt
print(f"Fitted center (mean): {mean:.2f} ns")
print(f"Fitted width (sigma): {sigma:.2f} ns")

# Plot
plt.hist(tdc_filtered, bins=100, range=(tdc_filtered.min(), tdc_filtered.max()), alpha=0.6, label="TDC Data")
x_fit = np.linspace(tdc_filtered.min(), tdc_filtered.max(), 500)
y_fit = gaussian(x_fit, *popt)
plt.plot(x_fit, y_fit, 'r--', label="Gaussian Fit")
plt.xlabel("TDC Time (ns)")
plt.ylabel("Counts")
plt.title(f"TDC Time for Detector {detector_id_focus}")
plt.legend()
plt.grid(True)

# Save or show
plt.savefig(f"tdc_distribution_detector_{detector_id_focus}_gaussian_fit.png", dpi=300)
plt.close()
