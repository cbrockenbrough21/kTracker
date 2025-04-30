# kTracker Hit Reduction Pipeline

This project implements a modular pipeline for cleaning and reducing physics hits stored in ROOT files. It supports filters for out-of-time hit removal, declustering, and deduplication, along with tools for noise injection, visualization, and evaluation.

## 📦 Dependencies

Before running the pipeline, ensure the following dependencies are installed:

- [ROOT](https://root.cern/) (with PyROOT enabled)
- Python 3 (tested with Python 3.11)
- NumPy
- Uproot
- Pytest

To install the Python packages:

```bash
pip install numpy uproot pytest


## 🚀 Running the Hit Reduction Pipeline

Run the main pipeline script:

python3 reduce_event/run_reduce_event.py

You can modify parameters inside the script to control:

- input_file and output_file
- Filters to apply: outoftime=True, decluster=True, dedup=True, etc.
- TDC timing cuts: tdc_center, tdc_width

Set the branches that are HIT vectors to be filtered when rewriting the file so all HIT vectors remain the same size

## 🧪 Running Tests

To run unit tests for the declustering logic:

pytest reduce_event/tests/

## 🔬 Analyzing Reduction Effectiveness

To compare original, noisy, and reduced files:

python3 reduce_event/scripts/analyze_reduction.py

This will report:

- Total real and noise hits
- Hits removed
- Real hit preservation and noise removal rates

## 🔊 Adding Noise for Testing

To inject noise into a ROOT file:

python3 noisy_data_gen/noisy_gen.py path/to/input.root

This creates a noisy_output.root file with synthetic noise hits (electronic + clustered).

## 📊 Plotting TDC Distributions

Plot per-detector TDC time histograms:

python3 tdctime_histogram/tdctime_histogram.py tdctime_histogram/yourfile.root

Fit a Gaussian to the TDC distribution:

python3 tdctime_histogram/gaussian_fit.py

