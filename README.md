# kTracker Hit Reduction Pipeline

This project implements a modular pipeline for cleaning and reducing physics hits stored in ROOT files. It supports filters for out-of-time hit removal, declustering, and deduplication, along with tools for noise injection, visualization, and evaluation.

## 🛠️ Environment Setup

If you're using Rivanna or need custom environment setup:
```bash
source setup.sh
```

This sets up required environment variables and paths for ROOT.


## 📦 Dependencies

Once your environment is set up, install the necessary Python packages:

- [ROOT](https://root.cern/) (with PyROOT enabled)
- Python 3 (tested with Python 3.11)
- NumPy
- Uproot
- Pytest

To install the Python packages:

```bash
pip install numpy uproot pytest
```



## 🔊 Adding Noise for Testing

To inject noise into a ROOT file:

```bash
python3 noisy_data_gen/noisy_gen.py path/to/input.root
```

Example file to add noise to in Rivanna: ```/project/ptgroup/Catherine/kTracker/data/small_raw/small_combined.root```

This creates a noisy_output.root file with synthetic noise hits (electronic + clustered).


## 🚀 Running the Hit Reduction Pipeline

```bash
cd reduce_event
```


### Running the Main Pipeline Script

Run the main pipeline script:

```bash
python3 run_reduce_event.py
```


You can modify parameters inside the script to control:

- input_file and output_file
- Filters to apply: outoftime=True, decluster=True, dedup=True, etc.

Set the branches that are HIT vectors to be filtered when rewriting the file so all HIT vectors remain the same size

### 🧪 Running Tests

To run unit tests for the declustering logic:

```bash
pytest tests/
```


### 🔬 Analyzing Reduction Effectiveness

To compare original, noisy, and reduced files (if currently in reduce_event folder):

```bash
python3 ../scripts/analyze_reduction.py
```


This will report:

- Total real and noise hits
- Hits removed
- Real hit preservation and noise removal rates


## 📊 Plotting TDC Distributions

Plot per-detector TDC time histograms:

```bash
python3 scripts/tdctime_histogram/tdctime_histogram.py tdctime_histogram/yourfile.root
```


Fit a Gaussian to the TDC distribution:

```bash
python3 scripts/tdctime_histogram/gaussian_fit.py
```

