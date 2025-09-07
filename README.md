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

To install the Python packages:

```bash
pip install numpy uproot
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
- Filters to apply: outoftime=True, decluster=True, dedup=True, hodomask=True etc.

Set the branches that are HIT vectors to be filtered when rewriting the file so all HIT vectors remain the same size



### 🔬 Analyzing Reduction Effectiveness

To compare original, noisy, and reduced files (if currently in reduce_event folder):

```bash
python3 ../analysis/analyze_python_reduction.py
```


This will report:

- Total real and noise hits
- Hits removed
- Real hit preservation and noise removal rates


## Running `Fun4Sim.C` Module

```bash
cd fun4sim
```

This directory contains:

- `Fun4Sim.C`  → the main C macro (ROOT)
- `filter_hit_info.py`  → the Python filter script to post-process the output

### How to Run the C Macro with ROOT

You can run the macro:

```bash
root -b -q 'Fun4Sim.C(1000, "input.root", "output.root")'
```

Reducer Options
The event reducer options can be set in the line:

`reco->set_evt_reducer_opt("h"); `

You can change the string passed to set_evt_reducer_opt to enable different combinations of reduction options. The following flags are supported:

a (afterhit), h (hodomask), o (outoftime), c (decluster), m (mergehodo), t (triggermask), s (sagitta), g (hough), r (realization), n (difnim).

For example, to enable both the hodomask and decluster reducers, you would set:

`reco->set_evt_reducer_opt("hc");`

### How to Use the Python Filter Script

The Python filter script processes the ROOT output, extracting only the hits preserved after running the reducer, and writes them to `filtered_hit_output.txt`.
You can run it with:

```bash
python3 filter_hit_info.py --n_events 5000 --input_file newinput.root --output_file newout.root --filtered_file filtered_hit_output.txt

```

### Analyze the Results

The script uses the `filtered_hit_output.txt` from `filter_root_output.py` from above.
Before running the analysis, make sure to update the input file paths in the main block of the script as needed.
Then run:

```bash
python3 ../analysis/analyze_c_reduction.py
```

### Compare the C++ Code Reduction and this Repo's Python Code Reduction

Again, be sure to adjust any file paths in the main block of the script to point to your data.
Then run:

```bash
python3 ../analysis/analyze_c_vs_python_reduction.py
```

