# README: Running the `Fun4Sim.C` Module

This directory contains:

- `Fun4Sim.C`  → the main C macro (ROOT)
- `filter_script.py`  → the Python filter script to post-process the output

## How to run the C macro with ROOT

You can run the macro using ROOT and save the output to a text file:

```bash
root -b -q 'Fun4Sim.C(1000, "input.root", "output.root")' > output.txt 2>&1
```

Or run it normally using: 

```bash
root -b -q 'Fun4Sim.C(1000, "input.root", "output.root")'
```



## How to use the Python filter script

The Python filter script processes the ROOT output, extracting only the hits preserved after running the reducer, and writes them to `filtered_output.txt`.
You can run it with:

```bash
python3 filter_root_output.py --n_events 5000 --input_file newinput.root --output_file newout.root --filtered_file filtered_output.txt

```

## Analyze the Results

The script uses the `filtered_output` from `filter_root_output.py` from above.
Before running the analysis, make sure to update the input file paths in the main block of the script as needed.
Then run:

```bash
python3 ../analysis/analyze_c_reduction.py
```

## Compare the C++ Code Reduction and this Repo's Python Code Reduction

Again, be sure to adjust any file paths in the main block of the script to point to your data.
Then run:

```bash
python3 ../analysis/c_vs_python_reduction_comparison.py
```