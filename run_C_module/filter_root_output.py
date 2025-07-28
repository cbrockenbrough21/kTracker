import subprocess
import re
import argparse

def filter_root_output(command_str, output_file="filtered_output.txt"):
    # Define patterns to capture LUT block
    start_lut_re = re.compile(r"\[HODO LUT DEBUG\] Starting LUT init")
    end_lut_re = re.compile(r"\[HODO LUT DEBUG\] Finished LUT init")

    in_lut_block = False

    with open(output_file, "w") as fout:
        proc = subprocess.Popen(
            command_str,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True,
            bufsize=1
        )

        for line in proc.stdout:
            line = line.rstrip()
            print(line)

            # Start of LUT block
            if start_lut_re.search(line):
                in_lut_block = True
                fout.write("\n" + line + "\n")
                continue

            # End of LUT block
            if end_lut_re.search(line):
                fout.write(line + "\n")
                in_lut_block = False
                continue

            # Inside LUT block → write to file
            if in_lut_block:
                fout.write(line + "\n")

        proc.wait()
        if proc.returncode != 0:
            print(f"Command failed with code {proc.returncode}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Fun4Sim and filter the ROOT output.")
    parser.add_argument(
        "--n_events", type=int, default=1000,
        help="Number of events to run in Fun4Sim (default: 1000)"
    )
    parser.add_argument(
        "--input_file", type=str, default="/project/ptgroup/Catherine/kTracker/data/noisy/MC_JPsi_Pythia8_Target_April17_10000_noisy_onlyElectronic.root",
        help="Input ROOT file path for Fun4Sim (default: my_input.root)"
    )
    parser.add_argument(
        "--output_file", type=str, default="MC_JPsi_Pythia8_Target_April17_10000_onlyElectronic_cleaned.root",
        help="Output ROOT file path for Fun4Sim (default: cleaned_output.root)"
    )
    parser.add_argument(
        "--filtered_file", type=str, default="filtered_output.txt",
        help="Where to save the filtered text output (default: filtered_output.txt)"
    )

    args = parser.parse_args()

    root_cmd = (
        f'root -b -q \'Fun4Sim.C({args.n_events}, "{args.input_file}", "{args.output_file}")\''
    )

    filter_root_output(root_cmd, args.filtered_file)
