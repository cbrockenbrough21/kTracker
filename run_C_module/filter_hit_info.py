import subprocess
import re
import argparse

def filter_hit_info(command, output_file="filtered_output.txt"):
    # Compile regex for speed
    run_event_re = re.compile(r"RunID:\s*(\d+),\s*EventID:\s*(\d+)")
    hit_dump_re = re.compile(r"^\s*\d+\s*:\s*\d+\s*:\s*\d+\s*:")

    with open(output_file, "w") as fout:
        # start the process
        proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True,
            bufsize=1
        )

        seen_run_event = False

        for line in proc.stdout:
            line = line.rstrip()
            if run_event_re.search(line):
                if seen_run_event:
                    continue  # skip duplicates
                seen_run_event = True
                fout.write(line + "\n")
            elif hit_dump_re.search(line):
                fout.write(line + "\n")
                seen_run_event = False

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
            "--filtered_file", type=str, default="filtered_hit_output.txt",
            help="Where to save the filtered text output (default: filtered_hit_output.txt)"
        )

        args = parser.parse_args()

        root_cmd = (
            f'root -b -q \'Fun4Sim.C({args.n_events}, "{args.input_file}", "{args.output_file}")\''
        )

        filter_hit_info(root_cmd, args.filtered_file)