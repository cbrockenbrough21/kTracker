import subprocess
import re

def filter_root_output(command, output_file="filtered_output.txt"):
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
    # edit to your ROOT command:
    root_cmd = 'root -b -q Fun4Sim.C'
    filter_root_output(root_cmd, "filtered_output.txt")
