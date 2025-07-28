import re

def parse_wire_lines(filepath):
    lines_by_key = {}
    with open(filepath, "r") as f:
        lines = f.readlines()

    current_key = None
    for line in lines:
        match_header = re.search(r'chamID: (\d+), paddle: (\d+), wire: (\d+)', line)
        if match_header:
            current_key = tuple(map(int, match_header.groups()))
        match_line = re.search(r'Wire line: \(([-\d.]+), ([-\d.]+)\) → \(([-\d.]+), ([-\d.]+)\)', line)
        if match_line and current_key:
            lines_by_key[current_key] = tuple(map(float, match_line.groups()))
            current_key = None  # reset after storing
    return lines_by_key

def compare_wire_lines(file_cpp, file_py, tolerance=1e-3):
    cpp_lines = parse_wire_lines(file_cpp)
    py_lines = parse_wire_lines(file_py)

    for key in sorted(set(cpp_lines) | set(py_lines)):
        cpp_vals = cpp_lines.get(key)
        py_vals = py_lines.get(key)

        if cpp_vals is None:
            print(f"[MISSING in C++] {key}")
            continue
        if py_vals is None:
            print(f"[MISSING in Python] {key}")
            continue

        diffs = [abs(c - p) for c, p in zip(cpp_vals, py_vals)]
        if any(d > tolerance for d in diffs):
            print(f"[DIFFERENT] {key}")
            print(f"  C++:    ({cpp_vals[0]:.3f}, {cpp_vals[1]:.3f}) → ({cpp_vals[2]:.3f}, {cpp_vals[3]:.3f})")
            print(f"  Python: ({py_vals[0]:.3f}, {py_vals[1]:.3f}) → ({py_vals[2]:.3f}, {py_vals[3]:.3f})")
        else:
            print(f"[MATCH] {key}")

# Example usage:
compare_wire_lines("/project/ptgroup/Catherine/kTracker/run_C_module/filtered_output.txt", "/project/ptgroup/Catherine/kTracker/reduce_event/output.txt")
