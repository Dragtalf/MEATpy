import sys
import csv

if len(sys.argv) != 3:
    print("Usage: fix_headers_helper.py <file_path> <log_path>")
    sys.exit(1)

filepath = sys.argv[1]
log_path = sys.argv[2]

with open(log_path) as f:
    reader = csv.DictReader(f)
    rename_map = {
        row['original_filename'].replace('.fasta', ''): row['new_filename'].replace('.fasta', '')
        for row in reader
    }

with open(filepath) as f:
    lines = f.readlines()

with open(filepath, "w") as f:
    for line in lines:
        if line.startswith(">"):
            for old, new in rename_map.items():
                if line.startswith(f">{old}~~"):
                    line = line.replace(f">{old}~~", f">{new}~~")
        f.write(line)