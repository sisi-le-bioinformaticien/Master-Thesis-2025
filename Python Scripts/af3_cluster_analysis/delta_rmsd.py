import csv

# Load WT RMSD values
wt_rmsd = {}
with open("wt_rmsd.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["chain"], int(row["resi"]))
        wt_rmsd[key] = float(row["rmsd"])

# Load MUT RMSD values
mut_rmsd = {}
with open("mut_rmsd.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["chain"], int(row["resi"]))
        mut_rmsd[key] = float(row["rmsd"])

# Compute ΔRMSD = MUT - WT
all_keys = set(wt_rmsd.keys()).union(mut_rmsd.keys())
delta_rmsd = {
    key: mut_rmsd.get(key, 0.0) - wt_rmsd.get(key, 0.0)
    for key in all_keys
}

# Save to CSV
with open("delta_rmsd.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["chain", "resi", "delta_rmsd"])
    for (chain, resi), val in sorted(delta_rmsd.items()):
        writer.writerow([chain, resi, val])


