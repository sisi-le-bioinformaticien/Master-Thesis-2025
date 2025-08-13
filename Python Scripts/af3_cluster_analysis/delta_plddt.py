# Load WT plddt values
wt_rmsd = {}
with open("wt_plddt.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["chain"], int(row["resi"]))
        wt_rmsd[key] = float(row["plddt"])

# Load MUT plddt values
mut_rmsd = {}
with open("mut_plddt.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["chain"], int(row["resi"]))
        mut_rmsd[key] = float(row["plddt"])

# Compute Δplddt = MUT - WT
all_keys = set(wt_rmsd.keys()).union(mut_rmsd.keys())
delta_rmsd = {
    key: -abs(mut_rmsd.get(key, 0.0) - wt_rmsd.get(key, 0.0))
    for key in all_keys
}

# Save to CSV
with open("delta_plddt.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["chain", "resi", "delta_plddt"])
    for (chain, resi), val in sorted(delta_rmsd.items()):
        writer.writerow([chain, resi, val])
