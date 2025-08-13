import csv
from pymol import cmd
from collections import defaultdict
import numpy as np
def compute_avg_plddt(prefix, n_models, out_csv):


    for i in range(1, n_models + 1):
        cmd.load(f"{prefix}{i}.cif", f"{prefix}{i}")

    # Collect pLDDT from each model
    plddt_map = defaultdict(list)

    for i in range(1, n_models + 1):
        model = cmd.get_model(f"{prefix}{i} and name CA")
        for atom in model.atom:
            key = (atom.chain, int(atom.resi))
            plddt_map[key].append(atom.b)

    # Average over models
    result = [(chain, resi, np.mean(plddts)) for (chain, resi), plddts in plddt_map.items()]

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chain", "resi", "plddt"])
        writer.writerows(result)