import numpy as np, csv
from pymol import cmd

def compute_residue_rmsd(prefix, n_models, out_csv):
    for i in range(1, n_models + 1):
        cmd.load(f"{prefix}{i}.cif", f"{prefix}{i}")

    # Align all models to the first
    for i in range(2, n_models + 1):
        cmd.align(f"{prefix}{i}", f"{prefix}1")

    # Collect all unique (chain, resi) from the reference model
    seen = set()
    residues = []
    model = cmd.get_model(f"{prefix}1 and name CA")
    for atom in model.atom:
        key = (atom.chain, int(atom.resi))
        if key not in seen:
            seen.add(key)
            residues.append(key)  # (chain, resi)

    # Compute RMSD per (chain, resi)
    res_rmsd = []
    for chain, resi in residues:
        coords = []
        for i in range(1, n_models + 1):
            try:
                coord = cmd.get_atom_coords(f"{prefix}{i} and name CA and resi {resi} and chain {chain}")
                coords.append(coord)
            except:
                coords.append([float('nan')] * 3)
        coords = np.array(coords)
        coords = coords[~np.isnan(coords).any(axis=1)]
        if len(coords) > 1:
            std = np.std(coords, axis=0)
            rmsd = np.linalg.norm(std)
        else:
            rmsd = 0.0
        res_rmsd.append((chain, resi, rmsd))

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chain", "resi", "rmsd"])
        writer.writerows(res_rmsd)

        