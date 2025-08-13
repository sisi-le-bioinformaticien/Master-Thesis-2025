import os
import re
import csv
import json

output_dir = "./af3_templates"
os.makedirs(output_dir, exist_ok=True)

gene_name = "BRAF"
gene_sequence = "TVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELL"
drug_name = "dabrafenib"
drug_smile = "CC(C)(C)C1=NC(=C(S1)C2=NC(=NC=C2)N)C3=C(C(=CC=C3)NS(=O)(=O)C4=C(C=CC=C4F)F)F"
offset = 458

mutation_table = [
    "Wildtype", "p.V600E", "p.D594G", "p.K601N", "p.G464E", "p.G466V", "p.G469A", "p.V600M",
    "p.A712T", "p.G464V", "p.G469V", "p.G596R", "p.K499T", 
    "p.L597R", "p.L597V", "p.M650I", "p.N581K", "p.N581Y", "p.R506G", "p.R671Q", 
    "p.S683N", "p.T529A", "p.V645I"
]

log_rows = []

def apply_mutations(seq, mutations, offset):
    sequence = seq
    applied = []
    status = "OK"
    for mut in mutations:
        mut = mut.strip()
        m = re.match(r"p\.([A-Z])(\d+)([A-Z*])", mut)
        if m:
            ref, pos, alt = m.groups()
            i = int(pos) - offset
            if i < 0 or i >= len(sequence):
                status = f"out-of-range: {mut}"
                continue
            if alt == "*":
                sequence = sequence[:i] + "*"
                sequence = sequence[:i+1]
                applied.append(f"{mut} (stop)")
            else:
                sequence = sequence[:i] + alt + sequence[i+1:]
                applied.append(mut)
            continue
        status = f"unrecognized: {mut}"
    return sequence, applied, status

def safe_filename(name):
    return name.replace("*", "stop").replace(";", "_").replace(",", "_").replace(" ", "").replace(".", "")

def generate_af3_json(sequence, name_suffix):
    data = {
        "name": f"{gene_name}_{name_suffix}_{drug_name}",
        "modelSeeds": [10, 20, 40],
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence}},
            {"ligand": {"id": "B", "smiles": drug_smile}},
        ],
        "dialect": "alphafold3",
        "version": 1
    }

    filename = f"{gene_name}_{name_suffix}_{drug_name}.json"
    filepath = os.path.join(output_dir, filename)
    with open(filepath, "w") as f:
        json.dump(data, f, indent=2)

    print(f"Saved: {filepath}")
    return filename

for mutation_entry in mutation_table:
    if mutation_entry.lower() == "wildtype":
        filename = generate_af3_json(gene_sequence, "WT")
        log_rows.append(["Wildtype", "", "OK", filename])
    else:
        raw_mutations = re.split(r"[;,]", mutation_entry)
        mutated_seq, applied_list, status = apply_mutations(gene_sequence, raw_mutations, offset)
        suffix = safe_filename("_".join(raw_mutations))
        filename = generate_af3_json(mutated_seq, suffix)
        log_rows.append([mutation_entry, "; ".join(applied_list), status, filename])

# Write to CSV
log_path = os.path.join(output_dir, f"{gene_name}_af3_mutation_log.csv")
with open(log_path, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["OriginalMutation", "AppliedMutations", "Status", "OutputFile"])
    writer.writerows(log_rows)

