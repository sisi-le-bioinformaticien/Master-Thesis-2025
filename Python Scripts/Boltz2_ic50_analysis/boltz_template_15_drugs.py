import os
import re
import csv
import yaml
import json

with open('./smiles_15.json', 'r') as f:
    smiles_data = json.load(f)

output_dir = "./boltz_15_drug_templates"
os.makedirs(output_dir, exist_ok=True)

gene_name = "BRAF"
gene_sequence = "TVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELL"
gene_sequence_mutated = "TVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELL"

mutation_table = ["Wildtype", "p.V600E"]

log_rows = []

def safe_filename(name):
    return name.replace("*", "stop").replace(";", "_").replace(",", "_").replace(" ", "").replace(".", "")

def generate_yaml(sequence, name_suffix, drug_name, drug_smile):
    data = {
        'sequences': [
            {'protein': {'id': 'A', 'sequence': sequence}},
            {'ligand': {'id': 'B', 'smiles': drug_smile}},
        ],
        'properties': [
            {'affinity': {'binder': 'B'}}
        ],
    }

    filename = f"{gene_name}_{name_suffix}_{safe_filename(drug_name)}.yaml"
    filepath = os.path.join(output_dir, filename)
    with open(filepath, "w") as f:
        yaml.dump(data, f, sort_keys=False)

    print(f"Saved: {filepath}")
    return filename

for mutation_entry in mutation_table:
    if mutation_entry.lower() == "wildtype":
        sequence = gene_sequence
        suffix = "WT"
    else:
        sequence = gene_sequence_mutated
        suffix = "V600E"

    for drug_name, drug_smile in smiles_data.items():
        filename = generate_yaml(sequence, suffix, drug_name, drug_smile)
        log_rows.append([mutation_entry, suffix, "OK", filename])

log_path = os.path.join(output_dir, f"{gene_name}_mutation_log.csv")
with open(log_path, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["OriginalMutation", "MutationApplied", "Status", "OutputFile"])
    writer.writerows(log_rows)

print(f"CSV log saved: {log_path}")
