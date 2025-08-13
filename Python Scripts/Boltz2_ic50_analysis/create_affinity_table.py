import os
import json
import csv

# USER INPUT
gene_name = "BRAF"  # Example: "FLT3"
mutations = ["WT", "V600E"]  # Example:  ["WT", "D835Y", "D835H", "D835E"] 
drugs = [
    "Dabrafenib", "PLX-4720", "SB590885", "Selumetinib", "Trametinib",
    "Foretinib", "Dasatinib", "Sapitinib", "Veliparib", "Dactolisib",
    "Alpelisib", "Taselisib", "Acetalax", "Osimertinib"
]
#drugs = ["Quizartinib"]

# Directories
input_dir = f'./{gene_name}_affinities'
output_dir = f'./affinity_tables_{gene_name}'
os.makedirs(output_dir, exist_ok=True)

# Initialize tables dynamically
affinity_keys = [
    "affinity_pred_value", "affinity_probability_binary",
    "affinity_pred_value1", "affinity_probability_binary1",
    "affinity_pred_value2", "affinity_probability_binary2"
]

tables = {key: {mutation: {} for mutation in mutations} for key in affinity_keys}

for filename in os.listdir(input_dir):
    if filename.endswith(".json"):
        filepath = os.path.join(input_dir, filename)
        with open(filepath, 'r') as f:
            data = json.load(f)
        mutation = None
        for mut in mutations:
            if mut in filename:
                mutation = mut
                break
        if mutation is None:
            print(f"Skipping {filename}: No matching mutation found.")
            continue
        drug_name = filename.replace("affinity_", "").replace(".json", "")
        drug_name = drug_name.replace(f"{gene_name}_{mutation}_", "")
        drug_name = drug_name.split("_")[-1]

        if drug_name not in drugs:
            print(f"Skipping {filename}: Drug {drug_name} not in drug list.")
            continue

        for key in affinity_keys:
            tables[key][mutation][drug_name] = data[key]

for table_name, table_data in tables.items():
    csv_path = os.path.join(output_dir, f"{table_name}.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(["Mutation"] + drugs)

        for mutation in mutations:
            row = [mutation]
            for drug in drugs:
                value = table_data[mutation].get(drug, "NA")
                row.append(value)
            writer.writerow(row)

    print(f"Saved: {csv_path}")

print("All CSV tables have been created.")
