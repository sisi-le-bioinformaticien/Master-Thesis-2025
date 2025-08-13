import os
import json
import csv

# USER INPUT
gene_name = "BRAF_WT"  # Example: "FLT3"

# Directories
input_dir = f'./{gene_name}_affinities'
output_dir = f'./affinity_tables_{gene_name}'
os.makedirs(output_dir, exist_ok=True)

# Initialize tables dynamically
affinity_keys = [
    "affinity_pred_value", "affinity_probability_binary"
]

tables = {key: {gene_name: {}} for key in affinity_keys}
detected_drugs = set()

# Process all JSON files
for filename in os.listdir(input_dir):
    if filename.endswith(".json"):
        filepath = os.path.join(input_dir, filename)
        with open(filepath, 'r') as f:
            data = json.load(f)

        # Extract drug name dynamically and fill table
        drug_name = filename.replace("affinity_template_", "").replace(".json", "")

        detected_drugs.add(drug_name)
        for key in affinity_keys:
            tables[key][gene_name][drug_name] = data[key]

# Sort detected drugs
detected_drugs = sorted(detected_drugs)

for table_name, table_data in tables.items():
    csv_path = os.path.join(output_dir, f"{table_name}.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(["Mutation"] + detected_drugs)

        row = [gene_name]
        for drug in detected_drugs:
            value = table_data[gene_name].get(drug, "NA")
            row.append(value)
        writer.writerow(row)

    print(f"Saved: {csv_path}")

print("All CSV tables have been created.")
