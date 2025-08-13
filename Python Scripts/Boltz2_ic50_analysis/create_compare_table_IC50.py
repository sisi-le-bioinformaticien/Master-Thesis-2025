import csv
import re

# Normalize drug names by removing non-alphanumeric characters and converting to lowercase
def normalize_drug_name(name):
    return re.sub(r'[^a-zA-Z0-9]', '', name).lower()

# Read GDSC_BRAF_IC50.csv and extract drug names and effect sizes
gdsc_data = {}
with open('affinity_tables_BRAF/GDSC_BRAF_IC50.csv', 'r', newline='') as gdsc_file:
    reader = csv.reader(gdsc_file, delimiter=',', quotechar='"', escapechar='\\')
    next(reader) 
    for row in reader:
        if not row:
            continue
        drug_name = row[0].strip('"')
        print(drug_name)
        try:
            effect_size = float(row[2])  # Effect size is in the third column
        except IndexError:
            print(f"Row skipped due to parsing error: {row}")
            continue
        normalized_name = normalize_drug_name(drug_name)
        gdsc_data[normalized_name] = effect_size

# Read affinity_pred_value.csv to extract drug names and Boltz predictions
boltz_data = {}
with open('affinity_tables_BRAF/affinity_pred_value.csv', 'r') as aff_file:
    reader = csv.reader(aff_file, delimiter=';')
    headers = next(reader) 
    values = next(reader)  
    drug_names = headers[1:]
    for i, drug in enumerate(drug_names):
        normalized_drug = normalize_drug_name(drug)
        boltz_value = float(values[i+1])  
        boltz_data[normalized_drug] = boltz_value

boltz_data_wt = {}
with open('affinity_tables_BRAF_WT/affinity_pred_value.csv', 'r') as aff_file:
    reader = csv.reader(aff_file, delimiter=';')
    headers = next(reader)  
    values = next(reader)   
    
    drug_names = headers[1:]
    for i, drug in enumerate(drug_names):
        normalized_drug = normalize_drug_name(drug)
        boltz_value = float(values[i+1]) 
        boltz_data_wt[normalized_drug] = boltz_value

common_drugs = []
for norm_drug in gdsc_data:
    if norm_drug in boltz_data and norm_drug in boltz_data_wt:
        common_drugs.append({
            'Drug': norm_drug,
            'GDSC_IC50': gdsc_data[norm_drug],
            'Boltz_IC50_BRAF': boltz_data[norm_drug],
            'Boltz_IC50_BRAF_WT': boltz_data_wt[norm_drug],
        })
    else:
        print(f"Skipping drug not found in all datasets: {norm_drug}")

with open('affinity_tables_BRAF_WT/drug_comparison_IC50.csv', 'w', newline='') as outfile:
    fieldnames = ['Drug', 'GDSC_IC50', 'Boltz_IC50_BRAF','Boltz_IC50_BRAF_WT']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for drug in common_drugs:
        writer.writerow(drug)

print(f"Generated drug_comparison_IC50.csv with {len(common_drugs)} common drugs.")
