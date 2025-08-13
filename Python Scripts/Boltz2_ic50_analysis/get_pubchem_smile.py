import requests
import csv 

SEPARATOR = ";"


def get_smiles(compound_name):
    # PubChem API URL
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/SMILES/JSON"
    
    # Make the request to the PubChem API
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['SMILES']
        return smiles
    else:
        return None



def get_drugs_name(file):
    drugs_names = {}
    with open(file, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        spamreader.__next__() # dirty way to skip header
        for row in spamreader:
            drugs_names[row[1].strip()]= [row[1].strip()] +  [d.strip() for d in row[2].split(',')]
    
    return drugs_names
         



drugs_smiles = { 'drug_name' : [], 'drug_smile': []}
drugs_names = get_drugs_name("./Drug_listFri Jun 27 14_00_22 2025.csv")
drugs_id = []


for drug_id, drug_names in drugs_names.items():
    for drug_name in drug_names: 
        smile = get_smiles(drug_name)
        if smile is not None:
            drugs_smiles['drug_name'].append(drug_id)
            drugs_smiles['drug_smile'].append(smile)
            break  


print(drugs_smiles)


with open("mycsvfile.csv", "w", newline="") as f:
    f.write(SEPARATOR.join(drugs_smiles.keys()))
    f.write('\n')
    for drug_name, drug_smile in zip(drugs_smiles['drug_name'], drugs_smiles['drug_smile']):
        f.write(f"{drug_name}{SEPARATOR}{drug_smile}")
        f.write('\n')
