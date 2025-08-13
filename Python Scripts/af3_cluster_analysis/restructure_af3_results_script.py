import shutil, os

source_dir = "./brafmut_dabrafenib_domain"
target_dir = "./formatted_brafmut_dabrafenib_domain"
os.makedirs(target_dir, exist_ok=True)

count = 1
for i in range(1, 11):
    for j in range(2):
        folder = f"seed-{i}_sample-{j}"
        source_pdb = os.path.join(source_dir, folder, "model.cif") 
        target_pdb = os.path.join(target_dir, f"mut{count}.cif")
        shutil.copy(source_pdb, target_pdb)
        count += 1