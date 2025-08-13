import csv
import matplotlib.pyplot as plt
from collections import defaultdict

chain_data = defaultdict(lambda: {"residues": [], "plddt": []})

with open("delta_plddt.csv") as f:
    reader = csv.reader(f)
    for row in reader:
        chain, resi, rmsd = row
        if chain == "chain":
            continue
        chain_data[chain]["residues"].append(f"{chain}{resi}")
        chain_data[chain]["plddt"].append(float(rmsd))

plt.figure(figsize=(14, 5))
colors = plt.get_cmap("tab10")

for i, (chain, data) in enumerate(sorted(chain_data.items())):
    plt.plot(data["residues"], data["plddt"], marker='o', linestyle='-', label=f"Chain {chain}", color=colors(i))

plt.xticks(rotation=90, fontsize=6)
plt.yticks(fontsize=8)
plt.xlabel("Residue", fontsize=10)
plt.ylabel("Δplddt", fontsize=10)
plt.title("Per-Residue plddt Difference", fontsize=12)
plt.legend(fontsize=8)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("delta_plddt_plot_colored.png", dpi=300)
plt.show()
