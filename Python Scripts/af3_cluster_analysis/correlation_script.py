import csv
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress

# Load delta_rmsd.csv
delta_rmsd = {}
with open("delta_rmsd.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["chain"], int(row["resi"]))
        delta_rmsd[key] = float(row["delta_rmsd"])

# Load delta_plddt.csv
delta_plddt = {}
with open("delta_plddt.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row["chain"], int(row["resi"]))
        delta_plddt[key] = float(row["delta_plddt"])

# Combine data
common_keys = sorted(set(delta_rmsd) & set(delta_plddt))
x_vals = [delta_plddt[k] for k in common_keys]  # X = ΔpLDDT
y_vals = [delta_rmsd[k] for k in common_keys]   # Y = ΔRMSD
labels = [f"{k[0]}{k[1]}" for k in common_keys]
chains = [k[0] for k in common_keys]

# Define chain colors
chain_colors = {"B": "blue", "M": "orange"}
colors = [chain_colors.get(chain, "gray") for chain in chains]

# Compute Spearman correlation
rho, pval = spearmanr(x_vals, y_vals)

# Compute linear regression for trendline
slope, intercept, r_value, p_value, std_err = linregress(x_vals, y_vals)
trend_y = [slope * x + intercept for x in x_vals]

# Plot
plt.figure(figsize=(7, 5))
scatter = plt.scatter(x_vals, y_vals, c=colors, edgecolors="k", alpha=0.75, label="Residues")
plt.plot(x_vals, trend_y, color="crimson", linewidth=2, label="Linear trend")

plt.xlabel("ΔpLDDT (mut - wt)", fontsize=12)
plt.ylabel("ΔRMSD (mut - wt)", fontsize=12)
plt.title(f"ΔRMSD vs ΔpLDDT (Spearman ρ = {rho:.2f}, p = {pval:.2g})", fontsize=13)
plt.grid(True, linestyle="--", alpha=0.5)

for i, label in enumerate(labels):
    if abs(x_vals[i]) > 10 or abs(y_vals[i]) > 1:
        plt.text(x_vals[i], y_vals[i], label, fontsize=6, ha='right')
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Chain B', markerfacecolor='blue', markersize=8, markeredgecolor='k'),
    Line2D([0], [0], marker='o', color='w', label='Chain M', markerfacecolor='orange', markersize=8, markeredgecolor='k'),
    Line2D([0], [0], color='crimson', lw=2, label='Linear trend')
]
plt.legend(handles=legend_elements, fontsize=9)

plt.tight_layout()
plt.savefig("spearman_by_chain_colored.png", dpi=300)
plt.show()
