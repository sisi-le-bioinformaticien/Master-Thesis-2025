import pandas as pd
import matplotlib.pyplot as plt

csv_file = './affinity_tables_BRAF/affinity_pred_value.csv' 
gene_label = "BRAF" 

df = pd.read_csv(csv_file, delimiter=';')

row = df[df['Mutation'] == gene_label].iloc[0]

drugs = row.index[1:]  
ic50_values = row.values[1:]  

ic50_values = pd.to_numeric(ic50_values, errors='coerce')

plt.figure(figsize=(8, len(drugs) * 0.02)) 

colors = []
for val in ic50_values:
    if val < -1:
        colors.append('green')
    elif val > 1:
        colors.append('red')
    else:
        colors.append('gray')

plt.scatter(ic50_values, range(len(drugs)), color=colors)

for i, (drug, val, color) in enumerate(zip(drugs, ic50_values, colors)):
    if color in ['green', 'red']:
        plt.text(val, i + 0.2, drug, fontsize=8, ha='center', color=color)

plt.axvline(x=-1, color='green', linestyle='--', alpha=0.5)
plt.axvline(x=1, color='red', linestyle='--', alpha=0.5)

plt.xlim(-3, 2)
plt.ylim(-0.5, len(drugs) - 0.5)

plt.xlabel('IC50 Value')
plt.ylabel('Drug')
plt.title(f'IC50 Distribution for {gene_label}')

plt.yticks(range(len(drugs)), [''] * len(drugs))  

plt.grid(axis='x', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(f"{gene_label}_volcano_plot.png")
plt.show()
