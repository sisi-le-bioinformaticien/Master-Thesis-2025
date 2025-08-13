import os
import pandas as pd
import matplotlib.pyplot as plt

# User-defined settings
gene_name = "BRAF"
tables_dir = f'./affinity_tables_{gene_name}'
output_dir = f'./affinity_graphs_{gene_name}'
os.makedirs(output_dir, exist_ok=True)

# Define table names and colors
table_files = [
    "affinity_pred_value.csv",
    "affinity_probability_binary.csv",
    "affinity_pred_value1.csv",
    "affinity_probability_binary1.csv",
    "affinity_pred_value2.csv",
    "affinity_probability_binary2.csv"
]

table_colors = {
    "affinity_pred_value": "red",
    "affinity_probability_binary": "red",
    "affinity_pred_value1": "green",
    "affinity_probability_binary1": "green",
    "affinity_pred_value2": "blue",
    "affinity_probability_binary2": "blue"
}

# Load all tables
tables = {}
for table_file in table_files:
    table_path = os.path.join(tables_dir, table_file)
    table_name = table_file.replace('.csv', '')
    tables[table_name] = pd.read_csv(table_path, delimiter=";")

# Extract list of drugs from columns
drugs = tables[table_name].columns[1:]  

# Group tables by type
binary_tables = ["affinity_probability_binary", "affinity_probability_binary1", "affinity_probability_binary2"]
pred_tables = ["affinity_pred_value", "affinity_pred_value1", "affinity_pred_value2"]

# Generate one graph per drug
for drug in drugs:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Top plot: Binary probabilities
    for table_name in binary_tables:
        df = tables[table_name]
        x = df['Mutation']
        y = df[drug]
        ax1.scatter(x, y, color=table_colors[table_name], label=table_name, s=100)

    ax1.set_title(f"{drug} - Binary Affinity Predictions ({gene_name})")
    ax1.set_ylabel("Binary Probability")
    ax1.set_ylim(0, 1)
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.5)

    # Bottom plot: Predicted values
    for table_name in pred_tables:
        df = tables[table_name]
        x = df['Mutation']
        y = df[drug]
        ax2.scatter(x, y, color=table_colors[table_name], label=table_name, s=100)

    ax2.set_title(f"{drug} - Predicted Affinity Values ({gene_name})")
    ax2.set_xlabel("Mutation Category")
    ax2.set_ylabel("Predicted Affinity")
    ax2.set_ylim(-3, 1)
    #ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.5)

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{gene_name}_{drug}_affinity_plot.png"))
    plt.close()

    print(f"Saved plot for {drug}")

print("All graphs have been created successfully.")
