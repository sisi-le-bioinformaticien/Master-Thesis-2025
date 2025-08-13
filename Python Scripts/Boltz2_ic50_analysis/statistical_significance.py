import os
import pandas as pd
from scipy.stats import wilcoxon

# Settings
gene_name = "FLT3"
#mutations = ["WT", "V600E"]
mutations =  ["WT", "D835Y", "D835H", "D835E"]
tables_dir = f'./affinity_tables_{gene_name}'
output_dir = f'./statistical_tests_RESULTS'
os.makedirs(output_dir, exist_ok=True)

# Replicate groups
replicate_groups = {
    "affinity_pred_value": ["affinity_pred_value.csv", "affinity_pred_value1.csv", "affinity_pred_value2.csv"],
    "affinity_probability_binary": ["affinity_probability_binary.csv", "affinity_probability_binary1.csv", "affinity_probability_binary2.csv"]
}

results = []

for group_name, file_list in replicate_groups.items():
    dfs = [pd.read_csv(os.path.join(tables_dir, file), sep=';') for file in file_list]

    drugs = dfs[0].columns[1:] 

    for drug in drugs:
        wt_values = []
        for df in dfs:
            wt_value = df[df['Mutation'] == "WT"][drug].values[0]
            wt_values.append(wt_value)

        for mutation in mutations:
            if mutation == "WT":
                continue  

            mut_values = []
            for df in dfs:
                mut_value = df[df['Mutation'] == mutation][drug].values[0]
                mut_values.append(mut_value)

            try:
                stat, p_value = wilcoxon(wt_values, mut_values)

                results.append({
                    "Table": group_name,
                    "Drug": drug,
                    "Mutation": mutation,
                    "WT_Mean": sum(wt_values) / len(wt_values),
                    "Mutation_Mean": sum(mut_values) / len(mut_values),
                    "Difference_Mean": (sum(mut_values) / len(mut_values)) - (sum(wt_values) / len(wt_values)),
                    "Wilcoxon_p_value": p_value
                })

            except Exception as e:
                print(f"Error processing {group_name} - {drug} - {mutation}: {e}")

# Save results
results_df = pd.DataFrame(results)
output_path = os.path.join(output_dir, f"{gene_name}_wilcoxon_results.csv")
results_df.to_csv(output_path, index=False, sep=';')

print(f"Wilcoxon test results saved to: {output_path}")