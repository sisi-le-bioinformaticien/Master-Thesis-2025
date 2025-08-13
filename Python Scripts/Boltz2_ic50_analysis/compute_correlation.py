import csv
from scipy.stats import spearmanr
import seaborn as sbs
import matplotlib.pyplot as plt

gdsc_values = []
boltz_values = []
boltz_values_wt = []


with open('affinity_tables_BRAF_WT/drug_comparison_IC50.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        gdsc_values.append(float(row['GDSC_IC50']))
        boltz_values.append(float(row['Boltz_IC50_BRAF']))
        boltz_values_wt.append(float(row['Boltz_IC50_BRAF_WT']))
    


# Compute Pearson correlation coefficient
correlation, p_value = spearmanr(gdsc_values, boltz_values)

fig, f = plt.subplots()
sbs.regplot(y=gdsc_values, x=boltz_values)
f.set_title("Spearman Correlation between GDSC and boltz (BRAF V600E) log(IC50) values")
f.set(xlabel='Boltz log(IC50) prediction', ylabel='GDSC log(IC50)')
f.legend()
f.get_figure().savefig('corr_gdsc_brafmut.png')
print(f"Spearman correlation between GDSC and boltz (BRAF V600E) log(IC50) values: {correlation:.4f}")
print(f"P-value: {p_value:.4e}")

correlation, p_value = spearmanr(gdsc_values, boltz_values_wt)


fig1, f1 = plt.subplots()
sbs.regplot(y=gdsc_values, x=boltz_values_wt)
f1.set_title("Spearman Correlation between GDSC and boltz (BRAF wildtype) log(IC50) values")
f1.set(xlabel='Boltz log(IC50) prediction', ylabel='GDSC log(IC50)')
f1.legend()
f1.get_figure().savefig('corr_gdsc_brafwt.png')
print(f"Spearman correlation between GDSC and boltz (BRAF wildtype) log(IC50) values: {correlation:.4f}")
print(f"P-value: {p_value:.4e}")


correlation, p_value = spearmanr(boltz_values, boltz_values_wt)

fig2, f2 = plt.subplots()
sbs.regplot(y=boltz_values, x=boltz_values_wt)
f2.set_title("Spearman Correlation between boltz (BRAF V600E)  and boltz (BRAF WT)  log(IC50) values")
f2.set(xlabel='Boltz_WT log(IC50) prediction', ylabel='Boltz_mut log(IC50) prediction')
f2.legend()
f2.get_figure().savefig('corr_brafmut_brafwt.png')
print(f"Spearman correlation between boltz (BRAF V600E)  and boltz (BRAF WT)  log(IC50) values: {correlation:.4f}")
print(f"P-value: {p_value:.4e}")