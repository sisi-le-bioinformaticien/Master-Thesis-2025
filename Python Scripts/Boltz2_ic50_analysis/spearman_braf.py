
import argparse
import pandas as pd
import re
from scipy.stats import spearmanr

def normalize_mut(m):
    if pd.isna(m):
        return None
    m = str(m).strip()
    m = m.replace('*', '*') 
    if m.startswith('p.'):
        m = 'p' + m[2:]
    elif m.startswith('p'):
        m = m
    else:
        m = 'p' + m
    m = m.replace('.', '')
    m = re.sub(r'\s+', '', m)
    return m

def parse_mut_from_filename(fn):
    m = re.search(r'affinity_BRAF_(p[^_]+)_', fn)
    if m:
        return m.group(1)
    m = re.search(r'(p[A-Za-z0-9_\*]+)', fn)
    return m.group(1) if m else None

def main():
    ap = argparse.ArgumentParser(description="Compute Spearman correlation across common BRAF mutations.")
    ap.add_argument("ic50_csv", help="Path to CSV with columns: BRAF_mutations, avg_IC50 (others ignored)")
    ap.add_argument("affinity_csv", help="Path to CSV with columns: filename, affinity_pred_value")
    ap.add_argument("-o", "--out", default="common_mutations_merged.csv", help="Output CSV for merged data")
    args = ap.parse_args()

    ic50 = pd.read_csv(args.ic50_csv)
    aff = pd.read_csv(args.affinity_csv)

    ic50['mut_key'] = ic50['BRAF_mutations'].apply(normalize_mut)

    # Parse mutation from filename and normalize
    aff['parsed_mut'] = aff['filename'].apply(parse_mut_from_filename)
    aff['mut_key'] = aff['parsed_mut'].apply(normalize_mut)
    merged = pd.merge(
        ic50[['BRAF_mutations', 'avg_IC50', 'n_cell_lines', 'mut_key']],
        aff[['filename', 'affinity_pred_value', 'mut_key']],
        on='mut_key',
        how='inner',
        validate='m:1' 
    )

    # Compute Spearman correlation between avg_IC50 and affinity_pred_value
    rho, pval = spearmanr(merged['avg_IC50'], merged['affinity_pred_value'])

    # Save merged and print stats
    merged = merged.sort_values('BRAF_mutations').reset_index(drop=True)
    merged.to_csv(args.out, index=False)

    print(f"Common mutations: {len(merged)}")
    print(f"Spearman rho: {rho:.4f}")
    print(f"P-value: {pval:.4g}")
    print()
    print("Head of merged data:")
    print(merged.head(20).to_string(index=False))

if __name__ == "__main__":
    main()
