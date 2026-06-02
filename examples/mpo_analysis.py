"""
Multi-objective optimization (MPO) scoring for D-peptide candidates.
Usage: python mpo_analysis.py --input candidates.csv --lead_tpsa 263.3 --output mpo_results.csv
"""
import argparse
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    return (round(Descriptors.MolLogP(mol), 2),
            round(Descriptors.TPSA(mol), 1),
            round(Descriptors.MolWt(mol), 1))

def mpo_score(LogP, TPSA, IC50_nM, lead_tpsa=263.3):
    """
    Composite MPO score (0-1, higher is better).
    Weights: LogP=0.3, TPSA=0.2, IC50=0.5
    """
    logp_s = max(0, 1 - abs(LogP - 2) / 2)           # peaks at LogP=2
    tpsa_s = max(0, 1 - TPSA / lead_tpsa)             # lower TPSA = better
    ic50_s = max(0, 1 - np.log10(max(IC50_nM, 1)) / np.log10(100000))
    return round(0.3*logp_s + 0.2*tpsa_s + 0.5*ic50_s, 4)

def run_mpo(df, lead_tpsa=263.3, smiles_col='SMILES', ic50_col='IC50_nM'):
    """Apply MPO criteria to a dataframe of candidates."""
    # Compute descriptors if needed
    if 'LogP' not in df.columns and smiles_col in df.columns:
        desc = df[smiles_col].apply(get_descriptors)
        df['LogP'] = [d[0] for d in desc]
        df['TPSA']  = [d[1] for d in desc]
        df['MW']    = [d[2] for d in desc]

    df = df.dropna(subset=['LogP', 'TPSA', ic50_col]).copy()

    # MPO criteria
    df['C1_LogP'] = df['LogP'].between(1, 3)
    df['C2_TPSA'] = df['TPSA'] < lead_tpsa
    df['C3_IC50'] = df[ic50_col] < 1000
    df['pass_all'] = df['C1_LogP'] & df['C2_TPSA'] & df['C3_IC50']
    df['MPO_score'] = df.apply(
        lambda r: mpo_score(r['LogP'], r['TPSA'], r[ic50_col], lead_tpsa), axis=1)

    return df.sort_values('MPO_score', ascending=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MPO analysis for D-peptide candidates')
    parser.add_argument('--input',      required=True, help='Input CSV with SMILES and IC50_nM columns')
    parser.add_argument('--lead_tpsa',  type=float, default=263.3, help='Lead peptide TPSA threshold')
    parser.add_argument('--output',     default='mpo_results.csv')
    parser.add_argument('--smiles_col', default='SMILES')
    parser.add_argument('--ic50_col',   default='IC50_nM')
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} candidates")

    results = run_mpo(df, args.lead_tpsa, args.smiles_col, args.ic50_col)
    results.to_csv(args.output, index=False)

    n_pass = results['pass_all'].sum()
    print(f"Pass all MPO criteria: {n_pass}/{len(results)}")
    print(f"Top 5:")
    cols = [c for c in ['sequence','IC50_nM','LogP','TPSA','MPO_score','pass_all'] if c in results.columns]
    print(results[cols].head(5).to_string(index=False))
    print(f"\nSaved: {args.output}")