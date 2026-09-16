"""
End-to-end D-peptide MPO analysis pipeline.
Usage: python run_pipeline.py --receptor SEQUENCE --peptides peptides.txt --output results/
"""
import argparse, os, subprocess, json, glob
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

try:
    from .affinity_utils import predicted_nominal_ic50_nm
except ImportError:  # Support direct execution: python examples/run_pipeline.py
    from affinity_utils import predicted_nominal_ic50_nm

def build_d_peptide_smiles(sequence):
    """Build SMILES for D-peptide from single-letter code (lowercase)."""
    aa_smiles = {
        'a': 'N[C@@H](C)C(=O)',           # D-Ala
        'r': 'N[C@@H](CCCNC(=N)N)C(=O)', # D-Arg
        'n': 'N[C@@H](CC(N)=O)C(=O)',     # D-Asn
        'd': 'N[C@@H](CC(=O)O)C(=O)',     # D-Asp
        'c': 'N[C@@H](CS)C(=O)',           # D-Cys
        'q': 'N[C@@H](CCC(N)=O)C(=O)',    # D-Gln
        'e': 'N[C@@H](CCC(=O)O)C(=O)',    # D-Glu
        'g': 'NCC(=O)',                    # Gly
        'h': 'N[C@@H](Cc1cnc[nH]1)C(=O)',# D-His
        'i': 'N[C@@H]([C@@H](C)CC)C(=O)',# D-Ile
        'l': 'N[C@@H](CC(C)C)C(=O)',      # D-Leu
        'k': 'N[C@@H](CCCCN)C(=O)',       # D-Lys
        'm': 'N[C@@H](CCSC)C(=O)',         # D-Met
        'f': 'N[C@@H](Cc1ccccc1)C(=O)',   # D-Phe
        'p': 'N1CCC[C@@H]1C(=O)',          # D-Pro
        's': 'N[C@@H](CO)C(=O)',           # D-Ser
        't': 'N[C@@H]([C@@H](O)C)C(=O)', # D-Thr
        'w': 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)', # D-Trp
        'y': 'N[C@@H](Cc1ccc(O)cc1)C(=O)',# D-Tyr
        'v': 'N[C@@H](C(C)C)C(=O)',        # D-Val
    }
    parts = []
    for i, aa in enumerate(sequence.lower()):
        if aa not in aa_smiles:
            raise ValueError(f"Unknown amino acid: {aa}")
        smi = aa_smiles[aa]
        if i == len(sequence)-1:
            # C-terminus: replace C(=O) with C(=O)O
            smi = smi.replace('C(=O)', 'C(=O)O', 1) if 'N1' not in smi else smi + 'O'
            # Remove trailing C(=O) if present
            if smi.endswith('C(=O)'):
                smi = smi[:-5] + 'C(=O)O'
        parts.append(smi)
    return ''.join(parts)

def get_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    return {
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'TPSA': round(Descriptors.TPSA(mol), 1),
        'MW':   round(Descriptors.MolWt(mol), 1),
        'HBD':  Descriptors.NumHDonors(mol),
        'HBA':  Descriptors.NumHAcceptors(mol),
    }

def create_boltz2_yaml(name, receptor_seq, ligand_smiles, output_dir):
    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {receptor_seq}
  - ligand:
      id: B
      smiles: "{ligand_smiles}"
properties:
  - affinity:
      binder: B
"""
    yaml_path = os.path.join(output_dir, f"{name}.yaml")
    with open(yaml_path, 'w') as f:
        f.write(yaml_content)
    return yaml_path

def run_boltz2(yaml_path, output_dir):
    cmd = [
        'boltz', 'predict', yaml_path,
        '--out_dir', output_dir,
        '--model', 'boltz2',
        '--affinity_checkpoint', os.path.expanduser('~/.boltz/boltz2_aff.ckpt'),
        '--sampling_steps', '200',
        '--diffusion_samples', '1',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0

def extract_ic50(result_dir):
    jsons = glob.glob(os.path.join(result_dir, '**/affinity*.json'), recursive=True)
    values = []
    for j in jsons:
        with open(j) as f:
            data = json.load(f)
        aff = data.get('affinity_pred_value')
        if aff is not None:
            values.append(predicted_nominal_ic50_nm(aff))
    return round(sum(values)/len(values), 1) if values else None

def mpo_score(LogP, TPSA, IC50_nM, lead_TPSA=263.3):
    c1 = 1 <= LogP <= 3
    c2 = TPSA < lead_TPSA
    c3 = IC50_nM < 1000
    logp_s = max(0, 1 - abs(LogP-2)/2)
    tpsa_s = max(0, 1 - TPSA/lead_TPSA)
    ic50_s = max(0, 1 - np.log10(max(IC50_nM,1))/np.log10(100000))
    score = 0.3*logp_s + 0.2*tpsa_s + 0.5*ic50_s
    return round(score, 4), c1 and c2 and c3

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--receptor', required=True, help='Receptor amino acid sequence')
    parser.add_argument('--peptides', required=True, help='Text file with one peptide sequence per line')
    parser.add_argument('--output', default='results/', help='Output directory')
    parser.add_argument('--lead_tpsa', type=float, default=263.3, help='Lead peptide TPSA threshold')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    results = []

    with open(args.peptides) as f:
        peptides = [l.strip() for l in f if l.strip()]

    print(f"Processing {len(peptides)} peptides...")
    for seq in peptides:
        print(f"  {seq}...", end=' ')
        try:
            smiles = build_d_peptide_smiles(seq)
            desc = get_descriptors(smiles)
            if desc is None:
                print("SMILES FAILED")
                continue

            yaml_path = create_boltz2_yaml(seq, args.receptor, smiles, args.output)
            boltz_dir = os.path.join(args.output, f'boltz_{seq}')
            success = run_boltz2(yaml_path, boltz_dir)

            ic50 = extract_ic50(boltz_dir) if success else None
            score, passes = mpo_score(desc['LogP'], desc['TPSA'],
                                       ic50 or 99999, args.lead_tpsa) if ic50 else (0, False)

            results.append({
                'sequence': seq,
                'SMILES': smiles,
                **desc,
                'IC50_nM': ic50,
                'MPO_score': score,
                'pass_all': passes,
            })
            print(f"IC50={ic50}nM MPO={score}")
        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results).sort_values('MPO_score', ascending=False)
    out_csv = os.path.join(args.output, 'mpo_results.csv')
    df.to_csv(out_csv, index=False)
    print(f"\nDone! Results: {out_csv}")
    print(f"Pass all MPO: {df['pass_all'].sum()}/{len(df)}")
