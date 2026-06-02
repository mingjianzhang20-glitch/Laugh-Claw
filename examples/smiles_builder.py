"""
Build D-peptide SMILES from single-letter sequence.
Usage: python smiles_builder.py --seq lgrmg
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
import argparse

# D-amino acid SMILES fragments (using @@ stereochemistry)
D_AA = {
    'a': ('N[C@@H](C)C(=O)',              'D-Ala'),
    'r': ('N[C@@H](CCCNC(=N)N)C(=O)',    'D-Arg'),
    'n': ('N[C@@H](CC(N)=O)C(=O)',        'D-Asn'),
    'd': ('N[C@@H](CC(=O)O)C(=O)',        'D-Asp'),
    'c': ('N[C@@H](CS)C(=O)',             'D-Cys'),
    'q': ('N[C@@H](CCC(N)=O)C(=O)',       'D-Gln'),
    'e': ('N[C@@H](CCC(=O)O)C(=O)',       'D-Glu'),
    'g': ('NCC(=O)',                       'Gly'),
    'h': ('N[C@@H](Cc1cnc[nH]1)C(=O)',   'D-His'),
    'i': ('N[C@@H]([C@@H](C)CC)C(=O)',   'D-Ile'),
    'l': ('N[C@@H](CC(C)C)C(=O)',         'D-Leu'),
    'k': ('N[C@@H](CCCCN)C(=O)',          'D-Lys'),
    'm': ('N[C@@H](CCSC)C(=O)',           'D-Met'),
    'f': ('N[C@@H](Cc1ccccc1)C(=O)',      'D-Phe'),
    'p': ('N1CCC[C@@H]1C(=O)',            'D-Pro'),
    's': ('N[C@@H](CO)C(=O)',             'D-Ser'),
    't': ('N[C@@H]([C@@H](O)C)C(=O)',    'D-Thr'),
    'w': ('N[C@@H](Cc1c[nH]c2ccccc12)C(=O)', 'D-Trp'),
    'y': ('N[C@@H](Cc1ccc(O)cc1)C(=O)',  'D-Tyr'),
    'v': ('N[C@@H](C(C)C)C(=O)',          'D-Val'),
}

def build_smiles(sequence, acetylated=False):
    """
    Build D-peptide SMILES.
    sequence: lowercase single-letter (e.g. 'lgrmg')
    acetylated: add N-terminal acetyl group (Ac-)
    Returns canonical SMILES string.
    """
    sequence = sequence.lower().strip()
    parts = []
    for i, aa in enumerate(sequence):
        if aa not in D_AA:
            raise ValueError(f"Unknown amino acid '{aa}'. Use lowercase single-letter D-amino acid codes.")
        smi, _ = D_AA[aa]
        if i == len(sequence) - 1:
            # C-terminus: add -OH
            if smi.endswith('C(=O)'):
                smi = smi[:-5] + 'C(=O)O'
            elif 'N1' in smi:  # Pro
                smi = smi + 'O'
        parts.append(smi)

    full_smiles = ''.join(parts)

    # N-terminal acetylation
    if acetylated:
        full_smiles = 'CC(=O)' + full_smiles[1:]  # replace N with CC(=O)N

    mol = Chem.MolFromSmiles(full_smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES generated for sequence '{sequence}'")
    return Chem.MolToSmiles(mol)

def get_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        'MW':   round(Descriptors.MolWt(mol), 1),
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'TPSA': round(Descriptors.TPSA(mol), 1),
        'HBD':  Descriptors.NumHDonors(mol),
        'HBA':  Descriptors.NumHAcceptors(mol),
        'RotB': Descriptors.NumRotatableBonds(mol),
    }

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build D-peptide SMILES')
    parser.add_argument('--seq', required=True, help='D-peptide sequence (lowercase)')
    parser.add_argument('--acetylated', action='store_true', help='Add N-terminal acetyl group')
    args = parser.parse_args()

    smiles = build_smiles(args.seq, args.acetylated)
    desc = get_descriptors(smiles)
    print(f"Sequence:  {args.seq}")
    print(f"SMILES:    {smiles}")
    print(f"MW:        {desc['MW']} Da")
    print(f"LogP:      {desc['LogP']}")
    print(f"TPSA:      {desc['TPSA']} A^2")