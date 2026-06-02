# Claude Code Instructions for D-Peptide Drug Design

## Setup (run once)
```bash
conda env create -f environment.yml
conda activate drug_design
```

## Key conventions
- D-amino acids: lowercase single-letter code (e.g., d-lgrmg = D-Leu-D-Gly-D-Arg-D-Met-D-Gly)
- D-peptide SMILES: use @@ stereochemistry at alpha-carbon
- L-peptide SMILES: use @ stereochemistry at alpha-carbon
- N-acetylated peptides: add CC(=O)N- prefix to SMILES
- All paths: use relative paths, not hardcoded server paths

## IC50 formula (Boltz2 output)
```python
IC50_nM = 10**(-affinity_pred_value) * 1000
```

## MPO scoring weights
```python
composite = 0.3 * logp_score + 0.2 * tpsa_score + 0.5 * ic50_score
```

## Common errors and fixes
- SMILES parsing fails: check stereochemistry @@ vs @
- Boltz2 NaN: receptor sequence too short, use full domain
- High TPSA: avoid Q,Y,W,K,R,H residues; prefer F,V,L,I,A
- OpenMM NaN: run restrained EM first before NPT
