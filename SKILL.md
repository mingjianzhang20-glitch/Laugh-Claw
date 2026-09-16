# Computational D-Peptide Drug Design Skill

## Manuscript reporting status

This file contains workflow examples from an end-to-end process that was previously tested in its original environment and produced repeatable outputs with matched inputs, dependencies, checkpoints, and run settings. The corrected manuscript's result identities, affinity interpretation, and evidence limits are documented in [README.md](README.md). This documentation update did not perform a new model rerun. Model scores and poses are computational hypotheses; neither the nominal affinity threshold nor workflow reproducibility demonstrates experimental potency.

## Description
This skill enables Claude Code to perform end-to-end computational D-peptide drug design including:
- Boltz2-based IC50 affinity prediction for D-peptide candidates
- Multi-objective optimization (MPO) with LogP, TPSA, IC50 criteria
- RDKit physicochemical descriptor calculation
- PyMOL structural visualization of D-peptide/receptor complexes
- Molecular dynamics simulation with OpenMM (AMBER14/CHARMM36m)

## Trigger conditions
Use this skill when asked to:
- Design or screen D-peptide inhibitors against a protein target
- Calculate IC50, LogP, TPSA, MW for peptide candidates
- Run MPO analysis on peptide libraries
- Visualize D-peptide binding poses
- Run MD simulations of peptide-protein complexes

## Requirements
- Python 3.10+
- boltz (for Boltz2 predictions)
- rdkit==2022.09.5
- openmm==8.5.1
- pdbfixer
- Install: `conda env create -f environment.yml`

## Workflow

### Step 1: Prepare input
- Receptor: FASTA sequence or PDB file
- Ligand: peptide sequence (D-amino acids as lowercase)

### Step 2: Generate D-peptide SMILES
```python
from rdkit import Chem
# Historical example: @/@@ symbols alone do not establish L/D identity.
# See examples/smiles_builder.py
```

### Step 3: Boltz2 IC50 prediction
```yaml
# Create YAML (see examples/template_affinity.yaml)
version: 1
sequences:
  - protein:
      id: A
      sequence: RECEPTOR_SEQUENCE
  - ligand:
      id: B
      smiles: "D_PEPTIDE_SMILES"
properties:
  - affinity:
      binder: B
```
```bash
boltz predict input.yaml \
    --out_dir results/ \
    --model boltz2 \
    --affinity_checkpoint ~/.boltz/boltz2_aff.ckpt \
    --sampling_steps 200 \
    --diffusion_samples 1
```

### Step 4: MPO analysis
```python
# Criteria:
# C1: 1 <= LogP <= 3
# C2: TPSA < lead_peptide_TPSA
# C3: IC50 < 1000 nM
# See examples/mpo_analysis.py
```

### Step 5: MD simulation (optional)
```python
# See examples/openmm_simulation.py
# Uses AMBER14 force field for L-peptides
# Uses CHARMM36m for D-peptides via GROMACS topology
```

## File structure
skill/
├── SKILL.md              # This file
├── CLAUDE.md             # Claude Code instructions
├── environment.yml       # Conda environment
├── requirements.txt      # Pip dependencies
├── examples/
│   ├── template_affinity.yaml    # Boltz2 YAML template
│   ├── smiles_builder.py         # Build D-peptide SMILES
│   ├── mpo_analysis.py           # MPO scoring pipeline
│   ├── openmm_simulation.py      # MD simulation
│   └── run_pipeline.py           # End-to-end pipeline
└── references/
├── boltz2-scoring.md
├── d-peptide-design.md
└── troubleshooting.md
