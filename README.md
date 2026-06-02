# Computational D-Peptide Drug Design Skill

A Claude Code skill for end-to-end computational D-peptide inhibitor design.

## Quick Start

```bash
# 1. Clone and setup
git clone https://github.com/mingjianzhang20-glitch/computational-drug-design-skill.git
cd computational-drug-design-skill
conda env create -f environment.yml
conda activate drug_design

# 2. Run pipeline
python examples/run_pipeline.py \
    --receptor YOUR_RECEPTOR_SEQUENCE \
    --peptides peptide_list.txt \
    --output results/
```

## What it does
1. Builds D-peptide SMILES from sequence
2. Predicts IC50 via Boltz2
3. Calculates LogP, TPSA, MW via RDKit
4. Scores candidates with MPO (multi-objective optimization)
5. Ranks and exports results

## MPO Criteria
- C1: 1 ≤ LogP ≤ 3 (membrane permeability)
- C2: TPSA < lead peptide TPSA (compactness)
- C3: IC50 < 1000 nM (potency)

## Tested on
- NDUFA9 target: lead d-LGRMG (IC50=45.2nM)
- NOTCH1 target: lead d-SSQCF (IC50=574.2nM)
- 2VSM target: lead d-GITLGGGS (IC50=44.8nM)

## Requirements
- Boltz2 model checkpoint (~/.boltz/boltz2_aff.ckpt)
- GPU recommended for Boltz2 inference
- See environment.yml for full dependencies
