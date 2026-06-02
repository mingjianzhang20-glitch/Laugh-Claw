# Computational D-Peptide Drug Design Skill

A Claude Code Skill for end-to-end computational D-peptide inhibitor design.

## Quick Start

```bash
# 1. Clone and setup
git clone https://github.com/mingjianzhang20-glitch/computational-drug-design-skill.git
cd computational-drug-design-skill
conda env create -f environment.yml
conda activate drug_design

# 2. Download Boltz2 checkpoint (~500MB)
boltz setup --model boltz2

# 3. Run full pipeline
python examples/run_pipeline.py \
    --receptor YOUR_RECEPTOR_SEQUENCE \
    --peptides peptide_list.txt \
    --output results/
```

## What it does

| Step | Script | Description |
|------|--------|-------------|
| 1. Build SMILES | `examples/smiles_builder.py` | D-peptide sequence → SMILES |
| 2. Predict IC50 | `examples/boltz2_predict.py` | Boltz2 affinity prediction |
| 3. MPO scoring  | `examples/mpo_analysis.py`   | LogP + TPSA + IC50 filter |
| 4. MD simulation| `examples/openmm_simulation.py` | OpenMM 100-500ns MD |
| 5. Full pipeline| `examples/run_pipeline.py`   | End-to-end automation |

## Individual script usage

```bash
# Build D-peptide SMILES
python examples/smiles_builder.py --seq lgrmg
python examples/smiles_builder.py --seq ffflggqpyw --acetylated

# Predict IC50 with Boltz2
python examples/boltz2_predict.py \
    --receptor APTLFRL \
    --smiles "N[C@@H](CC(C)C)C(=O)..." \
    --name d-lgrmg --output results/

# MPO analysis
python examples/mpo_analysis.py \
    --input candidates.csv \
    --lead_tpsa 263.3 \
    --output mpo_results.csv

# MD simulation
python examples/openmm_simulation.py \
    --pdb complex.pdb \
    --output md_results/ \
    --ns 100
```

## MPO Criteria

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| C1: LogP  | 1 – 3     | Membrane permeability |
| C2: TPSA  | < lead peptide TPSA | More compact than reference |
| C3: IC50  | < 1000 nM | Nanomolar potency |

## IC50 formula (Boltz2)

```python
IC50_nM = 10**(-affinity_pred_value) * 1000
```

## Tested targets

| Target | Lead D-peptide | IC50 |
|--------|---------------|------|
| NDUFA9 | d-LGRMG | 45.2 nM |
| NOTCH1 | d-SSQCF | 574.2 nM |
| 2VSM   | d-GITLGGGS | 44.8 nM |

## Requirements

- Python 3.10+
- Boltz2 checkpoint: `~/.boltz/boltz2_aff.ckpt`
- GPU recommended for Boltz2
- See `environment.yml` for full dependencies
