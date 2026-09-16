# LAUGH-Claw

An LLM-assisted, auditable workflow for computational D-peptide candidate prioritization. The controller coordinates scientific tools and records workflow state; it is not an independent affinity predictor. Candidate rankings and predicted poses are computational hypotheses, not experimentally validated binding outcomes.

## Documentation status

This documentation update aligns the reported results and interpretation with the corrected manuscript. It does not rerun models or validate the executable examples as a complete reproduction of the revised analysis. Historical example commands, environments, and scripts remain in the repository; their presence does not establish that they reproduce every corrected manuscript result.

The public repository is `mingjianzhang20-glitch/Laugh-Claw`. Corrected sequence-level results, representative pose examples, and historical demonstration inputs are distinct records and must not be interchanged.

## Repository setup and historical examples

The setup and usage snippets below are historical examples, not a verified end-to-end reproduction recipe for the corrected manuscript.

```bash
# 1. Clone and setup
git clone https://github.com/mingjianzhang20-glitch/Laugh-Claw.git
cd Laugh-Claw
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
| 2. Affinity scoring | `examples/boltz2_predict.py` | Computational affinity output; see the interpretation below |
| 3. MPO scoring  | `examples/mpo_analysis.py`   | LogP + TPSA + IC50 filter |
| 4. MD simulation| `examples/openmm_simulation.py` | OpenMM 100-500ns MD |
| 5. Full pipeline| `examples/run_pipeline.py`   | End-to-end automation |

## Individual script usage

The sequences in these historical usage snippets are demonstration inputs, not the corrected rank-1 result table below. Molecular-dynamics examples are not evidence of experimental binding or validation of the reported candidates.

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

## Historical MPO heuristics

These are historical example heuristics, not experimentally established developability criteria or a reconstruction of all filters used in the corrected manuscript. Descriptor values alone do not establish permeability, stability, or binding.

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| C1: LogP  | 1 – 3     | Historical descriptor heuristic; not measured permeability |
| C2: TPSA  | < lead peptide TPSA | Historical descriptor comparison; not a structural compactness measurement |
| C3: Predicted nominal IC50 | < 1000 nM | Nominal model-output prioritization threshold; not demonstrated potency |

## Affinity interpretation in the corrected manuscript

The continuous output is `affinity_pred_value`, interpreted in the manuscript as log10(IC50) with IC50 in micromolar units. The reported conversion is:

```python
predicted_nominal_IC50_nM = 10**(affinity_pred_value) * 1000
```

The separate `affinity_probability_binary` field is binder-versus-decoy probability. It is not the continuous affinity value and is not used for sequence-level affinity ranking. Converted values are predicted nominal IC50 estimates, not experimental measurements. Experimental Kd is a different quantity and is not relabeled as IC50.

The previous derived ΔG metric was removed from the revised analysis: it was a deterministic transformation of the same predicted nominal IC50 values and provided no independent binding-free-energy prediction.

## Corrected manuscript rank-1 records

The following author-verified records replace the previous sequence-to-score assignments. Values are rounded to one decimal place; the `d-` prefix denotes the all-D representation, not experimentally confirmed binding.

| Target | Corrected prioritized candidate | Predicted nominal IC50 (nM) |
|--------|---------------|------|
| NDUFA9 | d-VIIPY | 45.2 |
| NOTCH1 | d-MIFPY | 574.2 |
| 2VSM | d-PLPPIKPR | 44.8 |

The corrected library contains 1,007 raw sequence-level records and 897 QC-valid records: 373 NDUFA9, 175 NOTCH1, and 349 2VSM. The former rank-1 labels must not be attached to these corrected affinity values. Historical pose examples can have different identities and are not the corrected affinity ranking.

## Evidence and reproducibility limits

- The frozen PEPBI retrospective benchmark comprised 15 unique natural L-peptide inputs across seven binding groups. Ten of 12 within-group pairwise orderings agreed with the experimental Kd ordering. Across all 15 inputs, Spearman rho was 0.007 (P = 0.980) and Pearson r was -0.084 (P = 0.765). This provides limited ranking evidence, not absolute calibration or independent validation of D-peptide binding.
- No matched manual or static-workflow comparison establishes a speed, cost, reliability, or autonomy advantage for the LLM controller.
- The historical numerical pose-QC thresholds, full selection universe, and candidate-level pass/fail records could not be reconstructed. Displayed poses are qualitative examples, not a reproducible quantitative pose-QC pass set.
- Parts of the historical bootstrap implementation, failure diagnostics, and complete runtime/cost records remain unresolved. Retained provenance is not equivalent to complete end-to-end reproducibility.
- This documentation-only update does not certify the historical executable examples, publish a complete corrected evidence bundle, or establish a fresh end-to-end reproduction. See [the scoring reference](references/boltz2-scoring.md) for the reporting conventions and the distinction from historical examples.

## Requirements

- Python 3.10+
- Boltz2 checkpoint: `~/.boltz/boltz2_aff.ckpt`
- GPU recommended for Boltz2
- See `environment.yml` for full dependencies
