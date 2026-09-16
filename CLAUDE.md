# Claude Code Instructions for D-Peptide Drug Design

## Documentation scope

The reporting conventions below follow the corrected manuscript. The workflow and executable examples were previously tested end to end in the original environment and produced repeatable outputs with matched inputs, dependencies, checkpoints, and run settings. This documentation update did not perform a new model rerun. Preserve source records and distinguish demonstration inputs from corrected manuscript results. Report candidates as computational hypotheses, not validated inhibitors.

## Setup (run once)
```bash
conda env create -f environment.yml
conda activate drug_design
```

## Key conventions
- D-amino acids: lowercase single-letter code (e.g., d-lgrmg = D-Leu-D-Gly-D-Arg-D-Met-D-Gly)
- Literal `@`/`@@` symbols are not universal labels for L/D identity; parsing or inspecting those symbols alone does not validate the intended stereochemistry. Use independently verified records when reporting the manuscript's representative validation.
- All paths: use relative paths, not hardcoded server paths

## IC50 formula (Boltz2 output)
```python
predicted_nominal_IC50_nM = 10**(affinity_pred_value) * 1000
```

`affinity_pred_value` is the continuous log10 affinity output in the manuscript's micromolar convention. `affinity_probability_binary` is a separate binder-versus-decoy probability, not a continuous ranking substitute. Predicted nominal IC50 is not experimental Kd or measured potency. Derived ΔG was removed from the revised analysis. The corrected result identities are recorded in README.md; historical example names must not be reused as corrected rank-1 labels.

## MPO scoring weights

Historical example heuristic only; these weights are not evidence of experimentally established developability or a reconstruction of every manuscript selection stage.
```python
composite = 0.3 * logp_score + 0.2 * tpsa_score + 0.5 * ic50_score
```

## Reporting unresolved problems

Do not infer a biological or computational failure cause from a missing value alone. Preserve available diagnostics and report unresolved causes as unresolved. Do not infer binding, permeability, or a validated pose from a descriptor, model probability, or visual overlay alone.
