# Examples

Example input files for each stage of the D-peptide design pipeline.

---

## Boltz2 Affinity Scoring YAMLs

| File | Target | Protein length |
|------|--------|---------------|
| `ndufa9_boltz2_affinity.yaml` | NDUFA9 | 19 aa |
| `notch1_boltz2_affinity.yaml` | NOTCH1 | ~380 aa |

Run with:
```bash
conda activate drug_design

boltz predict examples/ndufa9_boltz2_affinity.yaml \
    --checkpoint ~/.boltz/boltz2_aff.ckpt \
    --out_dir ./affinity_output \
    --devices 2 \
    --accelerator gpu
```

---

## BoltzGen Design YAMLs

| File | Target | Design length |
|------|--------|--------------|
| `ndufa9_boltzgen_design.yaml` | NDUFA9 | 5-10 residues |

Run with:
```bash
boltzgen run examples/ndufa9_boltzgen_design.yaml \
    --output workbench/ndufa9 \
    --protocol peptide-anything \
    --num_designs 1000 \
    --budget 50
```

> Note: `.cif` file paths in BoltzGen YAMLs are relative to the YAML file location.
