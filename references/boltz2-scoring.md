# Boltz2 affinity reporting and historical scoring examples

The corrected manuscript reports computational affinity outputs as predicted nominal IC50 estimates. They are not measured affinities or experimentally validated potency.

The YAML and execution snippets below are retained historical examples. Their checkpoint flags, environment assumptions, and input identities have not been revalidated by this documentation-only update. They must not be represented as a complete, independently reproduced configuration for the revised manuscript. [README.md](../README.md) records the corrected results and remaining reproducibility limits.

---

## Input YAML Format

Each candidate requires one YAML file. The format uses `version: 1` with a protein chain (id: A) and a peptide ligand as SMILES (id: B).

### NDUFA9 example

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: DPVITNPDAKLSVQDTQVK
  - ligand:
      id: B
      smiles: CC(C)[C@@H](NC(=O)[C@H](NC(=O)[C@H](N)CCCCN)[C@H](C)O)C(=O)N[C@H](CC(N)=O)C(=O)N[C@H](Cc1ccc(O)cc1)C(=O)O
properties:
  - affinity:
      binder: B
```

### NOTCH1 example

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: CSPGAPQGPPFPGAPGPKQDYRCGGSFCMGPGCGKVGTCVDIWNSGPCLAVISGAAPDGFLSELGAGVPCRAVCSSVGCCRWGLPCSSTVCPCLVGAGFGFHVQWLCRTCVNTRDLRQWEERLGRPCPRVRDVDLLDRPWGPCLRWGPAPWRLSRDGPCGWDHWSCCASQLVGCCPVAGWRCVDIWDPCLRAGVCCRLVGANGFDRCWTCGAPGWGTCVWDWADCPRAADGPCRAVDWADPCLARGVCSPDGFCADLVCQDLWGPQGANGFVGCCDPGGTCVAGWDRDLSGPCADPWTCADGTCGTCVWDGSDWTPCHV
  - ligand:
      id: B
      smiles: CC[C@@H](C)[C@@H](NC(=O)[C@H](N)CO)C(=O)N[C@@H](C(=O)N[C@H](Cc1ccc(O)cc1)C(=O)N[C@H](CC(C)C)C(=O)O)C(C)C
properties:
  - affinity:
      binder: B
```

---

## Running Predictions

### Environment setup

```bash
source /home/ubuntu/miniconda/bin/activate boltz
```

### Single candidate

```bash
boltz predict candidate_00000.yaml \
    --checkpoint ~/.boltz/boltz2_aff.ckpt \
    --out_dir ./affinity_output \
    --devices 2 \
    --accelerator gpu \
    --diffusion_samples 1 \
    --sampling_steps 200
```

### Batch scoring

```bash
source /home/ubuntu/miniconda/bin/activate boltz

for yaml in yamls/*.yaml; do
    boltz predict "$yaml" \
        --checkpoint ~/.boltz/boltz2_aff.ckpt \
        --out_dir ./affinity_output \
        --devices 2 \
        --accelerator gpu \
        --diffusion_samples 1 \
        --sampling_steps 200 \
        --override
done
```

### Key flags

| Flag | Default | Description |
|------|---------|-------------|
| `--checkpoint` | — | Always use `~/.boltz/boltz2_aff.ckpt` |
| `--out_dir` | `./predictions` | Output directory |
| `--devices` | 1 | Number of GPUs; use 2 for dual RTX 4090 |
| `--accelerator` | gpu | Always use gpu |
| `--diffusion_samples` | 1 | Increase for diversity |
| `--sampling_steps` | 200 | Higher = more accurate, slower |
| `--override` | false | Overwrite existing results |

---

## Parsing Results

The previous snippet incorrectly assigned the raw `affinity_pred_value` directly to a field named `ic50_nM`. That snippet has been withdrawn.

In the corrected manuscript, `affinity_pred_value` is the continuous log10(IC50) output with IC50 expressed in μM. Predicted nominal IC50 in nM is reported as `10^(affinity_pred_value) × 1000`. The raw value is not itself an nM concentration. The separate `affinity_probability_binary` field represents binder-versus-decoy probability and is not used as the continuous affinity ranking value.

Missing or invalid output must not be reported as a successful affinity estimate. Experimental Kd remains a separate experimental quantity. No derived ΔG is used in the revised analysis.

---

## Corrected manuscript records

These author-verified sequence-to-score assignments are computational prioritization records, not experimentally validated binding results. Historical example inputs above are not evidence for these corrected identities.

| Target | Corrected prioritized candidate | Predicted nominal IC50 (nM) |
|--------|---------------|------|
| NDUFA9 | d-VIIPY | 45.2 |
| NOTCH1 | d-MIFPY | 574.2 |
| 2VSM | d-PLPPIKPR | 44.8 |

The corrected manuscript also reports a small natural L-peptide PEPBI benchmark. It supports only a limited statement about within-group ordering, not quantitative calibration of the nominal scale or experimental validation of these D-peptide candidates. See [the evidence summary](../README.md#evidence-and-reproducibility-limits).

---

## Historical troubleshooting notes

The statements below are retained historical notes, not verified failure diagnoses or a validated configuration for the revised manuscript. The manuscript explicitly reports gaps in historical failure diagnostics. Do not infer a root cause or silently substitute settings from these notes.

**Wrong checkpoint**: Always use `~/.boltz/boltz2_aff.ckpt`. The default checkpoint produces wrong results silently.

**CUDA out of memory**: Reduce `--diffusion_samples` to 1 and `--devices` to 1.

**Invalid SMILES**: Validate before submitting:
```python
from rdkit import Chem
assert Chem.MolFromSmiles(smiles) is not None, f"Invalid SMILES: {smiles}"
```
