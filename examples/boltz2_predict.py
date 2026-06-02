"""
Run Boltz2 affinity prediction for a D-peptide against a receptor.
Usage: python boltz2_predict.py --receptor SEQ --smiles SMILES --name myligand --output results/
"""
import argparse, os, subprocess, json, glob

def create_yaml(name, receptor_seq, ligand_smiles, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    yaml = f"""version: 1
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
    path = os.path.join(output_dir, f"{name}.yaml")
    with open(path, 'w') as f:
        f.write(yaml)
    return path

def run_boltz2(yaml_path, output_dir, checkpoint='~/.boltz/boltz2_aff.ckpt',
               sampling_steps=200, devices=1):
    cmd = [
        'boltz', 'predict', yaml_path,
        '--out_dir', output_dir,
        '--model', 'boltz2',
        '--affinity_checkpoint', os.path.expanduser(checkpoint),
        '--sampling_steps', str(sampling_steps),
        '--diffusion_samples', '1',
        '--devices', str(devices),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Boltz2 error: {result.stderr[-500:]}")
    return result.returncode == 0

def extract_ic50(result_dir):
    """Extract IC50 from Boltz2 affinity JSON output."""
    jsons = glob.glob(os.path.join(result_dir, '**/affinity*.json'), recursive=True)
    if not jsons:
        return None, None
    values = []
    for j in jsons:
        with open(j) as f:
            data = json.load(f)
        aff = data.get('affinity_pred_value')
        if aff is not None:
            # IC50 formula: 10^(-affinity_pred_value) * 1000
            values.append(10**(-aff) * 1000)
    if not values:
        return None, None
    mean_ic50 = round(sum(values)/len(values), 1)
    return mean_ic50, values

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boltz2 affinity prediction')
    parser.add_argument('--receptor', required=True, help='Receptor amino acid sequence')
    parser.add_argument('--smiles',   required=True, help='Ligand SMILES string')
    parser.add_argument('--name',     default='ligand', help='Name for output files')
    parser.add_argument('--output',   default='boltz2_results/', help='Output directory')
    parser.add_argument('--checkpoint', default='~/.boltz/boltz2_aff.ckpt')
    parser.add_argument('--steps',    type=int, default=200)
    parser.add_argument('--devices',  type=int, default=1)
    args = parser.parse_args()

    print(f"Running Boltz2 for {args.name}...")
    yaml_path = create_yaml(args.name, args.receptor, args.smiles, args.output)
    success = run_boltz2(yaml_path, args.output, args.checkpoint, args.steps, args.devices)

    if success:
        ic50, values = extract_ic50(args.output)
        print(f"IC50 = {ic50} nM (from {len(values)} predictions)")
    else:
        print("Boltz2 prediction failed.")