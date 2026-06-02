"""
OpenMM MD simulation for protein-peptide complex.
Usage: python openmm_simulation.py --pdb complex.pdb --output md_results/
Note: D-peptides are treated as L-peptides for force field purposes (AMBER14).
"""
import argparse, os, sys

def run_simulation(pdb_path, output_dir, production_ns=100, platform='OpenCL'):
    try:
        from openmm.app import PDBFile, ForceField, Modeller, Simulation, DCDReporter, StateDataReporter, PME, HBonds
        from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
        from openmm.unit import kelvin, picosecond, picoseconds, nanometer, bar, molar
        from pdbfixer import PDBFixer
    except ImportError:
        print("ERROR: OpenMM not installed. Run: conda install -c conda-forge openmm pdbfixer")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Fix structure
    print("[1] Fixing structure with PDBFixer...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    print(f"    Atoms: {fixer.topology.getNumAtoms()}")

    # Step 2: Solvate
    print("[2] Adding solvent...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(forcefield, padding=1.0*nanometer, ionicStrength=0.15*molar)
    print(f"    Total atoms: {modeller.topology.getNumAtoms()}")

    # Step 3: Create system
    print("[3] Creating system...")
    system = forcefield.createSystem(modeller.topology,
        nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer,
        constraints=HBonds, rigidWater=True)
    system.addForce(MonteCarloBarostat(1*bar, 303*kelvin))

    # Step 4: Platform
    try:
        plat = Platform.getPlatformByName(platform)
    except:
        plat = Platform.getPlatformByName('CPU')
        print(f"    {platform} unavailable, using CPU")

    # Step 5: Simulation
    integrator = LangevinMiddleIntegrator(303*kelvin, 1/picosecond, 0.002*picoseconds)
    sim = Simulation(modeller.topology, system, integrator, plat)
    sim.context.setPositions(modeller.positions)

    # Step 6: EM
    print("[4] Energy minimization...")
    sim.minimizeEnergy(maxIterations=5000)

    # Step 7: NVT 100ps
    print("[5] NVT equilibration (100 ps)...")
    sim.context.setVelocitiesToTemperature(303*kelvin)
    sim.step(50000)

    # Step 8: Production
    steps = int(production_ns * 500000)  # 0.002ps timestep
    print(f"[6] Production MD ({production_ns} ns = {steps} steps)...")
    sim.reporters.append(DCDReporter(f'{output_dir}/trajectory.dcd', 50000))
    sim.reporters.append(StateDataReporter(f'{output_dir}/production.csv', 5000,
        step=True, time=True, potentialEnergy=True, temperature=True,
        volume=True, density=True, progress=True, totalSteps=steps, speed=True))
    sim.reporters.append(StateDataReporter(sys.stdout, 500000,
        step=True, time=True, temperature=True, speed=True,
        progress=True, totalSteps=steps))
    sim.step(steps)
    print(f"Done! Trajectory: {output_dir}/trajectory.dcd")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='OpenMM MD simulation')
    parser.add_argument('--pdb',        required=True, help='Input PDB file (complex)')
    parser.add_argument('--output',     default='md_results/')
    parser.add_argument('--ns',         type=float, default=100, help='Production MD length in ns')
    parser.add_argument('--platform',   default='OpenCL', choices=['CUDA','OpenCL','CPU'])
    args = parser.parse_args()
    run_simulation(args.pdb, args.output, args.ns, args.platform)