cat > 7TN9_ddg_full_with_errors.py << 'EOF'
#!/usr/bin/env python3
"""
7TN9 Full Complex ddG Calculator with Error Bars
For the complete Inmazeb-GP complex (all mAbs together)
"""
import sys
import argparse
import random
import numpy as np
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.utility import vector1_bool

ALL_AAS = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

# Controls for full 7TN9 complex (chain T and V)
NEUTRAL_CONTROLS = [
    # Chain T controls (GP1 - REGN3470/3471 region)
    (112, "T", "E", "Conservative E→D (neutral)"),
    (116, "T", "P", "Conservative P→A (neutral)"),
    (263, "T", "S", "Conservative S→N (neutral)"),
    # Chain V controls (GP2 - REGN3479 region)
    (504, "V", "I", "Conservative I→V (neutral)"),
]

DESTABILIZING_CONTROLS = [
    # Chain T controls
    (265, "T", "K", "Known destabilizing mutation"),
    (280, "T", "E", "Known destabilizing mutation"),
    # Chain V controls
    (505, "V", "V", "Known destabilizing mutation"),
    (507, "V", "N", "Known destabilizing mutation"),
]

def set_random_seed(seed):
    """Set random seeds for reproducibility"""
    random.seed(seed)
    np.random.seed(seed)
    pyrosetta.rosetta.basic.random.init_random_generators(seed, "mt19937")

def mutate_residue(pose, resnum, mutant_aa, pack_radius, scorefxn):
    """Mutate a residue and repack side-chains within a given radius."""
    if not pose.is_fullatom():
        raise IOError("Full-atom pose required")
    
    new_pose = Pose()
    new_pose.assign(pose)
    task = TaskFactory.create_packer_task(new_pose)
    
    for i in range(1, new_pose.total_residue() + 1):
        if i == resnum:
            allowed = vector1_bool()
            for aa in ALL_AAS:
                allowed.append(False)
            mutant_index = aa_from_oneletter_code(mutant_aa)
            allowed[mutant_index] = True
            task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed)
        else:
            task.nonconst_residue_task(i).restrict_to_repacking()
    
    if pack_radius > 0:
        target_xyz = new_pose.residue(resnum).xyz("CA")
        for i in range(1, new_pose.total_residue() + 1):
            if i != resnum:
                dist = target_xyz.distance(new_pose.residue(i).xyz("CA"))
                if dist > pack_radius:
                    task.nonconst_residue_task(i).prevent_repacking()
    
    pack_mover = PackRotamersMover(scorefxn, task)
    pack_mover.apply(new_pose)
    return new_pose

def calculate_ddg_single(pose, resnum, mutant_aa, pack_radius, scorefxn, relax_rounds, seed=None):
    """Single trajectory ddG calculation"""
    if seed is not None:
        set_random_seed(seed)
    
    # Relax WT
    wt_pose = Pose()
    wt_pose.assign(pose)
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.max_iter(relax_rounds)
    relax.apply(wt_pose)
    native_score = scorefxn(wt_pose)
    
    # Mutate and relax
    mut_pose = Pose()
    mut_pose.assign(pose)
    mut_pose = mutate_residue(mut_pose, resnum, mutant_aa, pack_radius, scorefxn)
    relax.apply(mut_pose)
    mutant_score = scorefxn(mut_pose)
    
    ddG = mutant_score - native_score
    return ddG, native_score, mutant_score

def calculate_ddg_with_errors(pose, resnum, mutant_aa, pack_radius, scorefxn, 
                              relax_rounds, n_trajectories=35, base_seed=42):
    """Run multiple trajectories and return mean ± std dev"""
    ddGs = []
    
    print(f"  Running {n_trajectories} trajectories...")
    for i in range(n_trajectories):
        seed = base_seed + i * 100
        print(f"    Trajectory {i+1}/{n_trajectories} (seed={seed})...", end=" ", flush=True)
        
        ddG, nat_score, mut_score = calculate_ddg_single(
            pose, resnum, mutant_aa, pack_radius, scorefxn, relax_rounds, seed=seed
        )
        ddGs.append(ddG)
        print(f"ddG={ddG:.3f}")
    
    mean_ddG = np.mean(ddGs)
    std_ddG = np.std(ddGs)
    sem_ddG = std_ddG / np.sqrt(n_trajectories)
    
    print(f"    → Mean ddG = {mean_ddG:.3f} ± {std_ddG:.3f} (SD), SEM={sem_ddG:.3f}")
    
    return mean_ddG, std_ddG, sem_ddG, ddGs

def get_internal_resnum(pose, resnum, chain):
    """Convert PDB resnum+chain to internal index"""
    for i in range(1, pose.total_residue() + 1):
        pdb_num = pose.pdb_info().number(i)
        pdb_chain = pose.pdb_info().chain(i)
        if pdb_num == resnum and pdb_chain == chain:
            return i
    return None

def main():
    parser = argparse.ArgumentParser(description="7TN9 Full Complex ddG with error bars.")
    parser.add_argument("-p", "--pdb", required=True, help="Input PDB file (7TN9_clean.pdb)")
    parser.add_argument("-r", "--residue", type=int, required=True, help="PDB residue number")
    parser.add_argument("-c", "--chain", required=True, help="Chain ID (T for GP1, V for GP2)")
    parser.add_argument("-m", "--mutant", required=True, help="Mutant amino acid")
    parser.add_argument("--radius", type=float, default=8.0, help="Repack radius")
    parser.add_argument("--relax_rounds", type=int, default=5, help="FastRelax cycles")
    parser.add_argument("--n_trajectories", type=int, default=35, help="Number of trajectories")
    parser.add_argument("-o", "--output", default="ddg_results.csv", help="Output file")
    parser.add_argument("--run_controls", action="store_true", help="Include control mutations")
    parser.add_argument("--output_controls", default="controls.csv", help="Controls output file")
    
    args = parser.parse_args()
    
    init()
    scorefxn = get_fa_scorefxn()
    pose = pose_from_pdb(args.pdb)
    print(f"Loaded PDB: {args.pdb}, Total residues: {pose.total_residue()}")
    
    internal_resnum = get_internal_resnum(pose, args.residue, args.chain)
    if internal_resnum is None:
        print(f"ERROR: Residue {args.residue} on chain {args.chain} not found.")
        sys.exit(1)
    
    native_aa = pose.residue(internal_resnum).name1()
    print(f"\nTarget: Chain {args.chain} Residue {args.residue} -> {native_aa}")
    print(f"Mutating to: {args.mutant.upper()}")
    print(f"Trajectories: {args.n_trajectories}")
    
    mean_ddG, std_ddG, sem_ddG, all_ddGs = calculate_ddg_with_errors(
        pose, internal_resnum, args.mutant.upper(),
        args.radius, scorefxn, args.relax_rounds, 
        n_trajectories=args.n_trajectories,
        base_seed=42
    )
    
    with open(args.output, "w") as f:
        f.write("# resnum\tchain\twt\tmut\tn_trajectories\tmean_ddG\tstd_ddG\tsem_ddG\tall_ddGs\n")
        f.write(f"{args.residue}\t{args.chain}\t{native_aa}\t{args.mutant.upper()}\t"
                f"{args.n_trajectories}\t{mean_ddG:.3f}\t{std_ddG:.3f}\t{sem_ddG:.3f}\t")
        f.write(",".join([f"{x:.3f}" for x in all_ddGs]) + "\n")
    
    if args.run_controls:
        print("\n" + "="*60)
        print("RUNNING CONTROL MUTATIONS")
        print("="*60)
        
        with open(args.output_controls, "w") as f:
            f.write("# control_mutation\ttype\tdescription\tmean_ddG\tstd_ddG\tsem_ddG\n")
            
            print("\nNeutral controls:")
            for resnum, chain, mut_aa, desc in NEUTRAL_CONTROLS:
                internal_idx = get_internal_resnum(pose, resnum, chain)
                if internal_idx is None:
                    print(f"  ⚠️ Control residue {resnum}{chain} not found, skipping")
                    continue
                mean, std, sem, _ = calculate_ddg_with_errors(
                    pose, internal_idx, mut_aa, args.radius, scorefxn, 
                    args.relax_rounds, n_trajectories=min(3, args.n_trajectories)
                )
                f.write(f"{resnum}{chain}{mut_aa}\tneutral\t{desc}\t{mean:.3f}\t{std:.3f}\t{sem:.3f}\n")
            
            print("\nDestabilizing controls:")
            for resnum, chain, mut_aa, desc in DESTABILIZING_CONTROLS:
                internal_idx = get_internal_resnum(pose, resnum, chain)
                if internal_idx is None:
                    print(f"  ⚠️ Control residue {resnum}{chain} not found, skipping")
                    continue
                mean, std, sem, _ = calculate_ddg_with_errors(
                    pose, internal_idx, mut_aa, args.radius, scorefxn, 
                    args.relax_rounds, n_trajectories=min(3, args.n_trajectories)
                )
                f.write(f"{resnum}{chain}{mut_aa}\tdestabilizing\t{desc}\t{mean:.3f}\t{std:.3f}\t{sem:.3f}\n")
    
    print("\n" + "="*60)
    print(f"RESULT: {native_aa}{args.residue}{args.mutant.upper()} -> ΔΔG = {mean_ddG:.3f} ± {std_ddG:.3f} (SD)")
    print(f"        SEM = {sem_ddG:.3f}, n = {args.n_trajectories}")
    print("="*60)
    print(f"Saved to: {args.output}")
    if args.run_controls:
        print(f"Controls saved to: {args.output_controls}")

if __name__ == "__main__":
    main()
EOF

chmod +x 7TN9_ddg_full_with_errors.py
