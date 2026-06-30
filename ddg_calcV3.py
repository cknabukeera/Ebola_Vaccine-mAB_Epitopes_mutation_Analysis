#!/usr/bin/env python3
"""
Simple ddG calculator – exactly like your working notebook.
No constraints, no multiple trajectories.
Fixes restrict_absent_canonical_aas() by passing a vector1_bool.
"""
import sys
import argparse
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.utility import vector1_bool

ALL_AAS = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def mutate_residue(pose, resnum, mutant_aa, pack_radius, scorefxn):
    """Mutate a residue and repack side-chains within a given radius."""
    if not pose.is_fullatom():
        raise IOError("Full-atom pose required")
    
    new_pose = Pose()
    new_pose.assign(pose)
    task = TaskFactory.create_packer_task(new_pose)
    
    for i in range(1, new_pose.total_residue() + 1):
        if i == resnum:
            # Build allowed list: only the mutant amino acid is allowed
            allowed = vector1_bool()
            for aa in ALL_AAS:
                allowed.append(False)
            mutant_index = aa_from_oneletter_code(mutant_aa)  # 1-based index
            allowed[mutant_index] = True
            # Restrict to only this amino acid
            task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed)
        else:
            # All other residues can only repack (not mutate)
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

def calculate_ddg(pose, resnum, mutant_aa, pack_radius, scorefxn, relax_rounds, suffix=""):
    """Compute ddG exactly like your notebook: relax WT once, mutate, relax mutant once."""
    print("  Relaxing Wild Type...")
    wt_pose = Pose()
    wt_pose.assign(pose)
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.max_iter(relax_rounds)
    relax.apply(wt_pose)
    native_score = scorefxn(wt_pose)
    wt_pose.dump_pdb(f"wt_relaxed{suffix}.pdb")
    print(f"    Native score: {native_score:.4f} (saved as wt_relaxed{suffix}.pdb)")
    
    print("  Creating Mutant...")
    mut_pose = Pose()
    mut_pose.assign(pose)
    mut_pose = mutate_residue(mut_pose, resnum, mutant_aa, pack_radius, scorefxn)
    
    print("  Relaxing Mutant...")
    relax.apply(mut_pose)
    mutant_score = scorefxn(mut_pose)
    mut_pose.dump_pdb(f"mutant_{resnum}{mutant_aa}{suffix}.pdb")
    print(f"    Mutant score: {mutant_score:.4f} (saved as mutant_{resnum}{mutant_aa}{suffix}.pdb)")
    
    ddG = mutant_score - native_score
    return ddG, native_score, mutant_score

def main():
    parser = argparse.ArgumentParser(description="Simple PyRosetta ddG – notebook version.")
    parser.add_argument("-p", "--pdb", required=True, help="Input PDB file")
    parser.add_argument("-r", "--residue", type=int, required=True, help="PDB residue number")
    parser.add_argument("-c", "--chain", default=None, help="Chain ID (e.g., A).")
    parser.add_argument("-m", "--mutant", required=True, help="Mutant amino acid (one-letter code)")
    parser.add_argument("--radius", type=float, default=8.0, help="Repack radius (Å)")
    parser.add_argument("--relax_rounds", type=int, default=5, help="FastRelax cycles")
    parser.add_argument("-o", "--output", default="ddg_results.txt", help="Output file")
    
    args = parser.parse_args()
    
    init()                         # matches your notebook
    scorefxn = get_fa_scorefxn()   # matches your notebook
    pose = pose_from_pdb(args.pdb)
    print(f"Loaded PDB: {args.pdb}, Total residues: {pose.total_residue()}")
    
    # Map chain+residue to internal index
    internal_resnum = None
    if args.chain:
        for i in range(1, pose.total_residue() + 1):
            pdb_num = pose.pdb_info().number(i)
            pdb_chain = pose.pdb_info().chain(i)
            if pdb_num == args.residue and pdb_chain == args.chain:
                internal_resnum = i
                break
        if internal_resnum is None:
            print(f"ERROR: Residue {args.residue} on chain {args.chain} not found.")
            sys.exit(1)
    else:
        if args.residue > pose.total_residue():
            print(f"ERROR: Residue {args.residue} not found.")
            sys.exit(1)
        internal_resnum = args.residue
    
    native_aa = pose.residue(internal_resnum).name1()
    print(f"Target: Chain {args.chain or 'N/A'} Residue {args.residue} -> {native_aa}")
    print(f"Mutating to: {args.mutant.upper()}")
    print(f"Relax rounds: {args.relax_rounds}, Repack radius: {args.radius} Å")
    
    suffix = f"_{args.residue}_{args.chain}_{args.mutant.upper()}"
    
    ddG, nat_score, mut_score = calculate_ddg(
        pose, internal_resnum, args.mutant.upper(),
        args.radius, scorefxn, args.relax_rounds, suffix
    )
    
    with open(args.output, "w") as f:
        f.write("# resnum\tchain\twt\tmut\tddG\tnative_score\tmutant_score\n")
        f.write(f"{args.residue}\t{args.chain or 'X'}\t{native_aa}\t{args.mutant.upper()}\t{ddG:.4f}\t{nat_score:.4f}\t{mut_score:.4f}\n")
    
    print("")
    print("=" * 60)
    print(f"RESULT: {native_aa}{args.residue}{args.mutant.upper()} -> ddG = {ddG:.4f}")
    print("=" * 60)
    print(f"Saved to: {args.output}")

if __name__ == "__main__":
    main()
