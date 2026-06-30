#!/usr/bin/env python3

#This scripts  cleans the PDB structure, removes water
"""
clean_pdb.py - Remove water, ligands, and non-protein atoms from any PDB.
Optionally keep only specific chains.

Pure text-filtering approach -> works on ANY PDB, no rebuilding bugs.
"""
import argparse
import sys

# Standard 20 amino acids (3-letter codes)
STANDARD_AAS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

def main():
    parser = argparse.ArgumentParser(description="Clean any PDB for Rosetta")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file")
    parser.add_argument("-o", "--output", required=True, help="Output cleaned PDB")
    parser.add_argument("-c", "--keep_chains", default=None,
                        help="Comma-separated chains to keep (e.g., A,B). If not set, keeps all.")
    parser.add_argument("--remove_water", action="store_true", default=True,
                        help="Remove water (HOH)")
    parser.add_argument("--remove_ligands", action="store_true", default=True,
                        help="Remove non-canonical ligands/ions")
    parser.add_argument("--keep_hetatm", action="store_true", default=False,
                        help="Keep HETATM lines (by default, only ATOM lines are kept)")
    
    args = parser.parse_args()

    # Parse chain filter
    keep_chains = set()
    if args.keep_chains:
        keep_chains = set(args.keep_chains.upper().split(','))
        print(f"Keeping only chains: {', '.join(keep_chains)}")

    filtered_lines = []
    total_lines = 0
    kept_lines = 0

    with open(args.input, 'r') as f:
        for line in f:
            total_lines += 1
            
            # If not keeping HETATM, only accept ATOM lines
            if not args.keep_hetatm and not line.startswith("ATOM"):
                continue
            
            # Keep TER lines (they mark chain ends)
            if line.startswith("TER"):
                filtered_lines.append(line)
                kept_lines += 1
                continue

            # Only process ATOM or HETATM lines
            if not line.startswith(("ATOM", "HETATM")):
                continue

            # Extract chain (PDB column 22, 0-based index 21)
            if len(line) < 22:
                continue  # malformed line, skip
            chain = line[21]
            
            # Extract residue name (columns 18-20, 0-based index 17-20)
            resname = line[17:20].strip().upper()

            # --- Chain filter ---
            if keep_chains and chain not in keep_chains:
                continue

            # --- Water filter ---
            if args.remove_water and resname == "HOH":
                continue

            # --- Ligand filter ---
            if args.remove_ligands and resname not in STANDARD_AAS:
                # Check if it's a modified residue we might want to keep (e.g., MSE -> MET)
                # By default, we skip it to be safe.
                continue

            # If we made it here, keep this atom
            filtered_lines.append(line)
            kept_lines += 1

    # Write the cleaned PDB
    with open(args.output, 'w') as f:
        f.writelines(filtered_lines)

    print(f"Input: {args.input} ({total_lines} lines)")
    print(f"Kept: {kept_lines} lines")
    print(f"Cleaned PDB saved to: {args.output}")

if __name__ == "__main__":
    main()
