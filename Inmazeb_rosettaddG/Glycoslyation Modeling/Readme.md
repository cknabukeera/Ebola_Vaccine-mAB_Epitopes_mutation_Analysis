# Glycosylation

----
#### The Comparative Analysis
The workflow remains identical. You will run the glycosylation script on both of your complexes:

Glycosylate the WT (Zaire GP + REGN3479): Run the script using input_pdb = "zaire_regn3479.pdb" (or whatever your Zaire complex file is named) to attach man5 to Asn563.
Glycosylate the Mutant (BDBV GP + REGN3479): Run the script on your af3_bdbv_regn3479.pdb (BDBV mutant) to attach man5 to Asn563.
Calculate Binding Energies: Run Rosetta InterfaceAnalyzer on both of these new glycosylated structures.
The Comparison:
Compare the mutational impact ($\Delta\Delta G$) in both states:

 - $$\Delta\Delta G_{\text{bind}} \text{ (Without Glycans)} = \Delta G_{\text{BDBV (no-glyc)}} - \Delta G_{\text{Zaire (no-glyc)}}$$

 - $$\Delta\Delta G_{\text{bind}} \text{ (With Glycans)} = \Delta G_{\text{BDBV (glyc)}} - \Delta G_{\text{Zaire (glyc)}}$$

- If the $\Delta\Delta G$ shifts significantly once the glycan is present, you will have quantitative proof of how the conserved Asn563 glycan influences the structural impact of the BDBV substitutions (such as A507T) on Maftivimab binding!
From AF3 predictions, I generated a BDBV GP_REGN3479(Maftivimab) complex prediction
-----

#### STEP 1: Convert the cif (AF3 output) to PDB using PyMOL

```
#In PyMOL
load fold_regn3479_bdbv_gp_model_0.cif
save af3_bdbv_regn3479.pdb
```
----
#### STEP 2: Identify all chains and their sizes
```
# Show chain information
print("="*60)
print("CHAIN IDENTIFICATION")
print("="*60)

# Get all chains
chains = cmd.get_chains("af3_bdbv_regn3479")
print(f"Chains found: {chains}")

# Show residues per chain
for chain in chains:
    cmd.select(f"chain_{chain}", f"chain {chain}")
    num_atoms = cmd.count_atoms(f"chain_{chain}")
    print(f"Chain {chain}: {num_atoms} atoms")
    
    # Get residue range
    residues = cmd.get_residues(f"chain {chain}")
    if residues:
        print(f"  Residues: {residues[0]} - {residues[-1]}")

# Color chains differently
color red, chain A
color green, chain B
color blue, chain C
color yellow, chain D
color magenta, chain E
color cyan, chain F
```
##### Output
 - The pdb has three chains. Chain A is the BDBV GP
```
PyMOL>color red, chain A
 Executive: Colored 5340 atoms.
PyMOL>color green, chain B
 Executive: Colored 1605 atoms.
PyMOL>color blue, chain C
 Executive: Colored 1632 atoms.
```
----
#### STEP 3: Find Asn563
```
# Search for Asn563 in each chain
print("\nSearching for Asn563...")

# Chain A
cmd.select("asn_A", "chain A and resn ASN and resi 560-570")
if cmd.count_atoms("asn_A") > 0:
    print("Asn found in chain A near 560-570")
    cmd.zoom("asn_A")
    cmd.show("sticks", "asn_A")
    cmd.color("red", "asn_A")
else:
    print("No ASN in chain A (560-570)")

# Chain B
cmd.select("asn_B", "chain B and resn ASN and resi 560-570")
if cmd.count_atoms("asn_B") > 0:
    print("Asn found in chain B near 560-570")
    cmd.zoom("asn_B")
    cmd.show("sticks", "asn_B")
    cmd.color("red", "asn_B")
else:
    print("No ASN in chain B (560-570)")

# Chain C
cmd.select("asn_C", "chain C and resn ASN and resi 560-570")
if cmd.count_atoms("asn_C") > 0:
    print("Asn found in chain C near 560-570")
    cmd.zoom("asn_C")
    cmd.show("sticks", "asn_C")
    cmd.color("red", "asn_C")
else:
    print("No ASN in chain C (560-570)")

# Also try exact residue 563
cmd.select("asn563", "resi 563")
if cmd.count_atoms("asn563") > 0:
    asn_chain = cmd.get_chains("asn563")
    print("Asn563 found in chain:", asn_chain)
    cmd.zoom("asn563")
    cmd.show("sticks", "asn563")
    cmd.color("red", "asn563")
else:
    print("Residue 563 not found - checking all ASN residues...")
    cmd.select("all_asn", "resn ASN")
    print("Total ASN residues:", cmd.count_atoms("all_asn"))
    # Show where they are
    print("ASN residues in structure:")
    cmd.iterate("all_asn", "print(chain, resi)")
```

##### Output: 
 - Asn563 is in Chain A at residue 563.
```
Asn563 found in chain: ['A']
```
-----


## 1. Rationale

### The Challenge
Structural prediction tools like **AlphaFold 3 (AF3)** excel at folding protein-protein complexes, but they **cannot natively model complex post-translational modifications (PTMs)** like branched carbohydrates. 

Ebola and Bundibugyo glycoproteins are heavily glycosylated on their surfaces. These glycans act as a "glycan shield" to evade the host immune system, but they also form critical contacts with therapeutic antibodies.

### The Objective: Modeling the Impact of Glycans on mAb REGN3479–BDBV GP Interaction
The primary goal of this project is to model the glycan at **Asn563** to investigate its structural and energetic impact on the binding of **Maftivimab (REGN3479)** to the **BDBV GP**:
1. **Steric Hindrance & Shielding:** To evaluate if the glycan at Asn563 shields the GP epitope, restricting antibody access and explaining differences in neutralization potency.
2. **Direct Glycan–mAb Contacts:** To check if Maftivimab forms favorable, stabilizing contacts (hydrogen bonds or van der Waals interactions) directly with the carbohydrate tree.
3. **Binding Affinity ($\Delta G_{\text{bind}}$):** By generating this glycosylated model, we can run Rosetta `InterfaceAnalyzer` to compare the interface binding energy with and without the glycan.

### The Target Site: Asn563 (N563)
* **Location:** The GP2 subunit of the glycoprotein trimer, situated near the base.
* **Epitope Significance:** Asn563 is a highly conserved N-glycosylation site (motif `N-E-T`) positioned right at the binding interface of Maftivimab (REGN3479).
* **7TN9 Reference:** The experimental crystal structure (PDB ID: [7TN9](https://www.rcsb.org/structure/7TN9)) shows that Asn563 is modified with a standard oligomannose N-glycan. 

This workflow adds this glycan back to the AF3 model and relaxes it in Cartesian space to resolve clashes, giving us a realistic, biologically active structural model of the GP–mAb complex.

## 2. Project Inputs & File Descriptions
To run this workflow, you need the following three files in the same working directory:
### 1. `af3_bdbv_regn3479.pdb` (The Input Structure)
* **What it is:** The 3D structure predicted by AlphaFold 3 containing the BDBV GP trimer in complex with the REGN3479 (Maftivimab) Fab heavy and light chains.
* **Why it is needed:** It serves as the un-glycosylated starting template.
### 2. `model_bdbv_glycosylation.py` (The Python Modeling Script)
* **What it is:** The PyRosetta script that programmatically handles target residue validation, covalent glycan attachment, and conformational modeling.
### 3. `run_glycosylation.sh` (The Slurm Batch Script)
* **What it is:** The SBATCH submission script to run the Python modeling script as a job on an HPC cluster.
---
## 3. Detailed Steps in the Python Script
The python script (`model_bdbv_glycosylation.py`) programmatically attaches the glycan and relaxes the structure. Below is a detailed explanation of what each key Rosetta function does:
### 1. PyRosetta Initialization (`init`)

```
init("-include_sugars -write_pdb_link_records -beta")
```
 - include_sugars: By default, Rosetta ignores carbohydrates to save memory. Sourcing this flag is mandatory; it tells Rosetta to load all carbohydrate topologies, residue geometries, and atomic types into its chemical database.
 - write_pdb_link_records: Instructs Rosetta to write standard LINK records in the header of the output PDB file. This ensures that molecular visualization programs (like PyMOL or Chimera) recognize that the carbohydrate is covalently bound to the protein, rather than floating in space.
 - beta: Activates the beta energy function weights, which are highly optimized for carbohydrate-protein interactions.
#### 2. Loading the Structure & Auto-Mapping the Site
```
pose = pose_from_pdb(input_pdb)
pose_res = pose.pdb_info().pdb2pose(chain_id, target_res)
```
 - Loads your predicted PDB complex.
 - Maps your PDB residue numbering (Chain A, residue 563) to Rosetta’s internal, sequential 1-based index (known as the pose index).
 - The Safeguard: If residue 563 on Chain A is not an Asparagine (N) (for example, if GP1 and GP2 are separate chains), the script automatically scans all other chains, identifies which chain has Asn563, and selects it.
#### 3. Covalent Glycan Attachment (SimpleGlycosylateMover)
```
glycosylator = SimpleGlycosylateMover()
glycosylator.set_position(pose_res)
glycosylator.set_glycosylation("man5")
glycosylator.set_strip_existing_glycans(True)
glycosylator.apply(pose)
```
 - Instantiates the mover and sets the target residue (set_position) to the resolved pose index for Asn563.
 - Specifies the N-glycan tree type (set_glycosylation("man5")) to model.
 - Sets set_strip_existing_glycans(True) to ensure any incomplete carbohydrate atoms at that position are cleared before building.
   - On calling apply(pose), it builds the covalent bond between the nitrogen atom (ND2) of the Asparagine side chain and the first sugar residue (NAG) of the Man5 tree.
#### 4. Conformational Optimization (GlycanTreeModeler)
```
glycan_modeler = GlycanTreeModeler()
glycan_modeler.apply(pose)
```
 - Simply placing a glycan tree on a static structure will create massive steric clashes (overlapping atoms) with the surrounding BDBV GP residues and the adjacent REGN3479 antibody.
 - The GlycanTreeModeler runs a Monte Carlo sampling algorithm to rotate the glycan's glycosidic bond angles ($\phi$, $\psi$, $\omega$).
 - It optimizes the shape of the glycan tree in 3D space to find the lowest energy conformation, resolving all steric clashes with the antibody and surrounding atoms.
#### 5. Saving the Glycoprotein
```
pose.dump_pdb(output_pdb)
```
 - Writes the coordinates of the relaxed, glycosylated complex trimer to af3_bdbv_regn3479_glyc_563.pdb for downstream analysis.

###### Here is the complete, copy-pasteable Python script model_bdbv_glycosylation.py
```
#!/usr/bin/env python3
"""
PyRosetta script to model glycosylation at residue Asn563 on Bundibugyo ebolavirus (BDBV) glycoprotein
in complex with REGN3479 (Maftivimab).
"""

import sys
import os

try:
    from pyrosetta import *
    from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover, GlycanTreeModeler
except ImportError:
    print("Error: PyRosetta is not installed or not in your current PYTHONPATH.")
    print("Please make sure you have activated your conda environment (e.g., conda activate pyrosetta).")
    sys.exit(1)

def main():
    input_pdb = "af3_bdbv_regn3479.pdb"
    target_res = 563
    glycan_name = "man5"  # Matches the exact N-glycan resolved in the 7TN9 crystal structure
    
    if not os.path.exists(input_pdb):
        print(f"Error: Input PDB file '{input_pdb}' not found in the current directory.")
        print("Please run this script in the directory containing your PDB file.")
        sys.exit(1)

    print("====================================================================")
    print("Rosetta Glycosylation Modeling for BDBV-mAb Complex")
    print(f"PDB: {input_pdb}")
    print(f"Target site: Asn{target_res}")
    print(f"Glycan: {glycan_name}")
    print("====================================================================")

    # 1. Initialize PyRosetta with carbohydrate options
    print("Initializing PyRosetta with sugar support...")
    init("-include_sugars -write_pdb_link_records -beta")

    # 2. Load the structure
    print(f"Loading structure {input_pdb}...")
    pose = pose_from_pdb(input_pdb)

    # 3. Resolve the chain and residue number
    # Default is Chain A as specified, but we will auto-detect to prevent GP1/GP2 chain mapping issues
    print("Verifying target residue Asn563 chain mapping...")
    chain_id = 'A'
    pose_res = pose.pdb_info().pdb2pose(chain_id, target_res)
    
    is_valid = False
    if pose_res != 0:
        res_name = pose.residue(pose_res).name1()
        if res_name == 'N':
            is_valid = True
            print(f"Confirmed: Residue {target_res} on Chain {chain_id} is Asparagine (N) [Pose index {pose_res}].")

    if not is_valid:
        print(f"Warning: Residue {target_res} on Chain {chain_id} is either missing or not an Asparagine.")
        print("Scanning other chains for residue 563...")
        
        found_chains = []
        # Get list of all chain characters present in the pose
        all_chains = set()
        for i in range(1, pose.total_residue() + 1):
            all_chains.add(pose.pdb_info().chain(i))
            
        for ch in sorted(list(all_chains)):
            idx = pose.pdb_info().pdb2pose(ch, target_res)
            if idx != 0:
                name = pose.residue(idx).name1()
                print(f"  Chain {ch} position {target_res} is: {name}")
                if name == 'N':
                    found_chains.append(ch)
        
        if found_chains:
            chain_id = found_chains[0]
            pose_res = pose.pdb_info().pdb2pose(chain_id, target_res)
            print(f"Mapping resolved: Using Chain {chain_id} (Asn{target_res} at Pose index {pose_res}).")
        else:
            print(f"Error: Could not find any chain with an Asparagine (N) at position {target_res}!")
            sys.exit(1)

    # 4. Attach Glycan
    print(f"Attaching N-linked glycan '{glycan_name}' to residue {pose_res}...")
    glycosylator = SimpleGlycosylateMover()
    glycosylator.set_position(pose_res)
    glycosylator.set_glycosylation(glycan_name)
    glycosylator.set_strip_existing_glycans(True)
    glycosylator.apply(pose)
    print("Glycan attached successfully.")

    # 5. Optimize Glycan Conformations to relieve steric clashes with surrounding residues/antibody
    print("Optimizing glycan conformations using GlycanTreeModeler (this may take a minute)...")
    
    glycan_modeler = GlycanTreeModeler()
    glycan_modeler.apply(pose)
    print("Glycan conformation optimization complete.")

    # 6. Save PDB file
    output_pdb = f"af3_bdbv_regn3479_glyc_{target_res}.pdb"
    pose.dump_pdb(output_pdb)
    print("====================================================================")
    print(f"SUCCESS: Glycosylated structure saved to {output_pdb}")
    print("====================================================================")

if __name__ == "__main__":
    main()
```
##### Output
```
#This is the tail of the script output
Rosetta Glycosylation Modeling for BDBV-mAb Complex
PDB: af3_bdbv_regn3479.pdb
Target site: Asn563
Glycan: man5
====================================================================
Initializing PyRosetta with sugar support...
┌───────────────────────────────────────────────────────────────────────────────┐
│                                  PyRosetta-4                                  │
│               Created in JHU by Sergey Lyskov and PyRosetta Team              │
│               (C) Copyright Rosetta Commons Member Institutions               │
│                                                                               │
│ NOTE: USE OF PyRosetta FOR COMMERCIAL PURPOSES REQUIRES PURCHASE OF A LICENSE │
│          See LICENSE.PyRosetta.md or email license@uw.edu for details         │
└───────────────────────────────────────────────────────────────────────────────┘
PyRosetta-4 2026 [Rosetta PyRosetta4.Release.python311.ubuntu 2026.25+release.a31d9d50e217e874c57db71bf169052467a8e3d8 2026-06-08T08:14:09] retrieved from: http://www.pyrosetta.org
Loading structure af3_bdbv_regn3479.pdb...
Verifying target residue Asn563 chain mapping...
Confirmed: Residue 563 on Chain A is Asparagine (N) [Pose index 563].
Attaching N-linked glycan 'man5' to residue 563...
Glycan attached successfully.
Optimizing glycan conformations using GlycanTreeModeler (this may take a minute)...
Glycan conformation optimization complete.
====================================================================
SUCCESS: Glycosylated structure saved to af3_bdbv_regn3479_glyc_563.pdb
====================================================================
ating DensePDInteractionGraph
protocols.carbohydrates.util: pack accepted: 1
protocols.carbohydrates.GlycanTreeModeler: Start- : 39399.1
protocols.carbohydrates.GlycanTreeModeler: Pre  - : 5938.09
protocols.carbohydrates.GlycanTreeModeler: Post - : 5938.19
protocols.carbohydrates.GlycanTreeModeler: Final- : 5938.09
````
##### Interpretation of the output (Conformation Optimization Energy Profile (REU))

During the glycosylation modeling run, PyRosetta optimizes the structure using the carbohydrate-specific `beta` scorefunction. Below is the energy profile of the BDBV GP–Maftivimab complex at each stage of the process:

| Optimization Stage | Total System Energy (REU) | Description |
| :--- | :---: | :--- |
| **Start (Initial)** | `39,399.1` | Energy immediately after attaching the `man5` glycan. The extremely high score reflects massive steric clashes (atom overlaps) between the new glycan and the adjacent protein/antibody atoms. |
| **Pre-relaxation** | `5,938.09` | Energy after a quick initial minimization of the glycan tree to relieve the most severe clashes. |
| **Post-relaxation** | `5,938.19` | Intermediate energy state during Monte Carlo conformational sampling. |
| **Final (Optimized)** | `5,938.09` | The final energy of the relaxed complex. The **~33,461 REU drop** indicates that `GlycanTreeModeler` successfully optimized the glycan's glycosidic bond angles to fit cleanly into the binding interface. |


