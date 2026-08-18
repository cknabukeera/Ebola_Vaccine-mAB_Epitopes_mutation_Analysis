# Glycosylation

From AF3 predictions, I generated a BDBV GP_REGN3479(Maftivimab) complex prediction

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

# Setting up the scripts
#### add_glycan_to_af3.xml
  - add_glycan_to_af3.xml
```
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="glycan_score" weights="glycan.wts">
            <Reweight scoretype="fa_rep" weight="0.55"/>
            <Reweight scoretype="fa_sol" weight="0.65"/>
            <Reweight scoretype="mm_lj_intra_rep" weight="0.55"/>
            <Reweight scoretype="fa_elec" weight="0.75"/>
        </ScoreFunction>
    </SCOREFXNS>

    <RESIDUE_SELECTORS>
        <!-- Asn563 is in Chain A at residue 563 -->
        <Index name="asn563" resnums="563A"/>
        
        <!-- Select all glycan residues after modeling -->
        <GlycanSelector name="all_glycans" include_protein="false"/>
        
        <!-- GP is Chain A -->
        <Chain name="gp_chains" chains="A"/>
        
        <!-- REGN3479 Fab: Chain B = Heavy, Chain C = Light -->
        <Chain name="fab_chains" chains="B,C"/>
        
        <!-- Interface between GP (A) and Fab (B,C) -->
        <Interface name="gp_fab_interface" 
                   chain1="A" 
                   chain2="B,C"/>
    </RESIDUE_SELECTORS>

    <TASKOPERATIONS>
        <OperateOnCertainResidues name="repack_glycans">
            <ResidueSelector selector="all_glycans"/>
            <ResidueOperation>
                <RestrictToRepacking/>
            </ResidueOperation>
        </OperateOnCertainResidues>
        
        <OperateOnCertainResidues name="repack_nearby">
            <ResidueSelector selector="asn563"/>
            <ResidueOperation>
                <NeighborhoodResidueOperation neighborhood_distance="8.0">
                    <RestrictToRepacking/>
                </NeighborhoodResidueOperation>
            </ResidueOperation>
        </OperateOnCertainResidues>
    </TASKOPERATIONS>

    <MOVERS>
        <!-- Add glycan to Asn563 -->
        <SimpleGlycosylateMover name="add_glycan" 
                               residue_selector="asn563" 
                               glycan="man9"
                               attach_to_ndi="true"/>
        
        <!-- Model glycan conformation -->
        <GlycanModeler name="model_glycan" 
                      layer_size="2" 
                      window_size="1" 
                      rounds="3" 
                      refine="false"
                      scorefxn="glycan_score"/>
        
        <!-- Sample glycan conformations -->
        <GlycanSampler name="sample_glycan"
                      rounds="50"
                      kt="1.0"
                      scorefxn="glycan_score"/>
        
        <!-- Relax the glycan -->
        <GlycanRelax name="relax_glycan"
                    rounds="5"
                    temperature="1.0"
                    scorefxn="glycan_score"/>
        
        <!-- Minimize the full complex -->
        <MinMover name="minimize_complex"
                 scorefxn="glycan_score"
                 type="lbfgs_armijo_nonmonotone"
                 tolerance="0.001"
                 max_iter="500">
            <ResidueSelector>
                <Or selector1="all_glycans" selector2="asn563"/>
            </ResidueSelector>
        </MinMover>
        
        <!-- Analyze interface -->
        <InterfaceAnalyzerMover name="analyze_interface"
                               interface="gp_fab_interface"
                               scorefxn="glycan_score"
                               pack_separated="true"
                               compute_interface_sc="true"/>
    </MOVERS>

    <PROTOCOLS>
        <Add mover="add_glycan"/>
        <Add mover="model_glycan"/>
        <Add mover="sample_glycan"/>
        <Add mover="relax_glycan"/>
        <Add mover="minimize_complex"/>
        <Add mover="analyze_interface"/>
    </PROTOCOLS>

    <OUTPUT scorefxn="glycan_score"/>
</ROSETTASCRIPTS>
```
