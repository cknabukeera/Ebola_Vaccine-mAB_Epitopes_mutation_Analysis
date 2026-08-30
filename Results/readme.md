---
title: "Ebola Monoclonal Antibody Therapeutics vs. Bundibugyo Virus (BDBV) Glycoprotein (GP)"
---

# Overview

This document contains tables assessing the sequence conservation, binding affinity impact, and interface properties of clinical monoclonal antibodies (mAb114/Ebanga, Inmazeb, and MBP134) when mapped against the Bundibugyo virus (BDBV) GP sequence.

---

# Table 1: Epitopes & Contact Residues
**Glycoprotein (GP) epitopes targeted by monoclonal antibody therapeutics developed against Zaire ebolavirus (EBOV).**

Listed are the GP regions and contact residues recognized by the therapeutic antibodies mAb114 (Ebanga), Inmazeb (Maftivimab, Odesivimab, Atoltivimab), and MBP134 (ADI-15878, ADI-15946). These epitopes were mapped onto the BDBV GP sequence to determine whether key antibody-binding residues are conserved or mutated, enabling assessment of the potential impact of sequence variation on therapeutic binding and efficacy.

| Monoclonal Antibodies | Components | Target Region on GP | Contact Residues on EBOV (Epitopes) |
| :--- | :--- | :--- | :--- |
| **Ebanga** *(PMID: 30686586)* | mAb114 | GP1 | L111, E112, I113, K114, K115, P116, D117, G118, S119 |
| **Inmazeb** *(PMID: 36708708)* | Maftivimab (REGN3479) | Beta 1 - Beta 2 loops of GP1 | P34, L43, V45 |
| | | GP1 IFL | I527, G528, L529, A530, F535, G536 |
| | | GP2 HR1 | Q560, E564, Q567, and Glycan Asn563 |
| | | GP2 N-terminal | I504, V505, A507 |
| | Odesivimab (REGN3471) | GP1 | L111, E112, I113, K114, K115, P116, D117, G118 |
| | | GP1 | G143, T144, P146 |
| | Atoltivimab (REGN3470) | GP1 | I260, Y261, T262, S263, G264, K265, R266, S267, N268, T269, T270 |
| | | GP1 | N278, P279, E280, D282, N292 |
| **MBP134** *(PMID: 30184505, 28525755, 30629917)* | ADI-15878 | GP2 IFL | I527, G528, L529, A530, and Asn563 Glycan |
| | ADI-23774 / ADI-15946 | GP1 | E71, G72, N73, G74, V75 |

*Note: **IFL** = Internal Fusion Loop; **HR1** = Heptad Repeat 1; **GP** = Glycoprotein.*

---

# Table 2: Binding Energy Predictions
**Binding Energy Predictions upon introduction of BDBV epitope mutations into EBOV GP-Inmazeb complex.**

This table contains Rosetta ddG energy calculations for the full trimeric complex with mutations, as well as the independent Inmazeb components (Maftivimab/REGN3479, Atoltivimab/REGN3470, Odesivimab/REGN3471) interacting with the GP protein.

## 1. Inmazeb Trimeric Complex
| Position | Chain | Wild Type (WT) | Mutant (MUT) | ddG (kcal/mol) | Wild Type Score | Mutant Score |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 112 | GP1 (T) | E | D | -42.0800 | 964.3300 | 922.2500 |
| 116 | GP1 (T) | P | A | -128.9800 | 773.9900 | 645.0100 |
| 263 | GP1 (T) | S | N | -76.6400 | 552.6800 | 476.0400 |
| 265 | GP1 (T) | K | R | 25.6300 | 844.2700 | 869.9000 |
| 504 | GP2 (V) | I | T | -248.3424 | 995.2991 | 746.9567 |
| 505 | GP2 (V) | V | L | 83.2278 | 908.7707 | 991.9986 |
| 507 | GP2 (V) | A | T | -10.1850 | 849.9541 | 839.7691 |
| **Net Energy**| | | | **-397.3696** | | |

## 2. GP-REG3479 (Maftivimab) Interaction
| Position | Chain | Wild Type (WT) | Mutant (MUT) | ddG (kcal/mol) | Wild Type Score | Mutant Score |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 504 | GP2 (V) | I | T | 232.9000 | 39.1000 | 272.0000 |
| 505 | GP2 (V) | V | L | 302.1000 | 1.4000 | 303.5000 |
| 507 | GP2 (V) | A | T | -360.5846 | 391.2870 | 30.7024 |
| **Net Energy**| | | | **174.4154** | | |

## 3. GP-REGN3470 (Atoltivimab) Interaction
| Position | Chain | Wild Type (WT) | Mutant (MUT) | ddG (kcal/mol) | Wild Type Score | Mutant Score |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 263 | GP1 (T) | S | N | 247.3100 | 552.6600 | 799.9700 |
| 265 | GP1 (T) | K | R | -100.7900 | 687.6200 | 586.8300 |
| **Net Energy**| | | | **146.5200** | | |

## 4. GP-REGN3471 (Odesivimab) Interaction
| Position | Chain | Wild Type (WT) | Mutant (MUT) | ddG (kcal/mol) | Wild Type Score | Mutant Score |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 112 | GP1 (T) | E | D | -11.1413 | -61.5728 | -72.7142 |
| 116 | GP1 (T) | P | A | -39.6578 | 27.7207 | -11.9371 |
| **Net Energy**| | | | **-50.7991** | | |

---

# Table 3: Interface Properties
**Interface properties of BDBV GP complexed with Maftivimab (REGN3479) alone versus the full Inmazeb cocktail under glycosylated and un-glycosylated states.**

| Complex | State | Binding Energy (dG_separated, REU) | Interface Packing (packstat) | Glycan Effect (ΔΔG_glycan, REU) |
| :--- | :--- | :---: | :---: | :--- |
| **REGN3479_BDBV-GP** | Un-glycosylated | 606.54 | 0.634 | - |
| **REGN3479_BDBV-GP** | Glycosylated | 609.39 | 0.590 (Δ = -0.044) | +2.85 (+1.65 kcal/mol) *(Destabilizing)* |
| **Inmazeb_BDBV-GP** | Un-glycosylated | 946.52 | 0.721 | - |
| **Inmazeb_BDBV-GP** | Glycosylated | 946.06 (548.71 kcal/mol) | 0.708 (Δ = -0.013) | -0.46 (-0.27 kcal/mol) *(Neutral)* |

