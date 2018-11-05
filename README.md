# Utility Scripts

Most of these scripts are either quality-of-life scripts meant to automate some of the more tedious tasks, or they are data-processing and data-visualization scripts.

#### Hydrogen Bond Visualization
*By: Mark A. Hix*

The program **cpptraj** can be used to generate a list of hydrogen bonding interactions between atoms in a molecular dynamics system.  The generated list is sorted by the number of frames in which the interaction is present, and contains information on the hydrogen bond donor and acceptor atoms and residues.

In protein simulations, it is often necessary to examine these interactions carefully.  However, sifting through large tables can often be tedious, and can occasionally lead to missed points of interest.  The script provided under *HBond_Analysis.py* takes these data files and creates a matrix visualization showing the interactions between residue pairs along a logarithmic scale.  Generated outputs will show maximum values of 1 to indicate the pair with the most interaction through the simulation, and 0 for those with none, scaling logarithmically between them.

#### Correlation Decomposition Analysis
*By: Mark A. Hix*

**Cpptraj** can be used to produce an analysis of correlation between residue pairs in a protein or other biological system.  This correlation determines whether residues are correlated in their movement (positive value), anti-correlated (negative value), or non-correlated (near-zero).  Often, a comparison between two systems can be helpful to determine which residues are changing behavior and in what way.  The CDA script under *Corr_Decomp_Analysis.py* accepts two input files for comparison and outputs a single image file with four subplots.  These subplots show the following:

1. Change in correlated movement between two systems, both of which are still correlated and vary only by degree.
2. Change in anti-correlated movement between two systems, both of which are still anti-correlated and vary only by degree.
3. Systems which change from correlated to anti-correlated movement, and the degree of this change.
4. Systems which chanve from anti-correlated to correlated movement, and the degree of this change.

