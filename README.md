# 17KFishTree

This contains all data files, scripts, pipeline documentation, supplementary files of the manuscript "Unlocking the Evolutionary Relationships of 17,407 Ray-finned Fish Species: A Nested Tree Building Approach"

Re-usable nested-tree-pipeline-------

Phylogenetic Tree Construction Pipeline

Overview

This pipeline guides the preparation of datasets and the construction of a phylogenetic tree using a nested approach. It involves the following steps:
1.	Dataset Preparation: Generate genome, multi-gene, and placement sequence datasets using R scripts.
2.	Genome Tree Construction: Build a genome-based phylogenetic tree using RAxML-NG, IQ-TREE, or any preferred tree-building tool.
3.	Multi-Gene Tree Construction: Construct a multi-gene tree using the genome tree as a constraint.
4.	Placement of COI Barcode Sequences: Extend the tree by placing COI barcode sequences onto the multi-gene tree.
---------------------------------------
Step 1: Dataset Preparation

Execute the three R scripts to generate the necessary datasets:

1.1 Generate Genome Dataset
R script genome_dataset.R

1.2 Generate Multi-Gene Dataset
R script multigene_dataset.R

1.3 Generate Placement Sequences Dataset
R script placement_dataset.R

Ensure that each script outputs properly formatted sequence alignments and metadata files for subsequent analysis.

Output datasets: 
Genome dataset: genomeNew.fasta
Multi-gene dataset: FinalAlignmentChap2_12193New.fasta
Placement dataset: df_COI_5218_new.fasta

---------------------------------------
Step 2: Genome Tree Construction

Use RAxML-NG (IQ-TREE, or another tree-building tool) to generate the genome tree.

Using RAxML-NG
Note: The genome dataset is too large to be uploaded to GitHub. If you need access to the dataset, please contact me at thanuja@uoguelph.ca.

srun raxml-ng-mpi --data-type AA -seed 2 --msa genomeNew.fasta --model JTT+I+G4+F --outgroup Erpetoichthys_calabaricus,Erpetoichthys_calabaricus --tree pars{10},rand{10} --prefix GenomeTree

This step produces genome_tree.treefile, which serves as the constraint tree for the next step.

---------------------------------------
Step 3: Multi-Gene Tree Construction

Use the genome tree as a constraint to build the multi-gene tree.

Using RAxML-NG

srun raxml-ng-mpi -seed 35 --msa FinalAlignmentChap2_12193New.fasta --tree-constraint GenomeTree.nwk --threads 9 --force  perf_threads --model FinalAlignmentChap2_12193New.fasta.part.bic --outgroup Gallus_gallus,Homo_sapiens,Callorhinchus_milii --tree pars{10},rand{10} --prefix MultiGeneTree_12193tips

This results in multigene_tree.treefile, which will be used as the backbone tree for placement.

---------------------------------------
Step 4: Placement of COI Barcode Sequences

Use EPA-ng or another placement tool to extend the tree by placing COI barcode sequences.

Run EPA-ng for Placement

First, RAxML-ng evaluates the multi-gene reference tree.
raxml-ng --evaluate --msa FinalAlignmentCOI_ONLY_Chap2_12193New.fasta --tree MultiGeneTree_12193tips.nwk --prefix infochap2NewCOIONLY --model TPM2uf+I+G4

Then, EPA-ng places the COI sequences onto this reference tree, extending it.
srun epa-ng --ref-msa FinalAlignmentChap2_12193New.fasta --tree multigene_tree.nwk --query df_COI_5218_new.fasta --model TPM2uf{0.756132/8.888316/1.000000}+FU{0.356271/0.267503/0.117907/0.258319}+IU{0.164270}+G4m{0.401973} --preserve-rooting on --outdir PlacementTree_17410tips

The final output is a phylogenetic tree incorporating both multi-gene and COI-based placements.

---------------------------------------
Final Output
•	genome_tree.nwk: Genome-based tree
•	multigene_tree.nwk: Multi-gene tree constrained by genome tree
•	Placement_tree.nwk: Final tree with COI barcode sequence placements
This pipeline ensures a reproducible workflow for constructing phylogenetic trees using a nested approach. Modify parameters as needed based on dataset size and computational resources.
