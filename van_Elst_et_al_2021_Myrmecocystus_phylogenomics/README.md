# Comprehensive phylogeny of *Myrmecocystus* honey ants highlights cryptic diversity and infers evolution during aridification of the American Southwest

This repository holds scripts for the following analyses conducted for the publication [van Elst et al. (2021), *Mol. Phylogenet. Evol.*](https://doi.org/10.1016/j.ympev.2020.107036):
- Read cleaning
- Assembly
- Extraction and harvesting of UCE sequences
- Alignment and trimming
- Partitioning
- Phylogenetic inference
- Divergence time estimation
- Ancestral state inference

Each directory contains a script from which the respective pipeline is executed, labeled by the suffic "_sub". Input and output files can be found in the [Zenodo archive](https://doi.org/10.5281/zenodo.4061988). 

### Read cleaning
`./read_cleaning` contains scripts to clean demultiplexed raw sequencing reads with [Illumiprocessor](https://illumiprocessor.readthedocs.io/en/latest/). 

### Assembly
`./assembly` contains scripts to assemble clean sequencing reads with [metaSPAdes v3.13.1] (https://github.com/ablab/spades).

### Extraction and harvesting of UCE sequences
`./uces` contains scripts to harvest additional UCE sequences form published genomes and to extract UCEs for a specified sample set from assembled sequences with [phyluce v1.6.7](https://phyluce.readthedocs.io/en/latest/).

### Alignment and trimming
`./alignment_trimming` contains scripts to align, trim and concatenate single UCE loci with [MAFFT v7.429](https://mafft.cbrc.jp/alignment/software/), [trimAlv1.4.rev15](http://trimal.cgenomics.org/trimal)] and [AMAS](https://github.com/marekborowiec/AMAS), respectively.

### Partitioning
`./partitioning` contains scripts to create and optimize alternative partitioning schemes for the concatenated UCE alignment with [PartitionFinder v2.1.1](https://www.robertlanfear.com/partitionfinder/).

### Phylogenetic inference
`./phylogenetic_inference` contains scripts to infer maximum likelihood phylogenies with [IQ-TREE v1.6.11](http://www.iqtree.org/) as well as species trees with [ASTRAL-III v5.6.3](https://github.com/smirarab/ASTRAL) and [SVDquartets of PAUP* v4.0a](https://paup.phylosolutions.com/).

### Divergence time estimation
`./divergence_times` contains scripts to infer divergence times on a reduced sample set with [MCMCTree of PAML v4.8](http://abacus.gene.ucl.ac.uk/software/paml.html).

### Ancestral state inference
`./ancestral_states` contains scripts to estimate ancestral states of foraging behavior with SIMMAP of the R package [phytools v0.6–99](https://github.com/liamrevell/phytools).