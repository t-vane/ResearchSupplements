# Population genomic structure in Goodman's mouse lemur reveals long-standing separation of Madagascar's Central Highlands and eastern rainforests

This repository holds scripts for the following analyses conducted for the publication [Tiley et al. (2022), *Mol. Ecol.*](https://doi.org/10.1111/mec.16632):
- Read trimming
- Reference mapping
- Genotype calling
- VCF filtering
- Genotype likelihood estimation
- Phylogenetic inference
- Population structure inference
- Demographic modeling

Each directory contains a script from which the respective pipeline is executed, labeled by the suffic "_sub". Input and output files can be found in the [Dryad digital repository](https://doi.org/10.5061/dryad.dncjsxkvz). 

### Read trimming
`./read_trimming` contains scripts to trim demultiplexed raw sequencing reads with [Trimmomatic v0.39](https://github.com/usadellab/Trimmomatic). 

### Reference mapping
`./reference_mapping` contains scripts to align cleaned reads against the *Microcebus murinus* reference genome (Mmur 3.0; [Larsen et al. (2017), *BMC Biol.*](https://doi.org/10.1186/s12915-017-0439-6) with [BWA v7.0.17](https://github.com/lh3/bwa) and to filter BAM files with [SAMtools v1.11](http://www.htslib.org/).

### Genotype calling
`./genotype_calling` contains scripts to call genotypes with [Stacks v2.53](http://catchenlab.life.illinois.edu/stacks/).

### VCF filtering
`./vcf_filtering` contains scripts to apply various filters to VCF files with [VCFtools v0.1.17](https://vcftools.github.io/index.html) and [GATK v3.8.1/v4.1.9.0](https://gatk.broadinstitute.org/hc/en-us).

### Genotype likelihood estimation
`./genotype_likelihoods` contains scripts to infer genotype and site allele frequency likelihoods with [ANGSD v0.934](http://www.popgen.dk/angsd/index.php/ANGSD).

### Phylogenetic inference
`./phylogenetic_inference` contains scripts to infer maximum likelihood and quartet-based phylogenies with [RAxML-NG v1.0.2](https://github.com/amkozlov/raxml-ng) and [SVDquartets of PAUP* v4.0a](https://paup.phylosolutions.com/), respectively, as well as to conduct approximately unbiased tests with [IQ-TREE v2.2.0](http://www.iqtree.org/).

### Population structure inference
`./population_structure` contains scripts for the following population structure analyses:
- Principal component analysis with the R package ['adegenet' v2.1.3](https://cran.r-project.org/web/packages/adegenet/index.html) from genotype calls and with [PCAngsd v1.01](https://github.com/Rosemeis/pcangsd) from genotype likelihoods
- Ancestry inference with [Structure v2.3.4](https://web.stanford.edu/group/pritchardlab/structure.html) from genotype calls and with [NGSadmix v32](http://www.popgen.dk/software/index.php/NgsAdmix) from genotype likelihoods
- Testing for isolation-by-distance based on *F<sub>ST</sub>* values estimated with the R package ['hierfstat' v0.5-7](https://cran.r-project.org/web/packages/hierfstat/index.html) from genotype calls and with [realSFS](http://www.popgen.dk/angsd/index.php/RealSFS) from genotype likelihoods
- Analysis of molecular variance (AMOVA) with the R package ['poppr' v2.8.7](https://cran.r-project.org/web/packages/poppr/index.html) from genotype calls

### Demographic modeling
`./demographic_modeling` contains scripts to model demographic histories with [fastsimcoal v2.7](http://cmpg.unibe.ch/software/fastsimcoal27/) and to infer population size changes through time with [Stairway Plot v2.0](https://github.com/xiaoming-liu/stairway-plot-v2).
