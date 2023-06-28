# Diversification processes in Gerp's mouse lemur demonstrate the importance of rivers and altitude as biogeographic barriers in Madagascar's humid rainforests

This repository holds scripts for the following analyses conducted for the publication [van Elst et al. (2023), *Ecol. Evol.*](https://doi.org/10.1002/ece3.10254):
- Read trimming
- Reference mapping
- Genotype calling
- VCF filtering
- Genotype likelihood estimation
- Locus extraction
- Phylogenetic inference
- Population structure
- Coalescent modelling
- Ecological niche modelling

Each directory contains a script from which the respective pipeline is executed, labeled by the suffix "_sub". Input and output files can be found at [Dryad](https://doi.org/10.5061/dryad.9w0vt4bmr). 

### Read trimming
`./read_trimming` contains scripts to trim demultiplexed raw sequencing reads with [Trimmomatic v0.39](https://github.com/usadellab/Trimmomatic). 

### Reference mapping
`./reference_mapping` contains scripts to align cleaned reads against the *Microcebus murinus* reference genome (Mmur 3.0; [Larsen et al. (2017), *BMC Biol.*](https://doi.org/10.1186/s12915-017-0439-6) with [BWA v7.0.17](https://github.com/lh3/bwa) and to filter BAM files with [SAMtools v1.11](http://www.htslib.org/).

### Genotype calling
`./genotype_calling` contains scripts to call genotypes with [GATK v4.1.9.0](https://gatk.broadinstitute.org/hc/en-us).

### VCF filtering
`./vcf_filtering` contains scripts to apply various filters to VCF files with [VCFtools v0.1.17](https://vcftools.github.io/index.html) and [GATK v3.8.1/v4.1.9.0](https://gatk.broadinstitute.org/hc/en-us).

### Genotype likelihood estimation
`./genotype_likelihoods` contains scripts to infer genotype and site allele frequency likelihoods with [ANGSD v0.934](http://www.popgen.dk/angsd/index.php/ANGSD).

### Locus extraction
`./locus_extraction` contains scripts to extract FASTA sequences of RAD loci based on filtered genotype calls with [BEDtools v2.30.0](https://bedtools.readthedocs.io/en/latest/) and [SAMtools v1.11](http://www.htslib.org/).

### Phylogenetic inference
`./phylogenetic_inference` contains scripts to infer maximum likelihood and quartet-based phylogenies with [RAxML-NG v1.0.2](https://github.com/amkozlov/raxml-ng) and [SVDquartets of PAUP* v4.0a](https://paup.phylosolutions.com/), respectively, and to infer a split network with the NeighborNet method in [SplitsTree v4.17.1](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/).

### Population structure inference
`./population_structure` contains scripts for the following population structure analyses:
- Principal component analysis with [PCAngsd v1.01](https://github.com/Rosemeis/pcangsd) from genotype likelihoods
- Ancestry inference with [NGSadmix v32](http://www.popgen.dk/software/index.php/NgsAdmix) from genotype likelihoods
- Testing for isolation-by-distance based on *F<sub>ST</sub>* values estimated with [realSFS](http://www.popgen.dk/angsd/index.php/RealSFS) from genotype likelihoods and on mean genetic distances between populations estimated with the R package ['vcfR' v1.12](https://github.com/knausb/vcfR) from genotype calls.
- Estimation of effective mirgration surfaces with [EEMS](https://github.com/dipetkov/eems) from genotype calls

### Coalescent modeling
`./coalescent_modeling` contains scripts for coalescent modelling with [G-PhoCS](http://compgen.cshl.edu/GPhoCS/).

### Ecological niche modeling
`./ecological_niche_modeling` contains scripts to model ecological niches with the Maxent algorithm and a random forest model in the R packages ['ENMTools' v1.0.6](https://github.com/danlwarren/ENMTools) and ['biomod2' v3.5.1](https://biomodhub.github.io/biomod2/), respectively.