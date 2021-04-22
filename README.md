# Introduction

This is a Nextflow pipeline to run Genome Wide Association Studies (GWAS) in Arabidopsis.
It uses the [Limix framework](https://github.com/limix/limix) for linear mixed models and supports both single-trait and multi-trait analysis while correcting for population structure using an IBS (identity-by-state) kinship matrix.


# Summary

# Installation

The pipeline is bundled with a container image containing all necessary dependencies.
Therefore it is required to have Nextflow installed, along with one of the following Container engines: Docker, Podman, Singularity or Charliecloud.


To get the pipeline up an running, first clone the repository:

```shell
git clone git@gitlab.lrz.de:beckerlab/gwas-nf.git
```

> Note: The repository contains genotype information from the 1001genomes project (The 1001 Genomes Consortium, 2016), as well as the Regmap 250K panel (Horton et al., 2012). They are stored as `hdf5` files under `assets/genotypes`. If you wish to use these, it is required to have [`git-lfs`](https://git-lfs.github.com/) installed.


# Running the pipeline

To run the pipeline, it is required to provide at least the `--genotypes` and `--phenotypes` parameters.

## --genotypes

`hdf5` file containing the SNP information across all accessions as a binary matrix compared to the Col-0 reference allele. This has to be stored in a certain way inside a HDF5 file and is the main reason why this pipeline only works for Arabidopsis currently. In the future, I plan to add support for more standardized formats such as `PLINK` or `VCF`, so it can be used for other organisms as well without having to create a `hdf5` SNP-matrix first.

## --phenotypes

Plain CSV table containing *atleast* 2 columns. The first column has to contain the accession ID, which will be matched with the genotype information in the SNP-matrix. The remaining trait column(s) should contain the phenotype that was measured for the respective accession ID. (See `assets/phenotypes` for examples)

> Note: if there are more than 2 traits in the phenotype table, the pipeline will automatically split the table and submit parallel jobs for each phenotype.

If these 2 parameters are provided, the pipeline is ready to go (using default settings, test data within the repository and assuming that docker is running on the system):

```shell
nextflow gwas-nf/main.nf --phenotype gwas-nf/assets/phenotypes/FT10_mean.csv --genotype gwas-nf/assets/genotypes/regmap.hdf5 -profile docker
```

## Other Parameters

### --multitrait `(default: false)`

If specified, the pipeline runs in multitrait mode. In this case, *2* phenotype tables have to be provided and the path has to be enclosed in quotes. The pipeline will match them by the names of the trait columns.
Example:

```shell
nextflow gwas-nf/main.nf --phenotype "gwas-nf/assets/phenotypes/{FT10,FT16}_mean.csv" --genotype gwas-nf/assets/genotypes/regmap.hdf5 --multitrait -profile docker
```

### --maf `(default: 0.05)`

Minor allele frequency threshold for SNPs to be considered. Defaults to 5%.

### --transform `(default: 'mean_standardize')`

Transformation to be applied to the phenotype vector, valid options include 'no_transformation', 'mean_standardize', 'quantile_gaussianize' and 'boxcox'
Defaults to zscore scaling of the data to zero mean and unit variance.

### --trait `(default: false)`

Option to select a single column from a large phenotype table. Has to be a string matching the column name enclosed in quotes.

### --locus `(default: false)`

Option to specify a SNP location as a fixed effect. Needs to be passed as a comma-separated list of chromosome and position enclosed in quotes.
Example:

```shell
nextflow gwas-nf/main.nf --phenotype "gwas-nf/assets/phenotypes/FT10_mean.csv" --genotype gwas-nf/assets/genotypes/regmap.hdf5 --locus '3,2761231' -profile docker
```

### --outdir `(default: ./results)`

Directory where the results should be stored.

## Cluster specific profiles

## BioHPC Genomics cluster

The pipeline has a preconfigured profile to run on the BioHPC Genomics cluster using Charliecloud (requires Nextflow `21.03.0-edge` or later)

```shell
module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/x86_avx2/linux*
module load nextflow/21.03.0-edge charliecloud/0.22

nextflow gwas-nf/main.nf --phenotype "gwas-nf/assets/phenotypes/FT10_mean.csv" --genotype gwas-nf/assets/genotypes/regmap.hdf5 -profile biohpc_gen
```

## CBE cluster


The pipeline has a preconfigured profile to run on CBE using Singularity

```shell
module load nextflow/20.10.0

nextflow gwas-nf/main.nf --phenotype "gwas-nf/assets/phenotypes/FT10_mean.csv" --genotype gwas-nf/assets/genotypes/regmap.hdf5 -profile cbe
```
