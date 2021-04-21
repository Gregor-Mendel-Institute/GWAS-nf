class ParameterChecks {
  static void checkParams(params) {
    assert params.phenotype, "ERROR! A phenotype table has to be provided!"
    assert ['no_transformation', 'mean_standardize', 'quantile_gaussianize', 'boxcox'].contains(params.transform), "ERROR! Phenotype transformation has to be either no_transformation, mean_standardize, quantile_gaussianize or boxcox"
    assert params.genotype.endsWith('hdf5'), 'ERROR! SNP matrix in hdf5 format has to be provided'
    assert params.maf > 0 && params.maf < 1, 'ERROR! Minor allele frequence threshold has to be between zero (0%) and one (100%)'
    }
}