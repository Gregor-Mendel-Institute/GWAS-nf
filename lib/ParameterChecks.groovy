class ParameterChecks {
  static void checkParams(params) {
    assert params.phenotype, "ERROR! A phenotype table has to be provided!"
    assert ['no_transformation', 'mean_standardize', 'quantile_gaussianize', 'boxcox'].contains(params.transform), "ERROR! Phenotype transformation has to be either no_transformation, mean_standardize, quantile_gaussianize or boxcox"
    assert params.genotype.contains('hdf5'), 'ERROR! SNP matrix in hdf5 format has to be provided'
    assert params.mac instanceof Integer, 'ERROR! Minor allele count threshold has to be an Integer'
    }
}