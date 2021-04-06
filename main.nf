#!/usr/bin/env nextflow

/*
========================================================================================
                                G W A S - n f
========================================================================================
 Nextflow pipeline to run genome wide association studies using limix
 #### Homepage / Documentation
 https://gitlab.lrz.de/beckerlab/GWAS-nf
 #### Author
 Patrick Hüther <p.huether@lmu.de>
----------------------------------------------------------------------------------------
*/

// validate parameters
//ParameterChecks.checkParams(params)

Channel
    .fromPath(params.DMRs, checkIfExists: true)
    .splitCsv(sep: '\t')
    .unique()
    .map{ region -> tuple(region[0], region[1], region[2], region[-2]) }
    .set {ch_regions}

Channel
    .fromPath(params.arrowDataset, checkIfExists: true)
    .set {ch_arrow}

Channel
    .fromPath(params.genotype, checkIfExists: true)
    .set {ch_geno}

process retrieveRates {
    publishDir "${params.outdir}/traits", mode: 'copy'
    input:
        tuple val(chrom), val(start), val(end), val(context) from ch_regions
        path(dataset) from ch_arrow.collect()
    output:
        tuple val(chrom), val(start), val(end), val(context) , path('*.csv') into ch_pheno
    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import pyarrow.dataset as ds

        data = ds.dataset("${dataset}", format="arrow", partitioning="hive")
        
        reg = data.to_table(columns=['meth','unmeth','context','sample'], filter=(ds.field('chrom') == "${chrom}") & (ds.field('pos').isin(list(range(${start}, ${end}+1))))).to_pandas()

        reg = reg[reg['context'] == '${context}']
        reg['cov'] = reg['meth'] + reg['unmeth']
        reg['${chrom}_${start}_${end}'] = 100*(reg['meth']/(reg['meth'] + reg['unmeth']))

        reg = reg[reg['cov'] >= ${params.covthresh}]
        rates = reg.groupby(['sample'])['${chrom}_${start}_${end}'].mean()

        rates.to_csv("${context}_${chrom}_${start}_${end}.csv", na_rep='NA', float_format='%.3f', header=False)
        """
}

process filterGenotypes {
    tag "${chrom}_${start}_${end}"

    input:
        file geno from ch_geno.collect()
        tuple val(chrom), val(start), val(end), val(context), path(pheno) from ch_pheno
    output:
        tuple val(chrom), val(start), val(end), val(context), path('pheno.csv'), path('geno.pkl.xz'), path('kinship.pkl.xz') optional true into ch_filtered

    script:
        def kinship_mode = params.kinship_from_all_markers ? 'all' : 'filtered'
        """
        #!/usr/bin/env python

        import h5py
        import pandas as pd
        import numpy as np
        import logging

        from limix.stats import linear_kinship

        logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
        logger = logging.getLogger()

        pheno = pd.read_csv('${pheno}', index_col=[0], header=None).dropna() 
        pheno_acc_ids = np.array(pheno.index, dtype=np.uint32)

        # read SNP matrix
        with h5py.File('${geno}', 'r') as genofile:
            geno_acc_ids = np.array(genofile['accessions'][:], dtype=np.uint32)
            snps = np.array(genofile['snps'][:], dtype=np.bool)

            chr_names = genofile['positions'].attrs.get('chrs')
            chr_regions = np.array(genofile['positions'].attrs.get('chr_regions'))
            geno_chroms = []
            for i, reg in enumerate(chr_regions):
                geno_chroms.extend(np.repeat(chr_names[i].decode('utf-8'), reg[1]-reg[0]))
            pos = np.array(genofile['positions'][()], dtype=np.uint32)
            geno_chroms = np.array(geno_chroms)

        def get_kinship(snpmat, gower_norm):
            ibs = linear_kinship(snpmat.to_numpy().T)
            if gower_norm:
                from limix.qc import normalise_covariance
                ibs = normalise_covariance(ibs @ ibs.T)
            return pd.DataFrame(ibs, index=snpmat.columns, columns=snpmat.columns)
        
        genotypes = pd.DataFrame(snps,
                                 index=pd.MultiIndex.from_arrays([geno_chroms, pos]),
                                 columns=geno_acc_ids)

        if '${kinship_mode}' == 'all':
            kinship = get_kinship(genotypes, ${params.normalise_covariance.toString().capitalize()})

        pheno_geno_intersect = np.intersect1d(geno_acc_ids, pheno_acc_ids)
        
        phenotypes = pheno.loc[pheno_geno_intersect, :]
        genotypes = genotypes.loc[:, pheno_geno_intersect]
        
        logger.info('%i accessions with both genotype and phenotype. Removed %i accessions because of missing genotype.', len(phenotypes), len(pheno) - len(phenotypes))

        genotypes = genotypes[(genotypes.sum(axis=1) >= ${params.mac}) & (genotypes.sum(axis=1) <= genotypes.shape[1]-${params.mac})]

        logger.info('Removed SNPs below MAC threshold ${params.mac}. (Remaining SNPs: %i across %i accessions)', genotypes.shape[0], genotypes.shape[1])

        if '${kinship_mode}' == 'filtered':
            kinship = get_kinship(genotypes, ${params.normalise_covariance.toString().capitalize()})
        else:
            kinship = kinship.loc[genotypes.columns, genotypes.columns]

        if genotypes.shape[0] > 0:
            phenotypes.to_csv("pheno.csv", header=False)
            genotypes.to_pickle("geno.pkl.xz")
            kinship.to_pickle("kinship.pkl.xz")
        """
}

process runGWAS {
    tag "${chrom}_${start}_${end}"

    publishDir "${params.outdir}/pvals", mode: 'copy'
    input:
        tuple val(chrom), val(start), val(end), val(context), path(pheno), path(geno), path(kinship) from ch_filtered

    output:
        tuple val(chrom), val(start), val(end), val(context), path('*.csv') into ch_pvals mode flatten optional true

    script:
        def pheno_transform = params.transform == 'no_transformation' ? "" : ".apply(${params.transform}, raw=True)"
        def locus_fixed = params.locus ? "genotypes.xs( (${params.locus.tokenize(',')[0]},${params.locus.tokenize(',')[1]}), axis=0).to_numpy().ravel()" : "None"
        """
        #!/usr/bin/env python
        
        import pandas as pd
        import numpy as np
        import scipy.stats as stats

        from limix.qtl import scan
        from limix.qc import compute_maf, normalise_covariance, mean_standardize, quantile_gaussianize, boxcox
        from limix.stats import linear_kinship

        phenotypes = pd.read_csv('${pheno}', index_col=[0], dtype=np.float32, header=None)${pheno_transform}

        pheno = phenotypes.to_numpy(dtype=np.float32)

        genotypes = pd.read_pickle('${geno}')
        chromosomes = np.array(genotypes.index.get_level_values(0))
        positions = np.array(genotypes.index.get_level_values(1))

        geno = genotypes.to_numpy().T

        kinship = pd.read_pickle('${kinship}').to_numpy()

        # calculate maf and mac
        mafs = compute_maf(geno)
        macs = geno.sum(axis=0)

        freq = pd.DataFrame(data={'maf': np.array(mafs), 'mac': np.array(macs)},
                            index=pd.MultiIndex.from_arrays([chromosomes, positions]))

        st = scan(G=geno,
                Y=pheno,
                M=${locus_fixed},
                K=kinship,
                verbose=True)

        effsize = st.effsizes['h2'].loc[st.effsizes['h2']['effect_type'] == 'candidate']

        def phenotypic_variance_explained(beta, beta_se, mafs, n):
            '''Estimate phenotypic variance explained following Shim et al. (2015) https://doi.org/10.1371/journal.pone.0120758'''
            return (2 * beta**2 * mafs * (1 - mafs)) / (2 * beta**2 * mafs * (1 - mafs) + beta_se**2 * 2 * n * mafs * (1 - mafs))

        pve = phenotypic_variance_explained(effsize['effsize'].to_numpy(), effsize['effsize_se'].to_numpy(), mafs, pheno.shape[0])

        result = pd.DataFrame(data={'pv': st.stats['pv20'].to_numpy(), 'gve': effsize['effsize'].to_numpy(), 'pve': np.array(pve)},
                            index=pd.MultiIndex.from_arrays([chromosomes, positions]))


        chisq = stats.chi2(df=1)
        l = chisq.isf(result['pv'].median()) / chisq.median()

        #if (result['pv'].min() < 0.05/len(result)) & (0.9 < l < 1.1):
        if (result['pv'].min() < 0.05/len(result)):
            # only save when bonferroni threshold is passed, and lambda inflation is less than 10%
            #result['-log10pv'] = -np.log10(result['pv'])
            result = result.join(freq)
            result.to_csv("${context}_${chrom}_${start}_${end}_mac${params.mac}.csv", index_label=['chrom', 'pos'])
        """
}

process plotGWAS {
    tag "${chrom}_${start}_${end}"

    publishDir "${params.outdir}/plots", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("manhattan.png")) "manhattan/$filename"
            else if (filename.endsWith("qq.png")) "qq/$filename" }
    input:
        tuple val(chrom), val(start), val(end), val(context), path(pvals) from ch_pvals

    output:
        path('*png')

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import numpy as np
        import pylab as plt

        from limix.plot import manhattan, qqplot

        def correct_multiple_tests(pvals, alpha=0.05):
            c = len(pvals)
            bonferroni = alpha/c 
            for idx, value in enumerate(sorted(pvals)):
                benjamini_hochberg = ((idx + 1.0) / float(c)) * alpha
                if value > benjamini_hochberg:
                    break
            return bonferroni, benjamini_hochberg

        result = pd.read_csv('${pvals}')

        bf, bh = correct_multiple_tests(result['pv'])
        plt.figure(figsize=[15, 4])
        plt.title("${context}\\n${chrom}_${start}_${end}")
        manhattan(result,
                  colora='#AECF7B',
                  colorb='#09774D',
                  pts_kws=dict(markersize=10,
                               alpha=0.7))
        plt.axhline(-np.log10(bf), color='gray', linestyle='-', label='Bonferroni')
        plt.axhline(-np.log10(bh), color='gray', linestyle=':', label='Benjamini-Hochberg')
        plt.legend(bbox_to_anchor=(0,1), loc='lower left', ncol=2)
        plt.savefig("${pvals.baseName}_manhattan.png")
        plt.figure(figsize=[15, 4])
        plt.title("${context}\\n${chrom}_${start}_${end}")
        qqplot(result['pv'],
               band_kws=dict(color='#AECF7B',
                             alpha=0.5),
               pts_kws=dict(markersize=5,
                            markeredgecolor='black',
                            mew=0.5,
                            color='#09774D'))
        plt.savefig("${pvals.baseName}_qq.png")
        """
}
