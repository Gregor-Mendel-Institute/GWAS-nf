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
ParameterChecks.checkParams(params)

Channel
    .fromFilePairs(params.phenotype, checkIfExists: true, size: 1)
    .set {ch_pheno}

Channel
    .fromPath(params.genotype, checkIfExists: true)
    .set {ch_geno}

process scatterPhenotypes {
    tag "$env"

    echo true
    publishDir "${params.outdir}/traits", mode: 'copy'
    input:
        tuple val(env), path(pheno) from ch_pheno
    output:
        tuple val(env), path('*.csv') into traits mode flatten

    script:
        def selection = params.trait ? "['${params.trait.tokenize(',').join("','")}']" : "phenotype.columns"
        """
        #!/usr/bin/env python

        import pandas as pd

        phenotype = pd.read_csv("${pheno}", index_col=[0])
        for trait in ${selection}: 
            try:
                slice = phenotype[trait].dropna()
                assert slice.nunique() > 1
            except KeyError:
                print(f'Trait {trait} not found in ${pheno.name}. Skipping.')
            except AssertionError:
                print(f'Trait values for {trait} are all equal in ${pheno.name}. Skipping.')
            else:
                slice.sort_index().to_csv(f'{trait}.csv', header=False)
        """
}

traits
 .map { env, file -> [ env, file.baseName, file] }
 .groupTuple(by: params.multitrait ? 1 : 2, size: params.multitrait ? 0 : 1)
 .set {ch_traitsplit}

process filterGenotypes {
    tag "$traitname"

    input:
        path(geno) from ch_geno.collect()
        tuple val(env), val(traitname), path(traitfile, stageAs: 'trait*.csv') from ch_traitsplit
    output:
        tuple val(env), val(traitname), path('pheno.csv'), path('geno.pkl.xz'), path('kinship.pkl.xz') into ch_filtered

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

        pheno = pd.concat(map(lambda file: pd.read_csv(file, index_col=[0]), ['${traitfile.join("','")}']), axis=1).dropna()
        
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

        phenotypes.to_csv("pheno.csv")
        genotypes.to_pickle("geno.pkl.xz")
        kinship.to_pickle("kinship.pkl.xz")
        """
}

process runGWAS {
    tag "$traitname"

    publishDir "${params.outdir}/pvals", mode: 'copy'
    input:
        tuple val(env), val(traitname), path(pheno), path(geno), path(kinship) from ch_filtered

    output:
        tuple val(env), val(traitname), path('*.csv') into ch_pvals mode flatten optional true

    script:
        def pheno_transform = params.transform == 'no_transformation' ? "" : ".apply(${params.transform}, raw=True)"
        def locus_fixed = params.locus ? "genotypes.xs( (${params.locus.tokenize(',')[0]},${params.locus.tokenize(',')[1]}), axis=0).to_numpy().ravel()" : "None"
        if (!params.multitrait)
            """
            #!/usr/bin/env python
            
            import pandas as pd
            import numpy as np

            from limix.qtl import scan
            from limix.qc import compute_maf, mean_standardize, quantile_gaussianize, boxcox

            phenotypes = pd.read_csv('${pheno}', index_col=[0])${pheno_transform}

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

            if result['pv'].min() < ${params.pthresh}: 
                #result['-log10pv'] = -np.log10(result['pv'])
                result = result.join(freq)
                result.to_csv("${env}_${traitname}_mac${params.mac}.csv", index_label=['chrom', 'pos'])
            """
        else
            """
            #!/usr/bin/env python

            import pandas as pd
            import numpy as np
            
            from limix.qtl import scan
            from limix.qc import compute_maf, mean_standardize, quantile_gaussianize, boxcox

            phenotypes = pd.read_csv('${pheno}', index_col=[0])${pheno_transform}
            
            pheno = phenotypes.to_numpy(dtype=np.float32)
            
            genotypes = pd.read_pickle('${geno}')

            geno = genotypes.to_numpy().T

            kinship = pd.read_pickle('${kinship}')
            
            # calculate maf and mac
            mafs = compute_maf(geno)
            macs = geno.sum(axis=0)

            chromosomes = np.array(genotypes.index.get_level_values(0))
            positions = np.array(genotypes.index.get_level_values(1))

            freq = pd.DataFrame(data={'maf': np.array(mafs), 'mac': np.array(macs)},
                                index=pd.MultiIndex.from_arrays([chromosomes, positions]))

            n_pheno = pheno.shape[1]  # number of traits

            A = np.eye(n_pheno)  # p x p matrix of fixed effect sizes
            # common effects: 1 DoF
            Asnps0 = np.ones((n_pheno, 1))
            Asnps = np.eye(n_pheno)

            mtlmm = scan(G=geno,
                        Y=pheno,
                        K=kinship,
                        A=A,
                        A0=Asnps0,
                        A1=Asnps,
                        verbose=True)

            # specific (GxE)
            specific = pd.DataFrame(mtlmm.stats['pv21'].to_numpy(),
                                    index=pd.MultiIndex.from_arrays([chromosomes, positions]),
                                    columns=['pv'])                      

            # common (G)
            common = pd.DataFrame(mtlmm.stats['pv10'].to_numpy(),
                                  index=pd.MultiIndex.from_arrays([chromosomes, positions]),
                                  columns=['pv'])

            # common (G + GxE)
            any = pd.DataFrame(mtlmm.stats['pv20'].to_numpy(),
                               index=pd.MultiIndex.from_arrays([chromosomes, positions]),
                               columns=['pv'])

            results =  {'specific': specific, 'common': common, 'any': any}

            for name, result in results.items():
                if result['pv'].min() < ${params.pthresh}:
                    #result['-log10pv'] = -np.log10(result['pv'])
                    result = result.join(freq)
                    result.to_csv(f'${traitname}_mac${params.mac}_{name}.csv', index_label=['chrom', 'pos'])
            """
}

process plotGWAS {
    tag "$traitname"

    publishDir "${params.outdir}/plots", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("manhattan.png")) "manhattan/$filename"
            else if (filename.endsWith("qq.png")) "qq/$filename" }
    input:
        tuple val(env), val(traitname), path(pvals) from ch_pvals

    output:
        tuple val(env), val(traitname), path('*manhattan.png'), path('*qq.png') into ch_plots

    script:
        def effect = params.multitrait ? pvals.baseName.tokenize('_')[-1] : env
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
        plt.title("${effect}\\n${traitname}")
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
        plt.title("${effect}\\n${traitname}")
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
