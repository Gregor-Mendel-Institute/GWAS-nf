#!/usr/bin/env nextflow

/*
========================================================================================
                                G W A S - n f
========================================================================================
 Nextflow pipeline to run genome wide association studies using limix
 #### Homepage / Documentation
 https://github.com/Gregor-Mendel-Institute/GWAS-nf
 #### Author
 Patrick Hüther <patrick.huether@gmi.oeaw.ac.at>
----------------------------------------------------------------------------------------
*/

Channel
    .fromFilePairs(params.phenotype, checkIfExists: true, size: 1)
    .set {ch_pheno}

Channel
    .fromPath(params.genotype, checkIfExists: true)
    .set {ch_geno}

Channel
    .fromPath(params.kinship, checkIfExists: true)
    .set {ch_kinship}

process scatterPhenotypes {
    publishDir "${params.outdir}/traits", mode: 'copy'
    input:
        tuple val(env), file(pheno) from ch_pheno

    output:
        tuple val(env), file('*.csv') into traits mode flatten

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        
        phenotype = pd.read_csv("${pheno}", index_col=[0])
        for index, trait in enumerate(phenotype.columns):
            chunk = phenotype.iloc[:,[index]].dropna()
            chunk.to_csv("{}.csv".format(trait))
        """
}

traits
 .map { env, file -> [ env, file.baseName, file] }
 .groupTuple(by: 1)
 .set {ch_traitsplit}

process filterGenotypes {
    publishDir "${params.outdir}/genotypes", mode: 'copy'
    input:
        file geno from ch_geno.collect()
        file kinship from ch_kinship.collect()
        tuple val(env), val(traitname), file(traitfile:"?/*") from ch_traitsplit

    output:
        tuple val(traitname), file('pheno.csv'), file('geno.pkl.xz'), file('kinship.pkl.xz') into ch_filtered

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import numpy as np

        pheno = pd.concat(map(lambda file: pd.read_csv(file, index_col=[0]), ['${traitfile.join("','")}']), axis=1, join='inner')
        geno = pd.read_pickle('${geno}')

        pheno_geno_intersect, strain_idx, pheno_idx = np.intersect1d(np.array(geno.columns), 
                                                                     np.array(pheno.index),
                                                                     return_indices=True)

        accs_no_geno_info = np.array(pheno.index)[np.in1d(pheno.index, pheno_geno_intersect, invert=True)]
        print("No genotype information for accessions: {}".format(accs_no_geno_info))

        phenotypes = pheno.drop(accs_no_geno_info)

        print("Removed them from list of phenotypes, leaving {} accessions.".format(len(phenotypes)))

        geno = geno.filter(items=pheno_geno_intersect, axis=1)
        geno = geno[(geno.sum(axis=1) >= ${params.mac}) & (geno.sum(axis=1) <= geno.shape[1]-${params.mac})]
        genotypes = geno.reindex(sorted(geno.columns), axis=1)

        print("Removed {:d} snps because of MAC threshold ${params.mac}. (Remaining snps: {:d} across {:d} accessions)"
        .format(geno.shape[0]-genotypes.shape[0], genotypes.shape[0], genotypes.shape[1]))

        phenotypes.to_csv("pheno.csv")
        genotypes.to_pickle("geno.pkl.xz")

        geno_ids = np.array(genotypes.columns.astype(int))
        kinship = pd.read_pickle('${kinship}')
        kinship_filtered = kinship.loc[geno_ids,geno_ids]

        kinship_filtered.to_pickle("kinship.pkl.xz")
        """
}

process runGWAS {
    publishDir "${params.outdir}/pvals", mode: 'copy'
    input:
        tuple val(traitname), file(pheno), file(geno), file(kinship) from ch_filtered

    output:
        tuple val(traitname), file('*.csv') into ch_pvals mode flatten

    script:
        if (params.singletrait)
            """
            #!/usr/bin/env python

            import pandas as pd
            import numpy as np
            
            from limix.qtl import scan
            from limix.qc import compute_maf

            phenotypes = pd.read_csv('${pheno}')

            genotypes = pd.read_pickle('${geno}')
            kinship = pd.read_pickle('${kinship}')

            pheno_norm = phenotypes.values.astype(float)
            p1 = pheno_norm[:, 0]
            p1 = (p1 - p1.mean()) / p1.std()

            # calculate maf and mac
            mafs = compute_maf(genotypes.values.T)
            macs = genotypes.sum(axis=1)

            chromosomes = np.array(genotypes.index.get_level_values(0))
            positions = np.array(genotypes.index.get_level_values(1))

            freq = pd.DataFrame(data={'maf': np.array(mafs), 'mac': np.array(macs)},
                                index=pd.MultiIndex.from_arrays([chromosomes, positions],
                                names=['chrom','pos']))

            result = scan(G=genotypes.values.T,
                                    Y=p1,
                                    K=kinship.values,
                                    verbose=True)

            result = pd.DataFrame(result.stats['pv20'].values,
                                index=pd.MultiIndex.from_arrays([chromosomes, positions]),
                                names=['chrom','pos'],
                                columns=['pv']).join(freq)
            
            result.to_csv("mac${params.mac}_${traitname}_singletrait.csv")
            """
        else
            """
            #!/usr/bin/env python

            import pandas as pd
            import numpy as np
            
            from limix.qtl import scan
            from limix.qc import compute_maf

            phenotypes = pd.read_csv('${pheno}')

            genotypes = pd.read_pickle('${geno}')
            kinship = pd.read_pickle('${kinship}')

            pheno_norm = phenotypes.values.astype(float)
            p1 = pheno_norm[:, 0]
            p2 = pheno_norm[:, 1]
            p1 = (p1 - p1.mean()) / p1.std()
            p2 = (p2 - p2.mean()) / p2.std()
            pheno_norm = np.vstack([p1, p2]).T

            n_pheno = pheno_norm.shape[1]  # number of traits

            # calculate maf and mac
            mafs = compute_maf(genotypes.values.T)
            macs = genotypes.sum(axis=1)

            chromosomes = np.array(genotypes.index.get_level_values(0))
            positions = np.array(genotypes.index.get_level_values(1))

            freq = pd.DataFrame(data={'maf': np.array(mafs), 'mac': np.array(macs)},
                                index=pd.MultiIndex.from_arrays([chromosomes, positions],
                                names=['chrom','pos']))

            A = np.eye(n_pheno)  # p x p matrix of fixed effect sizes
            # common effects: 1 DoF
            Asnps0 = np.ones((n_pheno, 1))
            Asnps = np.eye(n_pheno)

            print("Calculating mtlmm ... ")
            mtlmm = scan(G=genotypes.values.T,
                        Y=pheno_norm,
                        K=kinship.values,
                        A=A,
                        A0=Asnps0,
                        A1=Asnps,
                        verbose=True)

            # specific (GxE)
            specific = pd.DataFrame(mtlmm.stats['pv21'].values,
                                    index=pd.MultiIndex.from_arrays([chromosomes, positions],
                                    names=['chrom','pos']),
                                    columns=['pv']).join(freq)
            # common (G)
            common = pd.DataFrame(mtlmm.stats['pv10'].values,
                                index=pd.MultiIndex.from_arrays([chromosomes, positions],
                                names=['chrom','pos']),
                                columns=['pv']).join(freq)
            # common (G + GxE)
            any = pd.DataFrame(mtlmm.stats['pv20'].values,
                            index=pd.MultiIndex.from_arrays([chromosomes, positions],
                            names=['chrom','pos']),
                            columns=['pv']).join(freq)

            results =  {'specific': specific, 'common': common, 'any': any}

            for name, result in results.items():
                result.to_csv("mac${params.mac}_${traitname}_{}.csv".format(name))
            """
}

process plotGWAS {
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:
        tuple val(traitname), file(pvals) from ch_pvals

    output:
        tuple val(traitname), file('*manhattan.png'), file('*qq.png') into ch_plots

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import numpy as np
        import pylab as plt

        from limix.plot import manhattan, qqplot

        def correct_multiple_tests(pvals, alpha=${params.alpha}):
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
        plt.title("${traitname} - ${pvals.baseName}")
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
        plt.title("${traitname} - ${pvals.baseName}")
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