import os
import pandas
import numpy
import subprocess as sp
import glob
import re

BASE_PATH = "/lab/data/gwas/2018_06_NHGRI_EBI_Catalog_download"

DATA = {
    'nhgri_gwas': "/lab/data/gwas/2018_06_NHGRI_EBI_Catalog_download/gwas_catalog_v1.0.2-associations_e92_r2018-05-29.tsv",
    '1000g' : "/lab/data/genomes/human/hg19/1000GenomesDownloads/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    '1000g_index' : "/lab/data/genomes/human/hg19/1000GenomesDownloads/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi",
    'vcf_sampleInfo' : "/lab/data/genomes/human/hg19/1000GenomesDownloads/igsr_3990samples_allPopulations.tsv",
}

DIRECTORIES = {
    'intermediateFiles' : "r2_0.2_EUR_pruning"
}

SCRIPTS = {
    'gwas_nhgri' : os.path.join(BASE_PATH, "configuration_for_pruning/scripts/get_GWAS.py"),
}

CHROM = list(range(1, 23))


rule final:
    """
    NHGRI GWAS 
    1. Prune by LD r2<0.2
    2. Run GREGOR enrichment
    """
    input:
        main = dynamic(os.path.join(DIRECTORIES['intermediateFiles'], "pruned.{trait}.txt"))
        

rule makeFullSnpList_forPruning:
    """Make full SNP list to subset vcf files once """
    input:
        nhgri_gwas = DATA['nhgri_gwas']
    output:
        main = os.path.join(DIRECTORIES['intermediateFiles'], "traits_pruning", "full_snp_list.txt"),
        full = os.path.join(DIRECTORIES['intermediateFiles'], "traits_pruning", "nhgri_gwas.dat"),
        reference_table = os.path.join(DIRECTORIES['intermediateFiles'], "traits_pruning", "gwas_reference.xlsx"),
    params:
        script = SCRIPTS['gwas_nhgri'],
        ancestry = "European",
        minN = 30,
    shell:
        ' python {params.script} '
        ' --gwasfile {input.nhgri_gwas} '
        ' --ancestry {params.ancestry} '
        ' --minN {params.minN} '
        ' --output_full {output.full} '
        ' --output_snpfile {output.main} '
        ' --output_reference_table {output.reference_rable} ' 


rule make_trait_files:
    input:
        full = rules.makeFullSnpList_forPruning.output.full,
    output:
        main = dynamic(os.path.join(DIRECTORIES['intermediateFiles'], "traits_pruning", "selected.{trait}.dat"))
    run:
        d_selected_traits = pandas.read_csv(input.full, sep='\t')
        for trait, group in d_selected_traits.groupby('trait'):
            filename = os.path.join(DIRECTORIES['intermediateFiles'], "traits_pruning", f"selected.{trait}.dat")
            group[['SNP', 'P']].to_csv(filename, sep='\t', index=False, na_rep="NA")

        
rule makeSampleFile:
    """Get samples according to population codes to subest 1000g vcf files for pruning"""
    input:
        sampleInfo = DATA['vcf_sampleInfo']
    output:
        samplefile = os.path.join(config['output_directory'], "subsetSamples.txt")
    run:
        d = pandas.read_csv(input.sampleInfo, sep='\t')
        d = d[d[config['population_type']] == config['population_code']]
        d[['Sample name']].to_csv(output[0], header=False, index=False)

        
rule subsetVCF:
    """Subset vcf for samples, remove indels and SNPs which are to be pruned among themselves."""
    input:
        snpfile = DATA['1000g'],
        samplefile = rules.makeSampleFile.output.samplefile,
        posfile = rules.makeFullSnpList_forPruning.output.main, #rules.setup_eqtl_for_pruning.output.snplist
    output:
        vcf = temp(os.path.join(config['output_directory'], "chr{chrom}.selected.recode.vcf.gz")),
        index = temp(os.path.join(config['output_directory'], "chr{chrom}.selected.recode.vcf.gz.tbi")),
    params:
        outstring = os.path.join(config['output_directory'], "chr{chrom}.selected"),
    shell:
        r"""
        vcftools --gzvcf {input.snpfile}  --keep {input.samplefile} \
        --remove-indels \
        --snps {input.posfile} \
        --out {params.outstring} --recode ;
        bgzip {params.outstring}.recode.vcf
        tabix {output.vcf}
        """

        
rule get_plink_files:
    """
    Make plink format input files after filtering 1000g vcf. IMP - vcf files are usually large so designated to be temp
    plink map and ped files only contain bi-allelic loci. 
    """
    input:
        vcf = expand(os.path.join(config['output_directory'], "chr{chrom}.selected.recode.vcf.gz"), chrom = CHROM),
        index = expand(os.path.join(config['output_directory'], "chr{chrom}.selected.recode.vcf.gz.tbi"), chrom = CHROM)
    output:
        vcf = temp(os.path.join(config['output_directory'], "myfile.selected.recode.vcf.gz")),
        mapfile = temp(os.path.join(config['output_directory'], "myfile.selected.map")),
        pedfile = temp(os.path.join(config['output_directory'], "myfile.selected.ped")),
    params:
        outstring = os.path.join(config['output_directory'], "myfile.selected"),
    shell:
        r"""
        vcf-concat {input.vcf} | bgzip -c > {output.vcf} ;
        vcftools --gzvcf {output.vcf} --plink --out {params.outstring}
        """

        
rule prune_plink:
    """
    Prune a list od SNPs using P value of association, using 1000g phase 3 vcf. Population codes or Superpopulation codes can be used to subset 1000g samples  
    Vcf files first subset by the selected population. All indels are removed. Also, plink files only contain bialleic SNPs. Pruning is done using --clump flags
    """
    input:
        mapfile = rules.get_plink_files.output.mapfile,
        pedfile = rules.get_plink_files.output.pedfile,
        inputfile = os.path.join(DIRECTORIES['intermediateFiles'], "traits_pruning", "selected.{trait}.dat") #rules.setup_eqtl_for_pruning.output.full
    output:
        clumpedfile = temp(os.path.join(config['output_directory'], "{trait}.selected.clumped")),
    params:
        instring = os.path.join(config['output_directory'], "myfile.selected"),#rules.get_plink_files.params.outstring,
        outstring =  os.path.join(config['output_directory'], "{trait}.selected"),
        r2 = config['prune_r2'],
        p1 = .99,
        p2 = .99,
    shell:
        r"""
        /lab/sw/modules/plink/1.9/bin/plink --file {params.instring} \
        --clump {input.inputfile} --clump-r2 {params.r2} --clump-p1 {params.p1} --clump-p2 {params.p2} \
        --out {params.outstring}
        """

        
rule organize_pruned_results:
    input:
        pruned = rules.prune_plink.output.clumpedfile,
    output:
        main = os.path.join(DIRECTORIES['intermediateFiles'], "pruned.{trait}.txt"),
    run:
        dpruned = pandas.read_csv(input.pruned, delim_whitespace=True)
        dpruned.loc[:,'out'] = dpruned.apply(lambda x: "chr{chrom}:{pos}".format(chrom = x['CHR'], pos = x['BP']), axis=1)
        print(output.main)
        dpruned[['out']].to_csv(output.main, index=False, header=False)
            

