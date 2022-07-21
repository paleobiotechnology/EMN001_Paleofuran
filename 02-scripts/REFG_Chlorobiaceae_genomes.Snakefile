################################################################################
# Project: Natural products from the Palaeolithic
# Part: Alignment against reference genomes
# Step: Genotype calling for all samples with strong evidence of Chlorobiaceae
#       DNA along the Chlorobiaceae MAG EMN001_21 using snpAD
#
# Dependent on:
#   - QUAL_refalignment_Chlorobiaceae.Snakefile
#
# Alex Huebner, 28/06/22
################################################################################

import os

import allel
import numpy as np
import pandas as pd
import pyfastx

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
alignment_stats = pd.read_csv("05-results/QUAL_dentalcalculus_Chlorobiaceae_refalignment.tsv", sep="\t")
SAMPLES = alignment_stats.query("alignedReads >= 50000")['sample'].tolist() + ['EMN001']
EMN001_CHL_MAG = "04-analysis/automatic_MAG_refinement/aDNA_samples_human/EMN001-megahit/bins/EMN001-megahit_021.fasta.gz"
CONTIGNAMES = {str(i): name for i, (name, _) in enumerate(pyfastx.Fasta(EMN001_CHL_MAG, build_index=False))}
TOTALLENGTH = sum([len(seq) for _, seq in pyfastx.Fasta(EMN001_CHL_MAG, build_index=False)])
################################################################################

wildcard_constraints:
    sample = "[A-Z]+[0-9]+",
    contig = "[0-9]+"

rule all:
    input:
        expand("04-analysis/refalignment/snpAD/{sample}.snpAD.fa", sample=SAMPLES),
        "05-results/REFG_Chlorobiaceae_genomes_snpAD.tsv"

#### Genotyping using snpAD ####################################################

rule download_snpad:
    output:
        "tmp/genotyping/snpAD-0.3.9.tar.gz"
    message: "Download the snpAD tarball"
    params:
        url = "https://bioinf.eva.mpg.de/snpAD/snpAD-0.3.9.tar.gz"
    shell:
        "wget -O {output} {params.url}"

rule extract_tarball:
    input:
        "tmp/genotyping/snpAD-0.3.9.tar.gz"
    output:
        bam2snpAD = "tmp/genotyping/snpAD-0.3.9/Bam2snpAD/Bam2snpAD",
        snpad = "tmp/genotyping/snpAD-0.3.9/snpAD/snpAD",
        snpADcall = "tmp/genotyping/snpAD-0.3.9/snpAD/snpADCall"
    message: "Extract the tarball"
    params:
        dir = "tmp/genotyping/snpAD-0.3.9"
    shell:
        """
        tar xvf {input} -C $(dirname {params.dir}) && cd {params.dir}/Bam2snpAD && make && cd ../snpAD && make
        """

rule mq_filter_bams:
    output:
        bam = temp("tmp/genotyping/{sample}.mqas_flt.bam"),
        bai = temp("tmp/genotyping/{sample}.mqas_flt.bam.bai")
    message: "Filter BAM file of sample {wildcards.sample} for MQ and AS"
    conda: "ENVS_samtools.yaml"
    params:
        bam = "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam"
    shell:
        """
        samtools view -hu -e 'mapq >= 20 || (mapq < 20 && [AS] >= -15)' \
                {params.bam} > {output.bam}
        samtools index {output.bam}
        """

rule decompress_fasta:
    output:
        temp("tmp/genotyping/EMN001-megahit_021.fasta")
    message: "Decompress FastA sequence of the Chlorobiaceae genome for use in snpAD"
    params:
        fasta = EMN001_CHL_MAG
    shell:
        "gunzip -c {params.fasta} > {output}"

rule convert_output:
    input:
        bam2snpAD = "tmp/genotyping/snpAD-0.3.9/Bam2snpAD/Bam2snpAD",
        fasta = "tmp/genotyping/EMN001-megahit_021.fasta",
        bam = "tmp/genotyping/{sample}.mqas_flt.bam",
        bai = "tmp/genotyping/{sample}.mqas_flt.bam.bai"
    output:
        temp("tmp/genotyping/{sample}.{contig}.snpad_input")
    message: "Convert BAM file into snpAD input format for sample {wildcards.sample} and contig {wildcards.contig}"
    group: "snpAD"
    params:
        contig = lambda wildcards: CONTIGNAMES[wildcards.contig]
    shell:
        """
        {input.bam2snpAD} \
            -Q 0 \
            -q 30 \
            -r {params.contig} \
            -f {input.fasta} \
            {input.bam} > {output}
        """

rule concat_snpAD_output:
    input:
        lambda wildcards: [f"tmp/genotyping/{wildcards.sample}.{contig}.snpad_input" for contig in CONTIGNAMES.keys()]
    output:
        temp("tmp/genotyping/{sample}.snpad_input")
    message: "Concatenate the snpAD output: {wildcards.sample}"
    group: "snpAD"
    shell:
        """
        cat {input} > {output}
        """

rule snpAD_estimation:
    input:
        snpad = "tmp/genotyping/snpAD-0.3.9/snpAD/snpAD",
        snpad_input = "tmp/genotyping/{sample}.snpad_input"
    output:
        priors = temp("tmp/genotyping/{sample}.priors.txt"),
        errors = temp("tmp/genotyping/{sample}.errors.txt")
    message: "Estimate the genotype likelihoods using snpAD for sample {wildcards.sample}"
    group: "snpAD"
    log: "tmp/genotyping/{sample}.snpAD.log"
    threads: 4
    shell:
        """
        {input.snpad} \
            --cpus={threads} \
            -o {output.priors} \
            -O {output.errors} \
            {input.snpad_input} > {log} 2>&1
        """

rule snpAD_modify_errors:
    input:
        "tmp/genotyping/{sample}.errors.txt"
    output:
        "tmp/genotyping/{sample}.errors.mod.txt"
    message: "Adapt the errors for haploid genotype calls: {wildcards.sample}"
    shell:
        """
        bioawk -t '$4 < 1e6{{print $1, $2, $3, $4}}' {input} > {output}
        """

rule snpAD_modify_priors:
    input:
        "tmp/genotyping/{sample}.priors.txt"
    output:
        "tmp/genotyping/{sample}.priors.mod.txt"
    message: "Adapt the priors for haploid genotype calls: {wildcards.sample}"
    run:
        with open(output[0], "wt") as outfile:
            with open(input[0], "rt") as infile:
                priors = infile.readline().rstrip().split(",")
                sum_homozygous_priors = sum([float(p) for p in priors[:4]])
                mod_priors = [f"{(float(p) / sum_homozygous_priors):.6f}"
                              if i < 4 else "1e-320"
                              for i, p in enumerate(priors)]
                outfile.write(",".join(mod_priors) + "\n")

rule snpAD_call:
    input:
        snpADcall = "tmp/genotyping/snpAD-0.3.9/snpAD/snpADCall",
        snpAD = "tmp/genotyping/{sample}.snpad_input",
        priors = "tmp/genotyping/{sample}.priors.mod.txt",
        errors = "tmp/genotyping/{sample}.errors.mod.txt"
    output:
        pipe("tmp/genotyping/{sample}.snpAD.vcf")
    message: "Call the genotypes using snpAD fixing the likelihood of a heterozygous genotype to a very small number: {wildcards.sample}"
    params:
        zeroprior = 1 / TOTALLENGTH
    shell:
        """
        {input.snpADcall} \
            -e {input.errors} \
            -p {input.priors} \
            --set_zero_priors {params.zeroprior} \
            --name={wildcards.sample} \
            {input.snpAD} > {output}
        """

rule snpAD_bgzip:
    input:
        "tmp/genotyping/{sample}.snpAD.vcf"
    output:
        "04-analysis/refalignment/snpAD/{sample}.snpAD.vcf.gz"
    message: "Compress the snpAD VCF file: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    shell:
        "bgzip --index --force {input} && mv {input}.gz* $(dirname {output})/"

rule snpAD_vcf2fa:
    input:
        "04-analysis/refalignment/snpAD/{sample}.snpAD.vcf.gz"
    output:
        "04-analysis/refalignment/snpAD/{sample}.snpAD.fa"
    message: "Convert snpAD VCF file into a FastA file for sample {wildcards.sample}"
    params:
        fasta = EMN001_CHL_MAG
    run:
        # Read VCF data
        vcf = allel.read_vcf(input[0],
                             fields=['variants/CHROM', 'variants/POS', 'variants/REF',
                                     'variants/ALT', 'variants/QUAL', 'calldata/*'])

        # Expected contigs
        contiglengths = {name: len(seq) for name, seq in pyfastx.Fasta(params.fasta, build_index=False)}
        consensus = {contig: ['N'] * contiglengths[contig]
                     for contig in contiglengths}

        # Replace Ns with genotypes when data fulfills quality requirements
        for contig in consensus:
            if contig in vcf['variants/CHROM']:
                sites_idx = np.where(vcf['variants/CHROM'] == contig)[0]
                for i in sites_idx:
                    if (vcf['variants/QUAL'][i] >= 30 and vcf['calldata/DP'][i][0] > 2 and
                            vcf['calldata/GQ'][i][0] >= 20):
                        gtfreq = np.bincount(vcf['calldata/GT'][i][0])
                        if np.max(gtfreq) == 2:  # homozygous call
                            if np.argmax(gtfreq) == 0:  # REF allele
                                gt = vcf['variants/REF'][i]
                            else:
                                gt = vcf['variants/ALT'][i][np.argmax(gtfreq) - 1]
                        else:  # heterozygous call
                            depth = np.zeros(4)
                            for j, allele in enumerate(['A', 'C', 'G', 'T']):
                                depth[j] = np.sum(vcf['calldata/' + allele][i])
                            if depth[np.argmax(depth)] / np.sum(depth) >= 0.66:
                                gt = ['A', 'C', 'G', 'T'][np.argmax(depth)]
                            else:
                                gt = 'N'
                        consensus[contig][vcf['variants/POS'][i] - 1] = gt

        # Write to FastA
        with open(output[0], "wt") as outfile:
            for contig in consensus:
                outfile.write(f">{contig}\n")
                outfile.write("".join(consensus[contig]) + "\n")

################################################################################

#### Summary statistics ########################################################

rule summary_stats:
    input:
        "04-analysis/refalignment/snpAD/{sample}.snpAD.vcf.gz"
    output:
        "04-analysis/refalignment/snpAD/{sample}.snpAD_summary"
    message: "Evaluate the fraction of the genome recovered and the number of SNPs: {wildcards.sample}"
    params:
        fasta = EMN001_CHL_MAG
    run:
        # Read VCF data
        vcf = allel.read_vcf(input[0],
                             fields=['variants/CHROM', 'variants/POS', 'variants/REF',
                                     'variants/ALT', 'variants/QUAL', 'calldata/*'])

        # Expected contigs
        contiglengths = {name: len(seq) for name, seq in pyfastx.Fasta(params.fasta, build_index=False)}
        total_length = sum([v for v in contiglengths.values()])

        # Coverage >= 3-fold
        mincov_sites = np.where((vcf['variants/QUAL'] >= 30) &
                                (vcf['calldata/DP'][:, 0] >= 3) &
                                (vcf['calldata/GQ'][:, 0] >= 20))

        # Type of genotypes
        genotypes = np.bincount(vcf['calldata/GT'][mincov_sites[0], 0, :].sum(axis=1))
        subst_sites = np.where(vcf['calldata/GT'][mincov_sites[0], 0, :].sum(axis=1) == 2)
        subst_sites = np.where((vcf['variants/QUAL'] >= 30) &
                                (vcf['calldata/DP'][:, 0] >= 3) &
                                (vcf['calldata/GQ'][:, 0] >= 20) &
                                (vcf['calldata/GT'][:, 0, :].sum(axis=1) == 2))
        substitutions = np.unique(vcf['variants/REF'][subst_sites[0]] + vcf['variants/ALT'][subst_sites[0], 0], return_counts=True)

        allele_df = pd.DataFrame.from_dict({'sample': [wildcards.sample],
                                            'sites_covered': [vcf['variants/POS'].shape[0] / total_length],
                                            'sites_genotyped': [mincov_sites[0].shape[0] / total_length],
                                            'HOMR': [genotypes[0] / genotypes.sum()],
                                            'HET': [genotypes[1] / genotypes.sum()],
                                            'HOMA': [genotypes[2] / genotypes.sum()]})
        subst_df = pd.DataFrame({'count': substitutions[1]}) \
            .transpose() \
            .reset_index(drop=True)
        subst_df.columns = substitutions[0]
        pd.concat([allele_df, subst_df], axis=1) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.5f")

rule concat_snpAD_summary:
    input:
        expand("04-analysis/refalignment/snpAD/{sample}.snpAD_summary", sample=SAMPLES[:-1])
    output:
        "05-results/REFG_Chlorobiaceae_genomes_snpAD.tsv"
    message: "Concatenate the summary stats files"
    run:
        pd.concat([pd.read_csv(fn, sep="\t")
                   for fn in input]) \
            .sort_values(['sites_covered'], ascending=False) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.5f")


################################################################################
