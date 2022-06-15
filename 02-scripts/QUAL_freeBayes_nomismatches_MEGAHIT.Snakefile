################################################################################
# Project: Natural products from the Palaeolithic
# Part: Quality analyses
# Step: Summarise the number of mismatches between the freeBayes genotypes and
#       the MEGAHIT contigs
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 13/06/22
################################################################################

from glob import glob
import os

import allel
import numpy as np
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("04-analysis/ancient_metagenome_assembly/consensus_correction/megahit/{sample}.filter.vcf.gz")
################################################################################

rule all:
    input:
        "05-results/QUAL_freeBayes_nSubsts.tsv",
        "05-results/QUAL_freeBayes_distributionMAF.tsv",
        "05-results/QUAL_freeBayes_nContigs.tsv"

rule count_observed_mutations:
    output:
        substs = "tmp/mismatches_freebayes/{sample}.substs.txt",
        maf = "tmp/mismatches_freebayes/{sample}.mafs.txt",
        ncontigs = "tmp/mismatches_freebayes/{sample}.ncontigs.txt"
    message: "Summarise the mismatches inferred by freeBayes: {wildcards.sample}"
    resources:
        mem = 4,
        cores = 1
    params:
        vcf = "04-analysis/ancient_metagenome_assembly/consensus_correction/megahit/{sample}.filter.vcf.gz"
    threads: 1
    run:
        content = allel.read_vcf(params.vcf, fields=['variants/CHROM', 'variants/REF',
                                                     'variants/ALT', 'calldata/AD',
                                                     'calldata/DP'])
        # Type of substitutions
        substitutions = {}
        for i in range(content['variants/REF'].shape[0]):
            ref_allele = content['variants/REF'][i]
            alt_allele = content['variants/ALT'][i][0]
            if len(ref_allele) == 1:
                if f"{ref_allele}{alt_allele}" in substitutions:
                    substitutions[f"{ref_allele}{alt_allele}"] += 1
                else:
                    substitutions[f"{ref_allele}{alt_allele}"] = 1
            else:
                for j in range(len(ref_allele)):
                    if ref_allele[j] != alt_allele[j]:
                        if f"{ref_allele[j]}{alt_allele[j]}" in substitutions:
                            substitutions[f"{ref_allele[j]}{alt_allele[j]}"] += 1
                        else:
                            substitutions[f"{ref_allele[j]}{alt_allele[j]}"] = 1
        pd.DataFrame.from_dict(substitutions, orient="index") \
            .sort_index() \
            .transpose() \
            .assign(sample=wildcards.sample) \
            .to_csv(output.substs, sep="\t", index=False)

        # Minor allele frequency
        minor_alleles = content['calldata/AD'][:, :, 1] / content['calldata/DP'][:, 0][:, np.newaxis]
        ma_bincounts = np.bincount(np.digitize(minor_alleles[:,0], np.arange(0, 1, 0.05)))
        pd.DataFrame.from_dict({'nMismatches': ma_bincounts}) \
            .transpose() \
            .assign(sample=wildcards.sample) \
            .to_csv(output.maf, sep="\t", index=False)

        # Number of contigs
        unique, counts = np.unique(content['variants/CHROM'], return_counts=True)
        with open(output.ncontigs, "wt") as outfile:
            outfile.write(f"{wildcards.sample}\t{len(unique)}\n")

rule summarise_substs:
    input:
        expand("tmp/mismatches_freebayes/{sample}.substs.txt", sample=SAMPLES)
    output:
        "05-results/QUAL_freeBayes_nSubsts.tsv"
    message: "Summarise the number of substitutions observed for each sample"
    run:
        df = pd.concat([pd.read_csv(fn, sep="\t")
                        for fn in input])
        cols = df.columns.tolist()
        df = df[cols[-1:] + cols[:-1]]
        df.fillna(0) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)

rule summarise_mafs:
    input:
        expand("tmp/mismatches_freebayes/{sample}.mafs.txt", sample=SAMPLES)
    output:
        "05-results/QUAL_freeBayes_distributionMAF.tsv"
    message: "Summarise the distribution of the minor allele frequencies observed for each sample"
    run:
        df = pd.concat([pd.read_csv(fn, sep="\t")
                        for fn in input])
        cols = df.columns.tolist()
        df = df[cols[-1:] + cols[:-1]]
        df.fillna(0) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)

rule summarise_ncontigs:
    input:
        expand("tmp/mismatches_freebayes/{sample}.ncontigs.txt", sample=SAMPLES)
    output:
        "05-results/QUAL_freeBayes_nContigs.tsv"
    message: "Summarise the number of affected contigs for each sample"
    run:
        values = [tuple(open(fn, "rt").readline().rstrip().split("\t"))
                  for fn in input]

        df = pd.DataFrame(values, columns=['sample', 'nContigs'])
        df['nContigs'] = df['nContigs'].astype(int)
        df.sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)
