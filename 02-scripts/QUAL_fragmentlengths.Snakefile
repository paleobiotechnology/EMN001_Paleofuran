################################################################################
# Project: Natural products from the Palaeolithic
# Part: Quality analyses
# Step: Infer the fragment lengths of DNA molecules for the dental calculus
#       samples
#
# Dependent on:
#   - PREP_preprocessing_dentalcalculus_sequencing_data.Snakefile
#
# Alex Huebner, 13/06/22
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SEQDATA = pd.read_csv("01-resources/overview_sequencingdata.tsv", sep="\t") \
    .query("sequencingSetup == 'PE75'")
LIBS = SEQDATA['libraryId']
SEQRUNS = {}
for s in LIBS:
    SEQRUNS[s] = {}
    for fn in glob(f"04-analysis/eager/results_pe/adapterremoval/output/{s}*.pe.settings"):
        SEQRUNS[s][os.path.basename(fn).replace(".settings", "")] = (fn.replace(".settings", ""))
################################################################################

rule all:
    input:
        "05-results/QUAL_fragmentlength_distribution.tsv"

rule collapse:
    output:
        pe1 = temp("tmp/fastp_collapse/{sample}-{seqrun}_1.fastq.gz"),
        pe2 = temp("tmp/fastp_collapse/{sample}-{seqrun}_2.fastq.gz"),
        coll = temp("tmp/fastp_collapse/{sample}-{seqrun}_merged.fastq.gz")
    message: "Merge overlapping reads using fastp: {wildcards.seqrun} for sample {wildcards.sample}"
    conda: "ENVS_fastp.yaml"
    resources:
        mem = 8
    params:
        pe1 = lambda wildcards: f"{SEQRUNS[wildcards.sample][wildcards.seqrun]}.pair1.truncated.gz",
        pe2 = lambda wildcards: f"{SEQRUNS[wildcards.sample][wildcards.seqrun]}.pair2.truncated.gz",
    shell:
        """
        fastp --in1 {params.pe1} \
              --in2 {params.pe2} \
              --merge \
              --out1 {output.pe1} \
              --out2 {output.pe2} \
              --merged_out {output.coll} \
              --overlap_len_require 11 \
              -A -G -Q -j /dev/null -h /dev/null
        """

rule count_collapsed:
    input:
        "tmp/fastp_collapse/{sample}-{seqrun}_merged.fastq.gz"
    output:
        "04-analysis/fragmentlength/{sample}-{seqrun}_collapsed.hist"
    message: "Infer the distribution of the insert size for collapsed reads: {wildcards.seqrun} for sample {wildcards.sample}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 4
    shell:
        """
        bioawk -c fastx '{{print length($seq)}}' {input} | \
            sort -k1,1n | \
            uniq -c > {output}
        """

rule count_noncollapsed:
    input:
        pe1 = "tmp/fastp_collapse/{sample}-{seqrun}_1.fastq.gz",
        pe2 = "tmp/fastp_collapse/{sample}-{seqrun}_2.fastq.gz"
    output:
        "04-analysis/fragmentlength/{sample}-{seqrun}_noncollapsed.n"
    message: "Count the number of DNA molecules that could not be collapsed: {wildcards.seqrun} for sample {wildcards.sample}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    shell:
        """
        bioawk -c fastx 'END{{print NR}}' {input} > {output}
        """

rule summarise:
    input:
        hist = [f"04-analysis/fragmentlength/{sample}-{seqrun}_collapsed.hist" for sample in SEQRUNS for seqrun in SEQRUNS[sample]],
        noncollapsed = [f"04-analysis/fragmentlength/{sample}-{seqrun}_noncollapsed.n" for sample in SEQRUNS for seqrun in SEQRUNS[sample]]
    output:
        "05-results/QUAL_fragmentlength_distribution.tsv"
    message: "Summarise the fragment length distributions"
    run:
        # Summarise the histograms
        histogram = pd.concat([pd.read_csv(fn, sep="\s+",
                                           header=None, names=['number_molecules', 'fragmentlength']) \
                                  .query("fragmentlength >= 30") \
                                  .assign(lib=os.path.basename(fn).split("-")[0])
                               for fn in input.hist])
        histogram['sample'] = histogram['lib'].str.split(".").str[0]
        histogram = histogram.groupby(['sample', 'fragmentlength'])['number_molecules'] \
            .agg('sum') \
            .unstack() \
            .reset_index()

        # Summarise non-collapsed molecules
        noncollapsed = pd.DataFrame([(os.path.basename(fn).split("-")[0], int(open(fn, "rt").readline().rstrip()))
                                     for fn in input.noncollapsed],
                                     columns=['lib', '> 140 bp'])
        noncollapsed['sample'] = noncollapsed['lib'].str.split(".").str[0]
        noncollapsed = noncollapsed.groupby(['sample'])['> 140 bp'].agg('sum') \
            .reset_index()

        # Calculate relative values
        histogram = histogram.merge(noncollapsed, how="left", on="sample")
        histogram['total'] = histogram.iloc[:, 1:].sum(axis=1)
        histogram.iloc[:, 1:113] = histogram.iloc[:, 1:113].apply(lambda c: c / histogram['total'].values, axis=0)
        histogram = histogram.iloc[:, [0, 113] + list(range(1, 113))]
        histogram.columns = ['sample', 'no. of DNA molecules'] + [f"{i}bp" for i in range(30, 141)] + ["> 140bp"]
        histogram.to_csv(output[0], sep="\t", index=False, float_format="%.4f")
