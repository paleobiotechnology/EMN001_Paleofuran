################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: Taxonomic profiling of the representative MAGs against the GTDBTK and
#       PhyloPhlAn3's metagenomic
#
# Dependent on:
#   - ASMB_automaticRefinement.Snakefile
#   - ASMB_dereplication.Snakefile
#
# Alex Huebner, 12/07/22
################################################################################

import os

import numpy as np
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
MAGS = pd.read_csv("05-results/ASMB_representativeMAGs.tsv", sep="\t")
REPRMAGS = MAGS['binID'].tolist()
################################################################################

#### Auxiiliary functions ######################################################

def return_fasta_fns(wildcards):
    if wildcards.mag.startswith("JAE") or wildcards.mag.startswith("VLC"):
        sampletype = "mDNA_samples_human"
        assembler = "metaspades"
    else:
        sampletype = "aDNA_samples_human"
        assembler = "megahit"
    sample = wildcards.mag.split("_")[0]
    binid = wildcards.mag.split("_")[1]
    return f"04-analysis/automatic_MAG_refinement/{sampletype}/{sample}-{assembler}/bins/{sample}-{assembler}_{binid}.fasta.gz"

################################################################################

rule all:
    input:
        "05-results/ASMB_taxonomic_profiling_MAGs.tsv"

rule decompress_fastas:
    output:
        "tmp/gtdbtk/fastas/{mag}.fa"
    message: "Decompress the FastA file: {wildcards.mag}"
    params:
        fa = lambda wildcards: return_fasta_fns(wildcards)
    threads: 1
    shell:
        """
        gunzip -c {params.fa} > {output}
        """

rule gtdbtk:
    input:
        expand("tmp/gtdbtk/fastas/{mag}.fa", mag=REPRMAGS)
    output:
        bact = "tmp/gtdbtk/gtdbtk/gtdbtk.bac120.summary.tsv",
        arch = "tmp/gtdbtk/gtdbtk/gtdbtk.ar53.summary.tsv"
    message: "Run the GTDBTK's classify workflow"
    resources:
        mem = 80,
        cores = 32
    params:
        fadir = "tmp/gtdbtk/fastas",
        outdir = "tmp/gtdbtk/gtdbtk",
        dbdir = "/mnt/archgen/users/huebner/refdbs/gtdbtk_r207_v2"
    threads: 32
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/gtdbtk/classify_wf"

rule extract_phylophlan:
    output:
        "tmp/gtdbtk/phylophlan3_repr_genomes.tsv"
    message: "Extract the PhyloPhlAn3 metagenomic results for the representative genomes"
    params:
        aDNA = "04-analysis/automatic_MAG_refinement/aDNA_samples_human/phylophlan_closestGenomes.tsv",
        mDNA = "04-analysis/automatic_MAG_refinement/mDNA_samples_human/phylophlan_closestGenomes.tsv"
    run:
        phylophlan = pd.concat([pd.read_csv(fn, sep="\t")
                                for fn in [params.aDNA, params.mDNA]])
        phylophlan['#input_bin'] = phylophlan['#input_bin'].str.replace(r'-(metaspades|megahit)', '', regex=True)
        phylophlan = phylophlan.loc[phylophlan['#input_bin'].isin(REPRMAGS)]
        phylophlan = phylophlan.iloc[:, :2]
        phylophlan.columns = ['binID', 'assignment']
        phylophlan['sgbtype'] = phylophlan['assignment'].str.split("_").str[0]
        phylophlan['sgb'] = phylophlan['assignment'].str.split(":").str[0]
        phylophlan['assignment_level'] = phylophlan['assignment'].str.split(":").str[1]
        phylophlan['taxon'] = phylophlan['assignment'].str.split(":").str[2]
        phylophlan['dist'] = phylophlan['assignment'].str.split(":").str[3].astype(float)
        phylophlan = phylophlan.drop(['assignment'], axis=1)

        taxon_levels = {'Species': -2,
                        'Genus': -3,
                        'Family': -4,
                        'Other': -5}
        phylophlan['taxon'] = ["|".join(mag.taxon.split("|")[taxon_levels[mag.assignment_level]:]) for mag in phylophlan.itertuples()]
        phylophlan['assignment_level'] = phylophlan['assignment_level'].str.lower()
            
        phylophlan.sort_values(['binID']) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.3f")

rule merge_gtdbtk_phylophlan3:
    input:
        bact = "tmp/gtdbtk/gtdbtk/gtdbtk.bac120.summary.tsv",
        arch = "tmp/gtdbtk/gtdbtk/gtdbtk.ar53.summary.tsv",
        phylophlan = "tmp/gtdbtk/phylophlan3_repr_genomes.tsv"
    output:
        "05-results/ASMB_taxonomic_profiling_MAGs.tsv"
    message: "Combine the taxonomic profiling results of GTDBTK and PhyloPhlAn3"
    run:
        # GTDBTK
        gtdbtk = pd.concat([pd.read_csv(fn, sep="\t")
                            for fn in [input.bact, input.arch]])
        notes = {'classification based on placement in class-level tree': 'class-level tree',
                 'topological placement and ANI have congruent species assignments': 'tree & ANI congruent',
                 'topological placement and ANI have incongruent species assignments': 'tree & ANI incongruent',
                 np.nan: ""}
        gtdbtk['placement'] = [notes[n] for n in gtdbtk['note'].tolist()]
        gtdbtk = gtdbtk[['user_genome', 'classification', 'red_value', 'placement']]
        gtdbtk.columns = ['binID', 'GTDB classification', 'GTDB RED value', 'GTDB placement']

        # PhyloPhlAn3
        phylophlan = pd.read_csv(input.phylophlan, sep="\t")
        phylophlan.columns = ['binID', 'SGB type', 'closest SGB',
                              'SGB assignment level', 'SGB classification', 'SGB MASH distance']

        gtdbtk.merge(phylophlan, how="left", on="binID") \
            .iloc[:, [0, 1, 7, 4, 5, 2, 8, 3, 6]] \
            .sort_values(['binID']) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.3f")
