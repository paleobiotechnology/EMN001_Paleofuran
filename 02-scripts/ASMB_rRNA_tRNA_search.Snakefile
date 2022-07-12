################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: Confirm the presence of rRNA and tRNA genes in HQ MAGs
#
# According to the MIMAG criteria, high-quality MAGs are reuqired to
# additionally have 5S, 16S, and 23S rRNA genes and at least 18 tRNA genes. I
# will follow Saheb Kashaf et al. (2021) to infer the presence of these genes.
#
# Dependent on:
#   - ASMB_automaticRefinement.Snakefile
#
# Alex Huebner, 11/07/22
################################################################################

import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
MAGS = pd.read_csv("05-results/ASMB_MAGS_metaWRAP_postfilter.tsv", sep="\t")
MAGS = MAGS.loc[MAGS['pass.MIMAG_high']]
HQMAGS = MAGS['binID'].tolist()
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
        "05-results/ASMB_rRNA_tRNA_presence.tsv"

rule rfam_download:
    output:
        "03-data/refdbs/rfam/Rfam.cm.gz"
    message: "Download the RFam database for the search for rRNA using INFERNAL"
    params:
        url = "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    shell:
        "wget -O {output} {params.url}"

rule decompress_rfam:
    input:
        "03-data/refdbs/rfam/Rfam.cm.gz"
    output:
        "03-data/refdbs/rfam/Rfam.cm"
    message: "Decompress the RFam database"
    shell:
        "gunzip -c {input} > {output}"

rule decompress_fasta:
    output:
        "tmp/infernal/{mag}.fasta"
    message: "Decompress the FastA file: {wildcards.mag}"
    params:
        fasta = lambda wildcards: return_fasta_fns(wildcards)
    shell:
        "gunzip -c {params.fasta} > {output}"

rule infernal:
    input:
        db = "03-data/refdbs/rfam/Rfam.cm.gz",
        fa = "tmp/infernal/{mag}.fasta"
    output:
        "04-analysis/ancient_metagenome_assembly/infernal_rrna/{mag}.txt"
    message: "Infer the presence of rRNA genes using INFERNAL: {wildcards.mag}"
    conda: "ENVS_infernal.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        cmsearch --cpu {threads} -Z 1000 --hmmonly --cut_ga --noali --tblout {output} -o /dev/null {input.db} {input.fa}
        """

rule trnascan_se:
    input:
        "tmp/infernal/{mag}.fasta"
    output:
        "04-analysis/ancient_metagenome_assembly/trnascan/{mag}.txt"
    message: "Search for tRNAs: {wildcards.mag}"
    conda: "ENVS_trnascanse.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        tRNAscan-SE -B -Q -o {output} --thread {threads} {input}
        """

rule summary:
    input:
        rrna = expand("04-analysis/ancient_metagenome_assembly/infernal_rrna/{mag}.txt", mag=HQMAGS),
        trna = expand("04-analysis/ancient_metagenome_assembly/trnascan/{mag}.txt", mag=HQMAGS)
    output:
        "05-results/ASMB_rRNA_tRNA_presence.tsv"
    message: "Summarise the results of scanning for rRNA and tRNA genes"
    run:
        # Infernal 
        rrna = pd.concat([pd.read_csv(fn, sep="\s+", header=None, comment="#", usecols=[0, 2, 3],
                                      names=['contig', 'query_name', 'accession'])
                          .assign(mag=os.path.basename(fn).replace(".txt", ""))
                         for fn in input.rrna])
        rrna = rrna.loc[rrna['query_name'].isin(['5S_rRNA', '23S-methyl', '6S'])]

        # tRNAscan-SE
        trna = pd.concat([pd.read_csv(fn, sep="\s+", header=None, comment="#", usecols=[0, 4],
                               names=['contig', 'tRNA'], skiprows=3) \
                          .assign(mag=os.path.basename(fn).replace(".txt", ""))
                          for fn in input.trna])

        # Combine
        rnas = trna.groupby(['mag', 'tRNA']).count() \
            .unstack() \
            .fillna(0) \
            .astype(int) \
            .reset_index()
        rnas.columns = ['mag'] + [rnas.columns[i][1] for i in range(1, rnas.shape[1])]
        rnas = rnas.drop(['Undet'], axis=1)
        rnas['total tRNAs'] = rnas.iloc[:, 1:].sum(axis=1)
        rnas = rnas.merge(rrna.groupby(['mag'])['query_name'].agg(lambda x: ", ".join(x))
                          .reset_index(), how="left", on="mag")
        rnas = rnas.rename({'mag': 'binID',
                            'query_name': 'rRNAs'}, axis=1)
        rnas['rRNAs'] = rnas['rRNAs'].fillna("")
        rnas.to_csv(output[0], sep="\t", index=False)
