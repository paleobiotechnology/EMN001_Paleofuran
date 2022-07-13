################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: De-replication of the MAGs
#
# Dependent on:
#   - ASMB_automaticRefinement.Snakefile
#
# Alex Huebner, 12/07/22
################################################################################

import os

import numpy as np
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
# Assembled MAGs
MAGS = pd.read_csv("05-results/ASMB_MAGS_metaWRAP_postfilter.tsv", sep="\t")
MAGS = MAGS.loc[(MAGS['checkM.completeness'] >= 50) &
                (MAGS['checkM.contamination'] < 10) &
                (MAGS['GUNC.contamination_portion'] < 0.10) &
                (MAGS['GUNC.CSS'] < 0.45) &
                (MAGS['polyrate'] < 0.01)]
HMQMAGS = MAGS['binID'].tolist()
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


def identify_homd_oral_genomes(wildcards):
    homd = pd.read_csv(checkpoints.download_homd_info.get(**wildcards).output[0], sep=",")
    homd = homd.loc[homd['Habitat'].isin(['Oral', 'Nasal', 'Nasal | Oral'])]
    return [f"tmp/HOMD/{seqid}.fasta" for seqid in homd['SEQ_ID'].tolist()]

################################################################################

localrules: decompress_fastas

rule all:
    input:
        "05-results/ASMB_representativeMAGs.tsv",
        "05-results/ASMB_oraltaxon_classification.tsv"

#### Pick representative genomes ###############################################

rule decompress_fastas:
    output:
        "tmp/drep/mags/{mag}.fasta"
    message: "Decompress FastA file: {wildcards.mag}"
    resources:
        mem = 2
    params:
        fasta = lambda wildcards: return_fasta_fns(wildcards)
    shell:
        "gunzip -c {params.fasta} > {output}"

rule genome_info:
    output:
        "tmp/drep/genome_info.csv"
    message: "Prepare a genome info table for dRep"
    run:
        mags = MAGS[['binID', 'checkM.completeness', 'checkM.contamination']] \
            .rename({'binID': 'genome',
                    'checkM.completeness': 'completeness',
                    'checkM.contamination': 'contamination'}, axis=1)
        mags['genome'] = mags['genome'] + ".fasta"
        mags.to_csv(output[0], sep=",", index=False)

rule drep_dereplicate_species:
    input:
        fas = expand("tmp/drep/mags/{mag}.fasta", mag=HMQMAGS),
        genomeinfo = "tmp/drep/genome_info.csv"
    output:
        "04-analysis/drep/species/data_tables/Cdb.csv"
    message: "Dereplicate the MAGs for non-redundant genomes"
    conda: "ENVS_dRep.yaml"
    resources:
        mem = 100,
        cores = 32
    params:
        outdir = "04-analysis/drep/species",
        fadir = "tmp/drep/mags"
    threads: 32
    shell:
        """
        dRep dereplicate {params.outdir} \
            -p {threads} \
            --genomeInfo {input.genomeinfo} \
            -pa 0.9 -sa 0.95 \
            -nc 0.30 -cm larger -comp 50 -con 10 \
            -g {params.fadir}/*.fasta
        """

rule summarise_dereplicate_species:
    input:
        "04-analysis/drep/species/data_tables/Cdb.csv"
    output:
        "05-results/ASMB_representativeMAGs.tsv"
    message: "Summarise the dRep replication results"
    params:
        wdb = "04-analysis/drep/species/data_tables/Wdb.csv"
    run:
        # Read the representative genomes
        repr_cluster = pd.read_csv(params.wdb, sep=",")
        repr_cluster['genome'] = repr_cluster['genome'].str.replace(".fasta", "", regex=False)

        # Combine with some genome facts
        clusters = MAGS.loc[MAGS['binID'].isin(repr_cluster['genome'].tolist())] \
            [['binID', 'checkM.genome_size', 'checkM.N50_contigs',
              'checkM.completeness', 'checkM.contamination']]
        clusters.columns = ['binID', 'genome size [Mb]', 'N50', 'completeness [%]', 'contamination [%]']
        clusters['genome size [Mb]'] /= 1e6
        clusters = clusters.merge(repr_cluster.rename({'genome': 'binID'}, axis=1),
                                  how="left", on="binID")

        # Add non-representative MAGs
        nonrepr_cluster = pd.read_csv(input[0], sep=",",
                                      usecols=['genome', 'secondary_cluster'])
        nonrepr_cluster['genome'] = nonrepr_cluster['genome'].str.replace(".fasta", "", regex=False)
        nonrepr_cluster['repr'] = nonrepr_cluster['genome'].isin(clusters['binID'].tolist())
        nonrepr_cluster = nonrepr_cluster.loc[~nonrepr_cluster['repr']]
        nonrepr_cluster = nonrepr_cluster.groupby(['secondary_cluster'])['genome'] \
            .agg([lambda x: len(x) + 1, lambda x: ", ".join(x)]) \
            .reset_index()
        nonrepr_cluster.columns = ['cluster', 'cluster size', 'members of cluster'] 

        clusters = clusters.merge(nonrepr_cluster, how="left", on="cluster")
        clusters['cluster size'] = clusters['cluster size'].fillna(1).astype(int)
        clusters['members of cluster'] = clusters['members of cluster'].fillna("")
        clusters['pc'] = clusters['cluster'].str.split("_").str[0].astype(int)
        clusters['sc'] = clusters['cluster'].str.split("_").str[1].astype(int)

        # Save
        clusters.sort_values(['pc', 'sc']) \
            .iloc[:, [0, 5, 1, 2, 3, 4, 6, 7, 8]] \
            .to_csv(output[0], sep="\t", index=False, float_format="%.2f")
            
################################################################################

#### Identify oral taxa ########################################################

checkpoint download_homd_info:
    output:
        "tmp/HOMD/SEQID_info.csv"
    message: "Download HOMD database genome overview"
    params:
        url = "https://homd.org/ftp/genomes/NCBI/V9.15a/SEQID_info.csv"
    shell:
        "wget -O {output} {params.url}"

rule download_homd_fasta:
    output:
        "tmp/HOMD/{seqid}.fasta"
    message: "Download the FastA file from the HOMD FTP server: {wildcards.seqid}"
    params:
        url = "https://homd.org/ftp/genomes/NCBI/V9.15a/fna/{seqid}.fna"
    shell:
        "wget -O {output} {params.url}"

rule link_repr_genomes:
    input:
        "04-analysis/drep/species/data_tables/Cdb.csv"
    output:
        touch("tmp/HOMD/link_repr_genomes.done")
    message: "Link all representative genomes into the folder of the HOMD genomes"
    params:
        wdb = "04-analysis/drep/species/data_tables/Wdb.csv"
    run:
        # Read the representative genomes
        repr_cluster = pd.read_csv(params.wdb, sep=",")
        for genome in repr_cluster['genome'].tolist():
            os.symlink(f"{os.getcwd()}/tmp/drep/mags/{genome}", f"tmp/HOMD/{genome}")

rule drep_compare:
    input:
        homd = identify_homd_oral_genomes,
        repr_clusters = "tmp/HOMD/link_repr_genomes.done"
    output:
        "04-analysis/drep/homd/data_tables/Cdb.csv"
    message: "Dereplicate the MAGs with the genomes of the HOMD to identify oral species"
    conda: "ENVS_dRep.yaml"
    resources:
        mem = 100,
        cores = 32
    params:
        outdir = "04-analysis/drep/homd",
        fadir = "tmp/HOMD"
    threads: 32
    shell:
        """
        dRep compare {params.outdir} \
            -p {threads} \
            -pa 0.8 -sa 0.95 \
            -nc 0.30 -cm larger \
            -g {params.fadir}/*.fasta
        """

rule annotate_homd:
    input:
        species = "04-analysis/drep/species/data_tables/Cdb.csv",
        homd = "04-analysis/drep/homd/data_tables/Cdb.csv"
    output:
        "05-results/ASMB_oraltaxon_classification.tsv"
    message: "Assign the representative MAGs to oral taxa based on the dRep results"
    params:
        ani = "04-analysis/drep/homd/data_tables/Ndb.csv",
        repr_cluster = "04-analysis/drep/species/data_tables/Wdb.csv"
    run:
        # Read the representative genomes
        repr_cluster = pd.read_csv(params.repr_cluster, sep=",")
        repr_cluster['genome'] = repr_cluster['genome'].str.replace(".fasta", "", regex=False)

        # Load cluster assignment
        drep = pd.read_csv(input.homd, sep=",", usecols=[0, 1])
        drep['genome'] = drep['genome'].str.replace(".fasta", "", regex=False)
        drep = drep.set_index(['genome'])
        drep['pc'] = drep['secondary_cluster'].str.split("_").str[0].astype(int)
        drep['HOMD'] = drep.index.str.startswith("SEQ")
        # Load ANI results
        ani = pd.read_csv(params.ani, sep=",")
        ani['reference'] = ani['reference'].str.replace(".fasta", "", regex=False)
        ani['querry'] = ani['querry'].str.replace(".fasta", "", regex=False)

        oraltaxon = []
        for mag in repr_cluster['genome'].tolist():
            pc = drep.loc[drep['pc'] == drep.at[mag, 'pc']]
            sc = drep.loc[drep['secondary_cluster'] == drep.at[mag, 'secondary_cluster']]
            if pc['HOMD'].sum() > 0:
                if sc['HOMD'].sum() > 0:
                    refs = sc.loc[sc.index.str.startswith("SEQ")].index.tolist()
                    pwani = ani.loc[(ani['reference'] == mag) & (ani['querry'].isin(refs)), 'ani'].max()
                    ref = ani.loc[(ani['reference'] == mag) &
                                  (ani['querry'].isin(refs)) &
                                  (ani['ani'] == pwani), 'querry'].values[0]
                    oraltaxon.append((mag, "secondary", pwani, ref))
                else:
                    refs = pc.loc[pc.index.str.startswith("SEQ")].index.tolist()
                    pwani = ani.loc[(ani['reference'] == mag) & (ani['querry'].isin(refs)), 'ani'].max()
                    ref = ani.loc[(ani['reference'] == mag) &
                                  (ani['querry'].isin(refs)) &
                                  (ani['ani'] == pwani), 'querry'].values[0]
                    oraltaxon.append((mag, "primary", pwani, ref))
            else:
                oraltaxon.append((mag, "none", np.nan, ""))
            
        pd.DataFrame(oraltaxon, columns=['binID', 'oral taxon', 'ANI', 'closest HOMD genome']) \
            .sort_values(['binID']) \
            .to_csv(output[0], sep="\t", index=False)

################################################################################
