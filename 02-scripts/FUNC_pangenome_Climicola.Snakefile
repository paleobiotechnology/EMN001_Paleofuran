################################################################################
# Project: Natural products from the Palaeolithic
# Part: Functional analysis
# Step: Pan-genome analysis between the Chlorobium MAGs and Chlorobium limicola
#
# Dependent on:
#   - ASMB_automaticRefinement.Snakefile
#
# Alex Huebner, 09/09/22
################################################################################

from glob import glob
import os
from urllib.error import HTTPError

from Bio.KEGG import REST
import numpy as np
import pandas as pd
import pyfastx

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
MAGS, = glob_wildcards("05-results/genomes/chlorobiaceae/{sample}.fa.gz")
MAGS = [mag for mag in MAGS if not mag.endswith("_RA")]
CLIMICOLA = {"GCF_001509575.1": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/509/575/GCF_001509575.1_ASM150957v1/GCF_001509575.1_ASM150957v1_genomic.fna.gz",
             "GCA_013335335.1": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/335/335/GCA_013335335.1_ASM1333533v1/GCA_013335335.1_ASM1333533v1_genomic.fna.gz",
             "GCF_000020465.1": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/465/GCF_000020465.1_ASM2046v1/GCF_000020465.1_ASM2046v1_genomic.fna.gz",
             "GCA_013335765.1": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/335/765/GCA_013335765.1_ASM1333576v1/GCA_013335765.1_ASM1333576v1_genomic.fna.gz"}
################################################################################

#### Auxilliary functions ######################################################

def return_path_gff(wildcards):
    sampleid, binid = wildcards.mag.split("_")
    return f"{os.getcwd()}/04-analysis/automatic_MAG_refinement/aDNA_samples_human/{sampleid}-megahit/bins/{sampleid}-megahit_{binid}.gff"

################################################################################

wildcard_constraints:
    genome = "GC[AF]_[0-9]+\.[0-9]",
    mag = "[A-Z]+[0-9]+_[0-9]+"

rule all:
    input:
        "05-results/FUNC_roary_pangenome.tsv",
        "05-results/FUNC_KEGG_specificKO.tsv"

#### Annotate C. limicola genomes with Prokka ##################################

rule decompress_fasta:
    output:
        temp("tmp/prokka_climicola/{genome}.fasta")
    message: "Decompress the FastA file: {wildcards.genome}"
    params:
        url = lambda wildcards: CLIMICOLA[wildcards.genome]
    shell:
        "wget -O - {params.url} | gunzip > {output}"

rule prokka:
    input:
        "tmp/prokka_climicola/{genome}.fasta"
    output:
        "tmp/prokka_climicola/prokka_{genome}/{genome}.gff"
    message: "Run Prokka on contigs: {wildcards.genome}"
    conda: "ENVS_prokka.yaml"
    resources:
        mem = 16,
        cores = 8
    params:
        tmpdir = "tmp/prokka_climicola/prokka_{genome}"
    threads: 8
    shell:
        """
        mkdir -p {params.tmpdir}
        prokka --outdir {params.tmpdir} \
               --prefix {wildcards.genome} \
               --force \
               --compliant \
               --cpus {threads} \
               --debug \
               {input}
        """

rule copy_genomes_gff:
    input:
        "tmp/prokka_climicola/prokka_{genome}/{genome}.gff"
    output:
        "04-analysis/roary/gffs/{genome}.gff"
    message: "Copy the GFF file: {wildcards.genome}"
    shell:
        "cp {input} {output}"

rule mag_decompress_gff:
    output:
        "04-analysis/roary/gffs/{mag}.gff"
    message: "Decompress GFF file: {wildcards.mag}"
    params:
        gff = lambda wildcards: return_path_gff(wildcards)
    shell:
        "ln -s {params.gff} {output}"

################################################################################

#### Pan-genome analysis with roary ############################################

checkpoint roary:
    input:
        genomes = expand("04-analysis/roary/gffs/{genome}.gff", genome=CLIMICOLA.keys()),
        mags = expand("04-analysis/roary/gffs/{mag}.gff", mag=MAGS)
    output:
        touch("04-analysis/roary/roary.done")
    message: "Perform the pan-genome analysis with Roary"
    conda: "ENVS_roary.yaml"
    params:
        dir = "04-analysis/roary",
        gffdir = "04-analysis/roary/gffs"
    threads: 16
    shell:
        """
        roary -p 8 \
            -f {params.dir}/ \
            -e -n -v \
            {params.gffdir}/*.gff
        """

################################################################################

#### Profile proteins with eggNOG-mapper #######################################

rule identify_climicola_specific_genes:
    input:
        "04-analysis/roary/roary.done"
    output:
        faa = "04-analysis/roary/climicola_specific_genes.faa",
        tsv = "04-analysis/roary/gene_presence_absence_annot.tsv"
    message: "Prepare FastA file with the aminoacid sequences of the genes specific to the C. limicola genomes"
    params:
        prokka_dir = "tmp/prokka_climicola/prokka",
        roary_dir = "04-analysis/roary/"
    run:
        roary_res_dir = glob(f"{params.roary_dir}/_*")[0]
        genes = pd.read_csv(f"{roary_res_dir}/gene_presence_absence.csv", sep=",")
        # Determine prevalence
        prev_mags = (~genes.iloc[:, [14, 19, 20, 21, 22, 23]].isnull()).sum(axis=1) / 6
        prev_genomes = (~genes.iloc[:, 15:19].isnull()).sum(axis=1) / 4
        # Identify the C. limicola specific proteins
        genome_specific = np.where((prev_genomes >= 0.5) & (prev_mags < 0.4))[0]
        genome_specific_genes = set(genes.iloc[genome_specific, 15:19]
                                    .apply(lambda r: r.loc[~r.isnull()].values[0], axis=1))
        with open(output.faa, "wt") as outfile:
            for genome in CLIMICOLA:
                for name, seq in pyfastx.Fasta(f"{params.prokka_dir}_{genome}/{genome}.faa", build_index=False):
                    if name in genome_specific_genes:
                        outfile.write(f">{name}\n{seq}\n")

        # Identify the Chlorobium MAG specific proteins
        mag_specific = np.where((prev_genomes < 0.5) & (prev_mags >= 0.6))[0]
        genes['specificity'] = "none"
        genes.iloc[genome_specific, -1] = "C. limicola"
        genes.iloc[mag_specific, -1] = "Chlorobium MAG"
        genes.to_csv(output.tsv, sep="\t", index=False)

rule download_eggnog_db:
    output:
        touch("03-data/refdbs/eggNOG/download.done")
    message: "Download the eggNOG database"
    conda: "ENVS_eggnogmapper.yaml"
    params:
        dir = "03-data/refdbs/eggNOG"
    shell:
        "download_eggnog_data.py --data_dir {params.dir}"

rule eggnog_mapper:
    input:
        faa = "04-analysis/roary/climicola_specific_genes.faa",
        db = "03-data/refdbs/eggNOG/download.done"
    output:
        "04-analysis/roary/eggnog/climicola_specific_genes.emapper.annotations"
    message: "Screen the C. limicola specific genes with eggnog-mapper"
    conda: "ENVS_eggnogmapper.yaml"
    params:
        db_dir = "03-data/refdbs/eggNOG",
        outdir = "04-analysis/roary/eggnog",
        tmp_dir = "/tmp"
    threads: 32
    shell:
        """
        emapper.py --cpu {threads} \
            -i {input.faa} \
            --data_dir {params.db_dir} \
            -m diamond \
            --output climicola_specific_genes \
            --output_dir {params.outdir} \
            --temp_dir {params.tmp_dir}
        """

rule annotate_roary_results:
    input:
        tsv = "04-analysis/roary/gene_presence_absence_annot.tsv",
        annot = "04-analysis/roary/eggnog/climicola_specific_genes.emapper.annotations"
    output:
        "05-results/FUNC_roary_pangenome.tsv"
    message: "Annotate the roary results with the eggNOG annotations"
    params:
        xlsx = "06-figures_tables/Dataset_S6.xlsx"
    run:
        # Roary results
        roary = pd.read_csv(input.tsv, sep="\t")
        # Annotations
        eggnog_genomes = pd.read_csv("04-analysis/roary/eggnog/climicola_specific_genes.emapper.annotations",
                                     sep="\t", skiprows=4, skipfooter=3, engine="python",
                                     usecols=[0, 6, 7, 8, 10, 11, 13, 14, 15, 18, 20])
        eggnog_mags = pd.concat([pd.read_excel(params.xlsx, sheet_name=i,
                                               skiprows=2, skipfooter=3)
                                 for i in [0, 1, 3, 5, 6, 7]])
        eggnog_mags = eggnog_mags[['#query', 'COG_category', 'Description', 'Preferred_name', 'EC',
                     'KEGG_ko', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'CAZy', 'PFAMs']]

        # Select representative gene
        roary['repr_protein'] = np.nan
        roary.iloc[np.where(roary['specificity'] == "C. limicola")[0], -1] = \
            roary.iloc[np.where(roary['specificity'] == "C. limicola")[0], 15:19] \
            .apply(lambda r: r.loc[~r.isnull()].values[0], axis=1)
        roary.iloc[np.where(roary['specificity'] == "Chlorobium MAG")[0], -1] = \
            roary.iloc[np.where(roary['specificity'] == "Chlorobium MAG")[0], [14, 19, 20, 21, 22, 23]] \
            .apply(lambda r: r.loc[~r.isnull()].values[0], axis=1)

        roary['repr_protein'] = roary['repr_protein'].str.replace("DOAKPEPN", "EMN001_021")
        roary['repr_protein'] = roary['repr_protein'].str.replace("JHDCJAEM", "GOY005_001")
        for i in range(14, 24):
            roary.iloc[:, i] = roary.iloc[:, i].str.replace("^[A-Z]+_", "", regex=True)

        roary.merge(pd.concat([eggnog_genomes, eggnog_mags]),
                    how="left", left_on="repr_protein", right_on="#query") \
            .drop(['#query', 'repr_protein'], axis=1) \
            .iloc[:, list(range(14)) + list(range(15, 19)) + [14] + list(range(19, 35))] \
            .to_csv(output[0], sep="\t", index=False)

rule kegg_pathway_analysis:
    input:
        "05-results/FUNC_roary_pangenome.tsv"
    output:
        "05-results/FUNC_KEGG_specificKO.tsv"
    message: "Infer the associated pathways for specific KOs using KEGG's REST API"
    run:
        # Read data
        roary = pd.read_csv(input[0], sep="\t")

        # Count the occurrences of a KO per group
        spec_genes_ko = roary.query('specificity != "none"') \
            .query('KEGG_ko != "-"') \
            .query('~KEGG_ko.isnull()')[['specificity', 'KEGG_ko']]
        spec_genes_ko['KEGG_ko'] = spec_genes_ko['KEGG_ko'].str.split(",")
        spec_genes_ko = spec_genes_ko.explode('KEGG_ko')
        ko_counts = spec_genes_ko.groupby(['specificity'])['KEGG_ko'].value_counts()

        # Determine the KOs specific for a group
        clim_kos = set(ko_counts.loc['C. limicola'].index.tolist())
        cmag_kos = set(ko_counts.loc['Chlorobium MAG'].index.tolist())
        kos_only_clim = clim_kos.difference(cmag_kos)
        kos_only_cmag = cmag_kos.difference(clim_kos)
        specific_kos = kos_only_clim.union(kos_only_cmag)

        kegg_ko_info = {}
        for n, k in enumerate(specific_kos):
            print(f"{n}:{k}")
            try:
                results = REST.kegg_get(k[3:]).read().split("\n")
                brite_start = [i for i, line in enumerate(results) if line.startswith("BRITE")][0]
                brite_end = [i for i, line in enumerate(results) if line.startswith("GENES")][0]
                kegg_ko_info[k] = "\n".join(results[brite_start:brite_end])
            except HTTPError:
                pass

        ko_tbl = [tuple(e.split("\n")[1:5]) for k, e in kegg_ko_info.items()]
        ko_tbl_df = pd.DataFrame(ko_tbl, columns = ['pathway_1stl', 'pathway_2ndl', 'pathway_3rdl', 'ko'])
        ko_tbl_df['pathway_1stl'] = ko_tbl_df['pathway_1stl'].str.replace(r'^PATHWAY +', '', regex=True)
        ko_tbl_df['pathway_2ndl'] = ko_tbl_df['pathway_2ndl'].str.replace(r'^ +', '', regex=True)
        ko_tbl_df['pathway_3rdl'] = ko_tbl_df['pathway_3rdl'].str.replace(r'^ +', '', regex=True)
        ko_tbl_df['ko'] = ko_tbl_df['ko'].str.replace(r'^ +', '', regex=True)
        ko_tbl_df['ko_desc'] = ko_tbl_df['ko'].str.split("  ").str[1]
        ko_tbl_df['ko'] = ko_tbl_df['ko'].str.split("  ").str[0]
        ko_tbl_df['specificity'] = ["C. limicola" if f"ko:{k}" in kos_only_clim else  "Chlorobium MAG" for k in ko_tbl_df['ko']]

        ko_tbl_df[['ko', 'pathway_1stl', 'pathway_2ndl', 'pathway_3rdl', 'ko_desc', 'specificity']] \
            .to_csv(output[0], sep="\t", index=False)


################################################################################
