################################################################################
# Project: Natural products from the Palaeolithic
# Part: Phylogenetic analyses
# Step: Protein-alignment tree of the family Anaerolineaceae using PhyloPhlAn3
#
# Dependent on:
#   - ASMB_automaticRefinement.Snakefile
#
# Alex Huebner, 26/07/22
################################################################################

from glob import glob
import os
import shutil

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("05-results/genomes/flexilinea/{sample}.fa.gz")
################################################################################

#### Auxilliary functions ######################################################

def list_flexilinea_accids(wildcards):
    filepaths = [line.rstrip()
                 for line in open(checkpoints.extract_genus_Flexilinea_urls.get(**wildcards).output[0], "rt")]
    accid_list = ["_".join(os.path.basename(fn).split("_")[:2]) for fn in filepaths]
    return [f"03-data/refgenomes/Flexilinea/{accid}.fna.gz" for accid in accid_list]


def determine_assembly_url(wildcards):
    filepaths = {"_".join(os.path.basename(line).split("_")[:2]): line.rstrip().replace("_assembly_report.txt", "_genomic.fna.gz")
                 for line in open(checkpoints.extract_genus_Flexilinea_urls.get(**wildcards).output[0], "rt")}
    return filepaths[wildcards.accid]


def return_fasta_fn(wildcards):
    mag_fasta_fns = []
    with open(checkpoints.fastani_samplelist.get(**wildcards).output[0], "rt") as infile:
        for line in infile:
            mag_fasta_fns.append(line.rstrip())
    return mag_fasta_fns


################################################################################

rule all:
    input:
        "05-results/PHYL_Flexilinea_proteintree_taxa.tsv",
        "05-results/PHYL_Flexilinea_proteintree_RAxML.tre",
        "05-results/PHYL_Flexilinea_proteintree_RAxML_bootstrap.tre"

checkpoint extract_genus_Flexilinea_urls:
    output:
        "04-analysis/phylogenetics/Flexilinea/ncbi_assembly_flexilinea_urls.txt"
    message: "Downlaod the list of genomes of the genus Flexilinea from NCBI Assembly"
    conda: "ENVS_entrez.yaml"
    shell:
        """
        esearch -db assembly -query '"Anaerolineaceae"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])' | \
        esummary | xtract -pattern DocumentSummary -element FtpPath_Assembly_rpt > {output}
        """

rule extract_genus_Flexilinea_taxids:
    output:
        "04-analysis/phylogenetics/Flexilinea/ncbi_assembly_flexilinea_taxids.txt"
    message: "Extract the information on the AssemblyAccession, Taxid, and SpeciesName from NCBI Assembly"
    conda: "ENVS_entrez.yaml"
    shell:
        """
        esearch -db assembly -query '"Anaerolineaceae"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])' | \
        esummary | xtract -pattern DocumentSummary -element AssemblyAccession,SpeciesName,Taxid > {output}
        """

rule download_fastas_Flexilinea_genomes:
    output:
        "03-data/refgenomes/Flexilinea/{accid}.fna.gz",
    message: "Download the genome from NCBI: {wildcards.accid}"
    params:
        url = lambda wildcards: determine_assembly_url(wildcards)
    shell:
        """
        wget -O {output} {params.url} || touch {output}
        """

rule download_fasta_outgroup:
    output:
        "03-data/refgenomes/Flexilinea/Cand_Desulfolinea_nitratireducens.fna.gz"
    params:
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/382/535/GCA_014382535.1_ASM1438253v1/GCA_014382535.1_ASM1438253v1_genomic.fna.gz"
    shell:
        """
        wget -O {output} {params.url}
        """

rule copy_refgenomes:
    input:
        refgenomes = list_flexilinea_accids,
        outgroup = "03-data/refgenomes/Flexilinea/Cand_Desulfolinea_nitratireducens.fna.gz"
    output:
        touch("tmp/Flexilinea_tree/refgenomes.done")
    message: "Copy the reference genomes of the genus Flexilinea and the outgroup"
    params:
        dir = "03-data/refgenomes/Flexilinea",
        outdir = "tmp/Flexilinea_tree/refgenomes"
    shell:
        """
        mkdir -p {params.outdir}
        for fn in $(find {params.dir} -name "*.fna.gz"); do
            gunzip -c ${{fn}} > {params.outdir}/$(basename ${{fn}} .fna.gz).fna
        done
        """

rule drep_dereplicate:
    input:
        "tmp/Flexilinea_tree/refgenomes.done"
    output:
        "tmp/Flexilinea_tree/drep/data_tables/Wdb.csv"
    message: "Dereplicate genomes of family Chlorobiaceae using dRep"
    conda: "ENVS_dRep.yaml"
    params:
        outdir = "tmp/Flexilinea_tree/drep",
        indir = "tmp/Flexilinea_tree/refgenomes"
    threads: 16
    shell:
        """
        dRep dereplicate {params.outdir} \
            -p {threads} \
            --ignoreGenomeQuality \
            -g {params.indir}/*.fna
        """

rule mag_genomes:
    output:
        touch("tmp/Flexilinea_tree/mags.done")
    message: "Copy the Flexilinea MAGs"
    params:
        dir = "../EMN001_Paleofuran/05-results/genomes/flexilinea",
        outdir = "tmp/Flexilinea_tree/mags"
    shell:
        """
        mkdir -p {params.outdir}
        for fn in {params.dir}/*.fa.gz; do
            gunzip -c ${{fn}} > {params.outdir}/$(basename ${{fn}} .gz)
        done
        """

rule link_genomes:
    input:
        mags = "tmp/Flexilinea_tree/mags.done",
        refgenomes = "tmp/Flexilinea_tree/drep/data_tables/Wdb.csv"
    output:
        touch("tmp/Flexilinea_tree/link_genomes.done")
    message: "Link all genomes into a single folder"
    params:
        mags = "tmp/Flexilinea_tree/mags",
        refgenomes = "tmp/Flexilinea_tree/drep/dereplicated_genomes",
        outdir = "tmp/Flexilinea_tree/genomes"
    shell:
        """
        mkdir -p {params.outdir}
        for fn in {params.mags}/*.fa; do
            ln -s ${{PWD}}/${{fn}} {params.outdir}/
        done
        for fn in {params.refgenomes}/*.fna; do
            ln -s ${{PWD}}/${{fn}} {params.outdir}/$(basename ${{fn}} .fna).fa
        done
        """

rule phylophlan3_write_config:
    output:
        "04-analysis/phylogenetics/Flexilinea/config.cfg"
    message: "Write the config file for PhyloPhlAn3"
    conda: "ENVS_PhyloPhlAn3.yaml"
    shell:
        """
        phylophlan_write_config_file \
            -d a \
            -o {output} \
            --db_aa diamond \
            --map_dna diamond \
            --map_aa diamond \
            --msa mafft \
            --trim trimal \
            --tree1 fasttree \
            --tree2 raxml
        """

rule setup_flexilinea_coregenes:
    output:
        "tmp/Flexilinea_tree/core_genes/s__Flexilinea_flocculi/s__Flexilinea_flocculi.faa"
    message: "Download the core genes of Flexilinea flocculi and generate a database"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        dir = "tmp/Flexilinea_tree/core_genes"
    shell:
        """
        mkdir -p {params.dir} && \
        phylophlan_setup_database \
            -g s__Flexilinea_flocculi --database_update \
            -o {params.dir} \
            --verbose
        """

rule phylophlan3_tree:
    input:
        decompress = "tmp/Flexilinea_tree/mags.done",
        database = "tmp/Flexilinea_tree/core_genes/s__Flexilinea_flocculi/s__Flexilinea_flocculi.faa",
        config = "04-analysis/phylogenetics/Flexilinea/config.cfg"
    output:
        "04-analysis/phylogenetics/Flexilinea/RAxML_bestTree.genomes_refined.tre"
    message: "Run PhyloPhlAn3 to generate a tree of the genus Flexilinea"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        genomes = "tmp/Flexilinea_tree/genomes",
        db_folder = "tmp/Flexilinea_tree/core_genes",
        outdir = "04-analysis/phylogenetics/Flexilinea"
    log: "04-analysis/phylogenetics/Flexilinea/phylophlan.log"
    threads: 36
    shell:
        """
        mkdir -p {params.genomes}/phylophlan
        phylophlan \
            -i {params.genomes} \
            -d s__Flexilinea_flocculi \
            --databases_folder {params.db_folder} \
            --diversity medium \
            --accurate \
            --genome_extension .fa \
            -f {input.config} \
            -o Flexilinea \
            --output_folder {params.outdir} \
            --nproc {threads} \
            --verbose 2>&1 | tee {log}
        for fn in genomes_concatenated.aln genomes.tre genomes_resolved.tre \
                  RAxML_bestTree.genomes_refined.tre RAxML_info.genomes_refined.tre RAxML_result.genomes_refined.tre; do
            mv {params.outdir}/Flexilinea/${{fn}} {params.outdir}/
        done
        """

rule summarise_results:
    input:
        taxids = "04-analysis/phylogenetics/Flexilinea/ncbi_assembly_flexilinea_taxids.txt",
        drep = "tmp/Flexilinea_tree/drep/data_tables/Wdb.csv",
        tree = "04-analysis/phylogenetics/Flexilinea/RAxML_bestTree.genomes_refined.tre"
    output:
        taxa = "05-results/PHYL_Flexilinea_proteintree_taxa.tsv",
        tree = "05-results/PHYL_Flexilinea_proteintree_RAxML.tre"
    message: "Summarise the dereplication of the taxa and provide the RAxML tree"
    run:
        # List of downloaded genomes
        taxids = pd.concat([pd.read_csv(input.taxids, sep="\t", header=None,
                                        names=['accession Id', 'name', 'NCBI taxonomy Id']),
                            pd.DataFrame.from_dict({'accession Id': "GCA_014382535.1",
                                                    'name': "Candidatus Desulfolinea nitratireducens",
                                                    'NCBI taxonomy Id': 2841698}, orient="index") \
                            .transpose()])

        # Dereplication results
        drep_report = pd.read_csv(input.drep, sep=",")
        drep_report['genome'] = drep_report['genome'].str.replace(".fna", "")
        taxids = taxids.loc[(taxids['accession Id'].isin(drep_report['genome'].tolist())) |
                            (taxids['accession Id'] == "GCA_014382535.1")]

        taxids.to_csv(output.taxa, sep="\t", index=False)
        shutil.copy2(input.tree, output.tree)

################################################################################

#### Calculate pairwise ANI ####################################################

checkpoint fastani_samplelist:
    input:
        "05-results/PHYL_Flexilinea_proteintree_taxa.tsv"
    output:
        "04-analysis/phylogenetics/Flexilinea/flexilinea_fasta.txt"
    message: "Prepare the input list for fastANI from all sequences belonging to the genus Flexilinea"
    params:
        amag = "05-results/genomes/flexilinea",
        refgenomes = "tmp/Flexilinea_tree/drep/dereplicated_genomes",
        outdir = "tmp/Flexilinea_tree/fastani"
    run:
        taxids = pd.read_csv("05-results/PHYL_Flexilinea_proteintree_taxa.tsv", sep="\t")
        accids = taxids.loc[taxids['NCBI taxonomy Id'] \
                .isin([1946754, 1946750, 2699749, 2017390,
                       1678840, 1946749, 1889813])]['accession Id'].tolist()
        with open(output[0], "wt") as outfile:
            for fn in glob(f"{params.amag}/*.fa.gz"):
                outfile.write(f"{params.outdir}/{os.path.basename(fn).replace('.fa.gz', '.fa')}\n")
            for fn in glob(f"{params.refgenomes}/*.fna"):
                if os.path.basename(fn).replace(".fna", "") in accids:
                    outfile.write(f"{params.outdir}/{os.path.basename(fn).replace('.fna', '.fa')}\n")

rule filter_fasta:
    output:
        "tmp/Flexilinea_tree/fastani/{genome}.fa"
    message: "Discard contigs < 2 kb: {wildcards.genome}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 8
    params:
        fa = lambda wildcards: f"tmp/Flexilinea_tree/drep/dereplicated_genomes/{wildcards.genome}.fna" if wildcards.genome.startswith("GC") else f"05-results/genomes/flexilinea/{wildcards.genome}.fa.gz"
    shell:
        """
        bioawk -c fastx '{{if (length($seq) >= 2000){{print ">" $name "\\n" $seq}}}}' {params.fa} > {output}
        """

rule fastani:
    input:
        fas = lambda wildcards: return_fasta_fn(wildcards),
        lst = "04-analysis/phylogenetics/Flexilinea/flexilinea_fasta.txt"
    output:
        "05-results/PHYL_Flexilinea_fastANI.tsv"
    message: "Calculate the pairwise ANI for the ancient MAGs and the neighbouring clades"
    conda: "ENVS_fastani.yaml"
    threads: 16
    shell:
        """
        fastANI --ql {input.lst} --rl {input.lst} -t {threads} -k 16 --matrix -o {output}
        """

################################################################################

#### Calculate the bootstrap analysis for the same subset as fastANI ###########

checkpoint raxml_samplelist:
    input:
        "05-results/PHYL_Flexilinea_proteintree_taxa.tsv"
    output:
        touch("04-analysis/phylogenetics/Flexilinea/raxml_bootstrap_subset.done")
    message: "Prepare the input list for fastANI from all sequences belonging to the genus Flexilinea"
    params:
        amag = "05-results/genomes/flexilinea",
        refgenomes = "tmp/Flexilinea_tree/drep/dereplicated_genomes",
        outdir = "tmp/Flexilinea_tree/bootstrap"
    run:
        # Prepare reference genomes
        taxids = pd.read_csv(input[0], sep="\t")
        accids = taxids.loc[taxids['NCBI taxonomy Id'] \
                .isin([1946754, 1946750, 2699749, 2017390,
                       1678840, 1946749, 1889813])]['accession Id'].tolist()
        accids.append("GCA_004294845.1")
        os.makedirs(params.outdir)
        for accid in accids:
            os.symlink(f"{os.getcwd()}/{params.refgenomes}/{accid}.fna",
                       f"{params.outdir}/{accid}.fa")
        # Prepare MAGs
        for fn in glob(f"{params.amag}/*.fa.gz"):
            accid = os.path.basename(fn).replace(".fa.gz", "")
            with open(f"{params.outdir}/{accid}.fa", "wb") as outfile:
                with gzip.open(fn, "rb") as fastafile:
                    for line in fastafile:
                        outfile.write(line)

rule phylophlan3_bootstrap_tree:
    input:
        database = "tmp/Flexilinea_tree/core_genes/s__Flexilinea_flocculi/s__Flexilinea_flocculi.faa",
        config = "04-analysis/phylogenetics/Flexilinea/config.cfg",
        genomes = "04-analysis/phylogenetics/Flexilinea/raxml_bootstrap_subset.done"
    output:
        "04-analysis/phylogenetics/Flexilinea/bootstrap/RAxML_bestTree.bootstrap_refined.tre"
    message: "Run PhyloPhlAn3 to generate a tree of the genus Flexilinea"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        genomes = "tmp/Flexilinea_tree/bootstrap",
        db_folder = "tmp/Flexilinea_tree/core_genes",
        outdir = "04-analysis/phylogenetics/Flexilinea/bootstrap"
    log: "04-analysis/phylogenetics/Flexilinea/bootstrap/phylophlan.log"
    threads: 36
    shell:
        """
        phylophlan \
            -i {params.genomes} \
            -d s__Flexilinea_flocculi \
            --databases_folder {params.db_folder} \
            --diversity medium \
            --accurate \
            --genome_extension .fa \
            -f {input.config} \
            -o Flexilinea \
            --output_folder {params.outdir} \
            --nproc {threads} \
            --verbose 2>&1 | tee {log}
        for fn in bootstrap_concatenated.aln bootstrap.tre bootstrap_resolved.tre \
                  RAxML_bestTree.bootstrap_refined.tre RAxML_info.bootstrap_refined.tre RAxML_result.bootstrap_refined.tre; do
            mv {params.outdir}/Flexilinea/${{fn}} {params.outdir}/
        done
        """

rule raxml_bootstrap:
    input:
        "04-analysis/phylogenetics/Flexilinea/bootstrap/RAxML_bestTree.bootstrap_refined.tre"
    output:
        "04-analysis/phylogenetics/Flexilinea/bootstrap/RAxML_bipartitionsBranchLabels.Flexilinea_RAxML_wbootstrap.tre"
    message: "Re-run RAxML with bootstrapping"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        tre = f"{os.getcwd()}/04-analysis/phylogenetics/Flexilinea/bootstrap/bootstrap_resolved.tre",
        aln = f"{os.getcwd()}/04-analysis/phylogenetics/Flexilinea/bootstrap/bootstrap_concatenated.aln",
        wdir = f"{os.getcwd()}/04-analysis/phylogenetics/Flexilinea/bootstrap"
    threads: 20
    shell:
        """
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG -p 1989 \
            -b 1 \
            -N 100 \
            -t {params.tre} \
            -w {params.wdir} \
            -s {params.aln} \
            -n bootstrap_refined_bootstrap.tre \
            -T {threads} &&
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG \
            -w {params.wdir} \
            -f b \
            -n Flexilinea_RAxML_wbootstrap.tre \
            -t {params.wdir}/RAxML_bestTree.bootstrap_refined.tre \
            -z {params.wdir}/RAxML_bootstrap.bootstrap_refined_bootstrap.tre
        """

rule copy_raxml_bootstrap:
    input:
        "04-analysis/phylogenetics/Flexilinea/bootstrap/RAxML_bipartitionsBranchLabels.Flexilinea_RAxML_wbootstrap.tre"
    output:
        "05-results/PHYL_Flexilinea_proteintree_RAxML_bootstrap.tre"
    message: "Copy the RAxML tree witht the bootstrap support values to the results"
    shell:
        "cp {input} {output}"

################################################################################
