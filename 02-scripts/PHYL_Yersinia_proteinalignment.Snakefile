################################################################################
# Project: Natural products from the Palaeolithic
# Part: Phylogenetic analyses
# Step: Protein-alignment tree of the genus Yersinia using PhyloPhlAn3
#
# Alex Huebner, 27/07/22
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("01-resources/yersinia_genomes/{sample}.fa.gz")
################################################################################

rule all:
    input:
        "05-results/PHYL_Yersinia_proteintree_RAxML_bootstrap.tre",
        "05-results/PHYL_Yersinia_fastANI.tsv"

#### Phylogenetic analysis #####################################################

rule phylophlan3_write_config:
    output:
        "04-analysis/phylogenetics/Yersinia/config.cfg"
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

rule setup_yersinia_coregenes:
    output:
        "tmp/Yersinia_tree/core_genes/s__Yersinia_pestis/s__Yersinia_pestis.faa"
    message: "Download the core genes of Yersinia pestis and generate a database"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        dir = "tmp/Yersinia_tree/core_genes"
    shell:
        """
        mkdir -p {params.dir} && \
        phylophlan_setup_database \
            -g s__Yersinia_pestis --database_update \
            -o {params.dir} \
            --verbose
        """

rule phylophlan3_tree:
    input:
        database = "tmp/Yersinia_tree/core_genes/s__Yersinia_pestis/s__Yersinia_pestis.faa",
        config = "04-analysis/phylogenetics/Yersinia/config.cfg"
    output:
        "04-analysis/phylogenetics/Yersinia/RAxML_bestTree.yersinia_genomes_refined.tre"
    message: "Run PhyloPhlAn3 to generate a tree of the genus Yersinia"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        genomes = "01-resources/yersinia_genomes",
        db_folder = "tmp/Yersinia_tree/core_genes",
        outdir = "04-analysis/phylogenetics/Yersinia"
    log: "04-analysis/phylogenetics/Yersinia/phylophlan.log"
    threads: 36
    shell:
        """
        mkdir -p {params.genomes}/phylophlan
        phylophlan \
            -i {params.genomes} \
            -d s__Yersinia_pestis \
            --databases_folder {params.db_folder} \
            --diversity medium \
            --accurate \
            --genome_extension .fa \
            -f {input.config} \
            -o Yersinia \
            --output_folder {params.outdir} \
            --nproc {threads} \
            --verbose 2>&1 | tee {log}
        for fn in yersinia_genomes_concatenated.aln yersinia_genomes.tre yersinia_genomes_resolved.tre \
                  RAxML_bestTree.yersinia_genomes_refined.tre RAxML_info.yersinia_genomes_refined.tre RAxML_result.yersinia_genomes_refined.tre; do
            mv {params.outdir}/Yersinia/${{fn}} {params.outdir}/
        done
        """

rule raxml_bootstrap:
    input:
        "04-analysis/phylogenetics/Yersinia/RAxML_bestTree.yersinia_genomes_refined.tre"
    output:
        "04-analysis/phylogenetics/Yersinia/RAxML_bipartitionsBranchLabels.Yersinia_RAxML_wbootstrap.tre"
    message: "Re-run RAxML with bootstrapping"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        tre = f"{os.getcwd()}/04-analysis/phylogenetics/Yersinia/yersinia_genomes_resolved.tre",
        aln = f"{os.getcwd()}/04-analysis/phylogenetics/Yersinia/yersinia_genomes_concatenated.aln",
        wdir = f"{os.getcwd()}/04-analysis/phylogenetics/Yersinia"
    threads: 20
    shell:
        """
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG -p 1989 \
            -b 1 \
            -N 100 \
            -t {params.tre} \
            -w {params.wdir} \
            -s {params.aln} \
            -n yersinia_genomes_refined_bootstrap.tre \
            -T {threads} &&
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG \
            -w {params.wdir} \
            -f b \
            -n Yersinia_RAxML_wbootstrap.tre \
            -t RAxML_bestTree.yersinia_genomes_refined.tre \
            -z RAxML_bootstrap.yersinia_genomes_refined_bootstrap.tre
        """

rule copy_raxml_bootstrap:
    input:
        "04-analysis/phylogenetics/Yersinia/RAxML_bipartitionsBranchLabels.Yersinia_RAxML_wbootstrap.tre"
    output:
        "05-results/PHYL_Yersinia_proteintree_RAxML_bootstrap.tre"
    message: "Copy the RAxML tree witht the bootstrap support values to the results"
    shell:
        "cp {input} {output}"

################################################################################

#### Calculate pairwise ANI ####################################################

rule decompress_fasta:
    output:
        "tmp/fastani_yersinia/{sample}.fa"
    message: "Decompress the FastA file: {wildcards.sample}"
    params:
        fa = "01-resources/yersinia_genomes/{sample}.fa.gz" 
    shell:
        "gunzip -c {params.fa} > {output}"

rule generate_input_list:
    input:
        expand("tmp/fastani_yersinia/{sample}.fa", sample=SAMPLES)
    output:
        "04-analysis/phylogenetics/Yersinia/yersinia_fasta.txt"
    message: "Generate list of FastA files for fastANI"
    params:
        tmpdir = "tmp/fastani_yersinia" 
    run:
        with open(output[0], "wt") as outfile:
            for sample in SAMPLES:
                outfile.write(f"{os.getcwd()}/{params.tmpdir}/{sample}.fa\n")

rule fastani:
    input:
        "04-analysis/phylogenetics/Yersinia/yersinia_fasta.txt"
    output:
        "05-results/PHYL_Yersinia_fastANI.tsv"
    message: "Calculate the pairwise ANI for the ancient MAGs and the neighbouring clades"
    conda: "ENVS_fastani.yaml"
    threads: 16
    shell:
        """
        fastANI --ql {input} --rl {input} -t {threads} -k 16 -o {output}
        """

################################################################################
