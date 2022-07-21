################################################################################
# Project: Natural products from the Palaeolithic
# Part: Quality analyses
# Step: Reference alignment and genome reconstruction of Chlorobiaceae genomes
#
# Dependent on:
#   - PREP_preprocessing_dentalcalculus_sequencing_data.Snakefile
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 22/06/22
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
# Dental claculus samples this study
SAMPLES = {os.path.basename(fn).replace("-alldata.concatenated", ""): fn.replace(".concatenated", "")
           for fn in glob("03-data/eager_fastqs/*-alldata.concatenated")}
# Published dental calculus samples
for fn in glob("03-data/eager_weyrich2017/*_1.fastq.gz"):
    SAMPLES[os.path.basename(fn).replace("_1.fastq.gz", "")] = fn.replace("_1.fastq.gz", "")
EMN001_CONTIGS = "04-analysis/ancient_metagenome_assembly/alignment/megahit/EMN001-megahit.fasta.gz"
EMN001_CHL_MAG = "04-analysis/automatic_MAG_refinement/aDNA_samples_human/EMN001-megahit/bins/EMN001-megahit_021.fasta.gz"
################################################################################

wildcard_constraints:
    sample = "[A-Za-z]+[0-9]+"

localrules: unzip_fasta

rule all:
    input:
        expand("04-analysis/refalignment/{sample}.{reads}", sample=SAMPLES, reads=['n_aligned', 'n_total'])
        "05-results/QUAL_dentalcalculus_Chlorobiaceae_refalignment.tsv",

rule unzip_fasta:
    output:
        "tmp/refalignment_emn001/EMN001.fa"
    message: "Decompress FastA with contigs of EMN001"
    params:
        fasta = EMN001_CONTIGS
    shell:
        "gunzip -c {params.fasta} > {output}"

rule index_contigs:
    input:
        "tmp/refalignment_emn001/EMN001.fa"
    output:
        "tmp/refalignment_emn001/EMN001.fa.1.bt2"
    message: "Index contigs of EMN001 for alignment with BowTie2"
    conda: "ENVS_bowtie2.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        bowtie2-build --threads {threads} {input} {input}
        """

rule align_sequences:
    input:
        "tmp/refalignment_emn001/EMN001.fa.1.bt2"
    output:
        pipe("tmp/refalignment_emn001/{sample}.sam")
    message: "Align sequences against EMN001: {wildcards.sample}"
    conda: "ENVS_bowtie2.yaml"
    resources:
        mem = 8,
        cores = 8
    params:
        pe1 = lambda wildcards: f"{SAMPLES[wildcards.sample]}_1.fastq.gz",
        pe2 = lambda wildcards: f"{SAMPLES[wildcards.sample]}_2.fastq.gz",
        pe0 = lambda wildcards: f"-U {SAMPLES[wildcards.sample]}_0.fastq.gz" if os.path.isfile(f"{SAMPLES[wildcards.sample]}_0.fastq.gz") else "",
        fasta = "tmp/refalignment_emn001/EMN001.fa",
    threads: 16
    shell:
        """
        bowtie2 -p {threads} -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x {params.fasta} \
            --rg-id {wildcards.sample} \
            -1 {params.pe1} -2 {params.pe2} {params.pe0} -S {output}
        """

rule samtools_view:
    input:
        "tmp/refalignment_emn001/{sample}.sam"
    output:
        temp("tmp/refalignment_emn001/{sample}.bam")
    message: "Convert SAM to BAM: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 2
    shell:
        "samtools view -Su {input} > {output}"

rule samtools_sort:
    input:
        "tmp/refalignment_emn001/{sample}.bam"
    output:
        pipe("tmp/refalignment_emn001/{sample}.sorted.bam")
    message: "Sort BAM file by coordinate: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8
    params:
        prefix = "/tmp/{sample}"
    shell:
        "samtools sort -u -T {params.prefix} -o {output} {input}"

rule samtools_calmd:
    input:
        bam = "tmp/refalignment_emn001/{sample}.sorted.bam",
        reffa = "tmp/refalignment_emn001/EMN001.fa"
    output:
        temp("tmp/refalignment_emn001/{sample}.calmd.bam")
    message: "Calculate the MD field: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8
    shell:
        "samtools calmd -u {input.bam} {input.reffa} > {output}"

rule samtools_reheader:
    input:
        "tmp/refalignment_emn001/{sample}.calmd.bam"
    output:
        temp("tmp/refalignment_emn001/{sample}.reheader.bam")
    message: "Add SM tag: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 4
    shell:
        """
        samtools view -H {input} | \
        sed "s/@RG\tID:{wildcards.sample}/@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}/" | \
        samtools reheader - {input} > {output}
        """

rule remove_duplicates:
    input:
        "tmp/refalignment_emn001/{sample}.reheader.bam"
    output:
        "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam"
    message: "Remove duplicated reads: {wildcards.sample}"
    conda: "ENVS_picard.yaml"
    resources:
        mem = 12,
        cores = 4
    log: "04-analysis/refalignment/{sample}.markdup.log"
    threads: 4
    shell:
        """
        picard -Xmx8g MarkDuplicates \
            -I {input} \
            -O {output} \
            -M {log} \
            --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate
        """

rule samtools_index:
    input:
        "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam"
    output:
        "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam.bai"
    message: "Index BAM file: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    shell:
        "samtools index {input}"

rule return_mag_contigs:
    output:
        pipe("tmp/refalignment_emn001/{sample}.contigs")
    message: "Return contigs belonging to EMN001's Chlorobiaceae MAG"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    params:
        mag = EMN001_CHL_MAG
    shell:
        "bioawk -c fastx '{{print $name}}' {params.mag} > {output}"

rule count_no_aligned_reads:
    input:
        bam = "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam",
        bai = "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam.bai",
        contigs = "tmp/refalignment_emn001/{sample}.contigs"
    output:
        "04-analysis/refalignment/{sample}.n_aligned"
    message: "Number of aligned reads against the Chlorobiaceae MAG: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 2
    shell:
        """
        cat {input.contigs} | \
        xargs samtools view -c {input.bam} > {output}
        """

rule total_reads:
    input:
        "04-analysis/refalignment/{sample}.sorted.calmd.markdup.bam"
    output:
        "04-analysis/refalignment/{sample}.n_total"
    message: "Total number of reads: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 2
    shell:
        """
        samtools flagstat {input} | head -1 | cut -f 1 -d" " > {output}
        """

rule summary:
    input:
        expand("04-analysis/refalignment/{sample}.{reads}", sample=[s for s in SAMPLES if s != "EMN001"], reads=['n_aligned', 'n_total'])
    output:
        "05-results/QUAL_dentalcalculus_Chlorobiaceae_refalignment.tsv"
    message: "Summarise the read counts"
    params:
        dir = "04-analysis/refalignment"
    run:
        read_counts = []
        for sample in SAMPLES:
            chlorobiaceae_reads = open(f"{params.dir}/{sample}.n_aligned", "rt") \
                .readline().rstrip()
            total_reads = open(f"{params.dir}/{sample}.n_total", "rt") \
                .readline().rstrip()
            read_counts.append((sample, chlorobiaceae_reads, total_reads))

        df = pd.DataFrame(read_counts, columns=['sample', 'alignedReads', 'totalReads'])
        df['alignedReads'] = df['alignedReads'].astype(int)
        df['totalReads'] = df['totalReads'].astype(int)

        df.sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)
