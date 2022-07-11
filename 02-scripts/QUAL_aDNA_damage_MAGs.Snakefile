################################################################################
# Project: Natural products from the Palaeolithic
# Part: Quality analyses
# Step: Profile the aDNA damage on specific MAGs using damageprofiler
#
# Dependent on:
#   - ASMB_automaticRefinement.Snakefile
#
# Alex Huebner, 10/07/22
################################################################################

import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
CONTIGS = "04-analysis/ancient_metagenome_assembly/alignment/megahit/EMN001-megahit.fasta.gz"
MAGS = ['EMN001-megahit_006',  # Desulfobulbus oralis
        'EMN001-megahit_021',  # Chlorobium limicola
       ]
SAMPLES = ['EMN001', 'PES001', 'PLV001', 'RIG001', 'GOY005']
################################################################################

wildcard_constraints:
    sample = "[A-Z]+[0-9]+",
    mag = "[A-Z]+[0-9]+-[a-z]+_[0-9]+"

rule all:
    input:
        "05-results/QUAL_damageprofile_Climicola_Doralis.tsv"

rule prepare_fasta:
    output:
        "tmp/damageprofiling_mags/EMN001_contigs.fasta.gz.fai"
    message: "Copy and index FastA file with all EMN001 contigs"
    conda: "ENVS_samtools.yaml"
    params:
        fa = "tmp/damageprofiling_mags/EMN001_contigs.fasta.gz"
    shell:
        """
        cp {CONTIGS} {params.fa}
        samtools faidx {params.fa}
        """

rule bowtie2_index_contigs:
    input:
        "tmp/damageprofiling_mags/EMN001_contigs.fasta.gz.fai"
    output:
        "tmp/damageprofiling_mags/EMN001_contigs.1.bt2"
    message: "Index MEGAHIT assembled contigs"
    conda: "ENVS_bowtie2.yaml"
    resources:
        mem = 8,
        cores = 8
    params:
        fa = "tmp/damageprofiling_mags/EMN001_contigs.fasta.gz",
        prefix = "tmp/damageprofiling_mags/EMN001_contigs"
    threads: 8
    shell:
        """
        bowtie2-build --threads {threads} \
                {params.fa} {params.prefix}
        """

rule bowtie2:
    input:
        "tmp/damageprofiling_mags/EMN001_contigs.1.bt2"
    output:
        pipe("tmp/damageprofiling_mags/{sample}.sam")
    message: "Align sequences against reference genomes using BowTie2: {wildcards.sample}"
    conda: "ENVS_bowtie2.yaml"
    resources:
        mem = 16,
        cores = 16
    params:
        index = "tmp/damageprofiling_mags/EMN001_contigs",
        pe1 = "03-data/eager_fastqs/{sample}-nonUDG_1.fastq.gz",
        pe2 = "03-data/eager_fastqs/{sample}-nonUDG_2.fastq.gz"
    threads: 16
    shell:
        """
        bowtie2 -p {threads} --very-sensitive -N 1 -x {params.index} \
            -1 {params.pe1} -2 {params.pe2} > {output}
        """

rule sam2bam:
    input:
        "tmp/damageprofiling_mags/{sample}.sam"
    output:
        "tmp/damageprofiling_mags/{sample}.bam"
    message: "Convert SAM to BAM format: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    shell:
        "samtools view -Su {input} > {output}"

rule samtools_fixmate:
    input:
        "tmp/damageprofiling_mags/{sample}.bam"
    output:
        pipe("tmp/damageprofiling_mags/{sample}.fixmate.bam")
    message: "Fix mate flags: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    shell:
        "samtools fixmate -mu {input} {output}"

rule samtools_sort:
    input:
        "tmp/damageprofiling_mags/{sample}.fixmate.bam"
    output:
        pipe("tmp/damageprofiling_mags/{sample}.sorted.bam")
    message: "Sort BAM file by coordinate: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 6
    shell:
        "samtools sort -u -o {output} {input}"

rule samtools_calmd:
    input:
        "tmp/damageprofiling_mags/{sample}.sorted.bam"
    output:
        "tmp/damageprofiling_mags/{sample}.calmd.bam"
    message: "Calculate the MD tag: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    params:
        fa = CONTIGS
    shell:
        "samtools calmd -b {input} {params.fa} > {output}"

rule samtools_index:
    input:
        "tmp/damageprofiling_mags/{sample}.calmd.bam"
    output:
        "tmp/damageprofiling_mags/{sample}.calmd.bam.bai"
    message: "Index the BAM file: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    shell:
        "samtools index {input}"

rule list_contigs:
    output:
        "tmp/damageprofiling_mags/{mag}.contiglist.txt"
    message: "List all contigs: {wildcards.mag}"
    conda: "ENVS_bioawk.yaml"
    params:
        mag = lambda wildcards: f"04-analysis/automatic_MAG_refinement/aDNA_samples_human/{wildcards.mag.split('_')[0]}/bins/{wildcards.mag}.fasta.gz"
    shell:
        """
        bioawk -c fastx '{{print $name}}' {params.mag} > {output}
        """

rule subset_to_mag:
    input:
        contigs = "tmp/damageprofiling_mags/{mag}.contiglist.txt",
        bam = "tmp/damageprofiling_mags/{sample}.calmd.bam",
        bai = "tmp/damageprofiling_mags/{sample}.calmd.bam.bai"
    output:
        "tmp/damageprofiling_mags/{sample}-{mag}.bam"
    message: "Subset the BAM file to contigs of the MAG: {wildcards.mag} for sample {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    shell:
        """
        cat {input.contigs} | xargs samtools view -bh {input.bam} > {output}
        """

rule decompress_fasta:
    output:
        "tmp/damageprofiling_mags/EMN001.fasta"
    message: "Decompress FastA file of all EMN001 contigs"
    params:
        fa = CONTIGS
    shell:
        "gunzip -c {params.fa} > {output}"

rule samtools_faidx:
    input:
        "tmp/damageprofiling_mags/EMN001.fasta"
    output:
        "tmp/damageprofiling_mags/EMN001.fasta.fai"
    message: "Index FastA file of all EMN001 contigs"
    conda: "ENVS_samtools.yaml"
    shell:
        "samtools faidx {input}"

rule damageprofiler:
    input:
        bam = "tmp/damageprofiling_mags/{sample}-{mag}.bam",
        fa = "tmp/damageprofiling_mags/EMN001.fasta",
        faidx = "tmp/damageprofiling_mags/EMN001.fasta.fai"
    output:
        "tmp/damageprofiling_mags/{sample}-{mag}/5p_freq_misincorporations.txt"
    message: "Profile the aDNA damage using damageprofiler: {wildcards.mag} for sample {wildcards.sample}"
    conda: "ENVS_damageprofiler.yaml"
    resources:
        mem = 8
    params:
        outdir = "tmp/damageprofiling_mags/{sample}-{mag}"
    shell:
        """
        damageprofiler -i {input.bam} \
            -o {params.outdir} \
            -r {input.fa}
        """

rule summarise_damageprofiler:
    input:
        expand("tmp/damageprofiling_mags/{sample}-{mag}/5p_freq_misincorporations.txt", sample=SAMPLES, mag=MAGS)
    output:
        "05-results/QUAL_damageprofile_Climicola_Doralis.tsv"
    message: "Summarise the substitution frequency at the 5' end"
    run:
        damage = pd.concat([pd.read_csv(fn, sep="\t", skiprows=3) \
                                .assign(sample_MAG=os.path.basename(os.path.dirname(fn)))
                            for fn in input])
        damage['sample'] = damage['sample_MAG'].str.split("-").str[0]
        damage['MAG'] = damage['sample_MAG'].str.split("_").str[1]
        damage['MAG'] = ["Doralis" if m == "006" else "Climicola" for m in damage['MAG'].tolist()]
        damage = damage.drop(["sample_MAG"], axis=1)
        damage.iloc[:, [-2, -1] + list(range(0, 13))] \
            .to_csv(output[0], sep="\t", index=False, float_format="%.4f")
