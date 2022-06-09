################################################################################
# Project: Natural products from the Palaeolithic
# Part: Data preparation
# Step: Pre-process the sequencing data using nf-core/eager
#
# This script prepares the files required to pre-process the sequencing data
# with the Nextflow pipeline nf-core/eager and subsequently merges the data of
# an individual. Due to the set-up of nf-core/eager, we will process ssDNA and
# dsDNA libraries in separate runs. Furthermore, we will combine the sequencing
# data of each individual in two ways, either using all sequencing data or just
# the non-UDG treated one. The latter will be used for the analysis of ancient
# DNA damage.
#
# Alex Huebner, 09/06/22
################################################################################

import os

import pandas as pd


if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
LIB_OVERVIEW = pd.read_csv("01-resources/overview_sequencingdata.tsv", sep="\t")
LIB_OVERVIEW['sequencingType'] = LIB_OVERVIEW['sequencingSetup'].str[:2]
LIBRARIES_SSDNA = [f'{lib.libraryId}_{lib.sequencingType}'
                   for lib in LIB_OVERVIEW.loc[LIB_OVERVIEW['libraryType'] == "ssDNA"].itertuples()]
SAMPLES = [f"{s}-alldata" for s in LIB_OVERVIEW['individualId'].unique()]
SAMPLES.extend([f"{s}-nonUDG"
                for s in LIB_OVERVIEW.loc[LIB_OVERVIEW['libraryTreatment'] == "none", 'individualId'].unique()])
LIBRARIES = [f"{lib.libraryId}_{lib.sequencingType}"
             for lib in LIB_OVERVIEW.itertuples()]
################################################################################

rule all:
    input:
        "05-results/PREP_Nextflow_EAGER_noReads_per_sample.tsv"

#### Prepare input tables for nf-core/eager ####################################

rule generate_eager_tsv:
    output:
        pe = "04-analysis/eager/eager_input_pe.tsv",
        se = "04-analysis/eager/eager_input_se.tsv"
    message: "Generate the EAGER input TSV"
    run:
        eager = LIB_OVERVIEW[['individualId', 'libraryId']] \
            .rename({'individualId': 'Sample_Name',
                     'libraryId': 'Library_ID'}, axis=1)
        # Add infos about the library preparation and sequencing
        eager['Lane'] = 1
        eager['Colour_Chemistry'] = [4 if platform == "HiSeq4000" else 2
                                     for platform in LIB_OVERVIEW['sequencingPlatform']]
        eager['SeqType'] = LIB_OVERVIEW['sequencingType']
        eager['Organism'] = "NA"
        eager['Strandedness'] = ["double" if libtype == "dsDNA" else "single"
                                 for libtype in LIB_OVERVIEW['libraryType']]
        eager['UDG_Treatment'] = [treatment.replace("-UDG", "")
                                  for treatment in LIB_OVERVIEW['libraryTreatment']]
        # Add infos about the sequencing
        seqdata = []
        for library in LIB_OVERVIEW.itertuples():
            if library.sequencingType == "PE":
                seqdata.append((f"{os.getcwd()}/03-data/raw_data/{library.run_accession}_1.fastq.gz",
                                f"{os.getcwd()}/03-data/raw_data/{library.run_accession}_2.fastq.gz"))
            elif library.sequencingType == "SE":
                seqdata.append((f"{os.getcwd()}/03-data/raw_data/{library.run_accession}_0.fastq.gz",
                                "NA"))
        seqdata_df = pd.DataFrame(seqdata, columns=['R1', 'R2'])
        eager = pd.concat([eager, seqdata_df], axis=1)
        eager['BAM'] = "NA"

        eager = eager[['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType',
               'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']] \
            .sort_values(['Sample_Name', 'Library_ID'])

        eager.loc[eager['Strandedness'] == "double"] \
            .to_csv(output["pe"], sep="\t", index=False)
        eager.loc[eager['Strandedness'] == "single"] \
            .to_csv(output["se"], sep="\t", index=False)

################################################################################

#### Run nf-core/eager per library type ########################################

rule run_eager_pe:
    input:
        pe = "04-analysis/eager/eager_input_pe.tsv",
        se = "04-analysis/eager/eager_input_se.tsv",
    output:
        touch("04-analysis/eager/eager_pe.done")
    message: "Run nf-core/EAGER to process the dsDNA sequencing data"
    params:
        outdir = "04-analysis/eager/results_pe"
    shell:
        """
        nextflow run nf-core/eager -r 2.4.2 \
            -profile eva,archgen \
            --input "{input.pe}" \
            --fasta '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa' \
            --fasta_index '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa.fai' \
            --bwa_index '/mnt/archgen/Reference_Genomes/Human/hs37d5' \
            --seq_dict '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.dict' \
            --skip_deduplication --skip_damage_calculation \
            --complexity_filter_poly_g \
            --skip_collapse \
            --outdir "{params.outdir}"
        """

rule run_eager_se:
    input:
        pe = "04-analysis/eager/eager_input_pe.tsv",
        se = "04-analysis/eager/eager_input_se.tsv",
    output:
        touch("04-analysis/eager/eager_se.done")
    message: "Run nf-core/EAGER to process the ssDNA sequencing data"
    params:
        outdir = "04-analysis/eager/results_se"
    shell:
        """
        nextflow run nf-core/eager -r 2.4.2 \
            -profile eva,archgen \
            --input "{input.se}" \
            --fasta '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa' \
            --fasta_index '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa.fai' \
            --bwa_index '/mnt/archgen/Reference_Genomes/Human/hs37d5' \
            --seq_dict '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.dict' \
            --skip_deduplication --skip_damage_calculation \
            --complexity_filter_poly_g \
            --skip_collapse \
            --clip_reverse_adaptor GGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
            --outdir "{params.outdir}"
        """

################################################################################

#### Extract the unaligned reads in FastQ format ###############################

rule samtools_sort_by_name:
    input:
        eager_pe = "04-analysis/eager/eager_pe.done",
        eager_se = "04-analysis/eager/eager_se.done"
    output:
        pipe("tmp/eager_extract_unmapped/{lib}.nsorted.bam")
    message: "Sort the BAM file by name: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 12,
        cores = 4
    params:
        bam = lambda wildcards: f"04-analysis/eager/results_se/mapping/bwa/{wildcards.lib}.mapped.bam" if wildcards.lib in LIBRARIES_SSDNA else f"04-analysis/eager/results_pe/mapping/bwa/{wildcards.lib}.mapped.bam"
    threads: 4
    shell:
        """
        samtools sort -n -o {output} {params.bam}
        """

rule samtools_fixmate:
    input:
        "tmp/eager_extract_unmapped/{lib}.nsorted.bam"
    output:
        pipe("tmp/eager_extract_unmapped/{lib}.fixmate.bam")
    message: "Fix mate flags: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        samtools fixmate -mcu -@ {threads} {input} {output}
        """

rule extract_unmapped_reads:
    input:
        "tmp/eager_extract_unmapped/{lib}.fixmate.bam"
    output:
        pe1 = "tmp/eager_extract_unmapped/{lib}_1.fastq.gz",
        pe2 = "tmp/eager_extract_unmapped/{lib}_2.fastq.gz",
        pe0 = "tmp/eager_extract_unmapped/{lib}_0.fastq.gz"
    message: "Extract all reads for which are not aligned in a proper pair and convert to fastq: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 2
    threads: 2
    shell:
        """
        samtools view -uh -e '(flag.paired && (flag.unmap || flag.munmap)) || (!flag.paired && flag.unmap)' {input} | \
        samtools fastq -1 {output.pe1} \
                       -2 {output.pe2} \
                       -0 {output.pe0} -
        """

rule concat_fastqs:
    input:
        lambda wildcards: [f"tmp/eager_extract_unmapped/{lib}_{i}.fastq.gz" for lib in [l for l in LIBRARIES if l.split(".")[0] == wildcards.sample.split("-")[0]] for i in range(3)]
    output:
        touch("03-data/eager_fastqs/{sample}.concatenated")
    message: "Concatenate the FastQ files: {wildcards.sample}"
    params:
        outdir = "03-data/eager_fastqs"
    run:
        # Identify libraries for merging
        sample_name, datatype = wildcards.sample.split("-")
        sample_overview = LIB_OVERVIEW.loc[LIB_OVERVIEW['individualId'] == sample_name]
        sample_overview = sample_overview.set_index(['individualId'])
        if datatype == "nonUDG":
            sample_overview = sample_overview.loc[sample_overview['libraryTreatment'] == "none"]

        # Create file lists
        filelist = {'SE': [], 'PE': []}
        for lib in sample_overview.itertuples():
            filelist[lib.sequencingType].append(f"{lib.libraryId}_{lib.sequencingType}")

        # Concatenate
        if len(filelist['SE']) > 0:
            with open(f"{params.outdir}/{wildcards.sample}_0.fastq.gz", "wb") as outfile:
                for lib in filelist['SE']:
                    with open(f"tmp/eager_extract_unmapped/{lib}_0.fastq.gz", "rb") as fastqfile:
                        for line in fastqfile:
                            outfile.write(line)
        if len(filelist['PE']) > 0:
            for i in range(1, 3):
                with open(f"{params.outdir}/{wildcards.sample}_{i}.fastq.gz", "wb") as outfile:
                    for lib in filelist['PE']:
                        with open(f"tmp/eager_extract_unmapped/{lib}_{i}.fastq.gz", "rb") as fastqfile:
                            for line in fastqfile:
                                outfile.write(line)

rule count_reads:
    input:
        "03-data/eager_fastqs/{sample}.concatenated"
    output:
        temp("03-data/eager_fastqs/{sample}.n")
    message: "Count the number of reads: {wildcards.sample}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    params:
        pe1 = "03-data/eager_fastqs/{sample}_1.fastq.gz",
        pe2 = "03-data/eager_fastqs/{sample}_2.fastq.gz",
        pe0 = "03-data/eager_fastqs/{sample}_0.fastq.gz"
    shell:
        """
        reads_PE1=$(bioawk -c fastx 'END{{print NR}}' {params.pe1})
        reads_PE2=$(bioawk -c fastx 'END{{print NR}}' {params.pe2})
        if [[ -f {params.pe0} ]]; then
            reads_PE0=$(bioawk -c fastx 'END{{print NR}}' {params.pe0})
        else
            reads_PE0=0
        fi
        echo -e "{wildcards.sample}\t${{reads_PE1}}\t${{reads_PE2}}\t${{reads_PE0}}" > {output}
        """

rule summarise_count_reads:
    input:
        expand("03-data/eager_fastqs/{sample}.n", sample=SAMPLES)
    output:
        "05-results/PREP_Nextflow_EAGER_noReads_per_sample.tsv"
    message: "Summarise the number of reads per sample"
    run:
        pd.concat([pd.read_csv(fn, sep="\t", header=None, names=['sample', 'R1', 'R2', 'R0'])
                   for fn in input]) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)

################################################################################
