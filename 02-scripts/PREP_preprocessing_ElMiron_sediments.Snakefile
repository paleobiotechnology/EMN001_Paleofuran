################################################################################
# Project: Natural products from the Palaeolithic
# Part: Data preparation
# Step: Prepare the sequencing data of the sediment samples collected at El
#       Miron
#
# Download the sequencing data of the El Miron sediment samples from ENA,
# pre-process them using nf-core/eager before extracting the non-human
# sequencing data.
#
# Alex Huebner, 09/06/22
################################################################################

import json

import pandas as pd

#### SAMPLES ###################################################################
SAMPLES = {'EMN002': 'ERSX0000001',
           'EMN003': 'ERSX0000002',
           'EMN004': 'ERSX0000003',
           'EMN005': 'ERSX0000004',
           'EMN007': 'ERSX0000005',
           'EMN008': 'ERSX0000006'}
################################################################################

#### Auxilliary functions ######################################################

def return_url(wildcards):
    with open(checkpoints.fetch_ena_info.get(**wildcards).output[0], "rt") as jsonfile:
        ena_rep = json.load(jsonfile)
    return ena_rep[int(wildcards.i) - 1]['url']

################################################################################

rule all:
    input:
        "05-results/PREP_Nextflow_EAGER_noReads_ElMiron_sediments.tsv"

#### Download data from ENA ####################################################

checkpoint fetch_ena_info:
    output:
        "tmp/ffq/{err}.json"
    message: "Download the FTP URL and MD5sum from ENA: {wildcards.err}"
    shell:
        "ffq --ftp {wildcards.err} > {output}"

rule download_fastq_file:
    input:
        "tmp/ffq/{err}.json"
    output:
        "03-data/raw_data/{err}_{i}.fastq.gz"
    message: "Download the FastQ file: {wildcards.err} for read {wildcards.i}"
    params:
        url = lambda wildcards: return_url(wildcards),
        cutdirs = lambda wildcards: return_url(wildcards).count("/") - 1
    shell:
        """
        wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent \
            --reject "index.html*" --cut-dirs={params.cutdirs} -O {output} {params.url}
        """

rule calculate_md5sum:
    input:
        "03-data/raw_data/{err}_{i}.fastq.gz"
    output:
        temp("03-data/raw_data/{err}_{i}.md5")
    message: "Calculate the md5sum: {wildcards.err} for read {wildcards.i}"
    shell:
        "md5sum {input} > {output}"

rule validate_md5sum:
    input:
        lambda wildcards: [f"03-data/raw_data/{wildcards.err}_{i}.md5" for i in [1, 2]]
    output:
        "03-data/raw_data/{err}.validated"
    message: "Validate md5sum: {wildcards.err}"
    run:
        with open(checkpoints.fetch_ena_info.get(**wildcards).output[0], "rt") as jsonfile:
            ena_rep = json.load(jsonfile)

        md5sums = [open(fn, 'rt').readline().split()[0]
                   for fn in input]

        if all([ena_rep[i]['md5'] == md5sums[i] for i in range(len(md5sums))]):
            Path(output[0]).touch()
        else:
            print("The md5sums don't match. The files have the md5sums\n" +
                  ", ".join(md5sums) +
                  "\nwhile ENA reports the following values\n" +
                  ", ".join([f['md5'] for f in ena_rep]))
            sys.exit(1)


################################################################################

#### Prepare input tables for nf-core/eager ####################################

rule generate_eager_tsv_profile:
    input:
        expand("03-data/raw_data/{err}.validated", err=SAMPLES.values())
    output:
        "04-analysis/eager/elmiron_sediments.tsv"
    message: "Generate the EAGER input TSV"
    params:
        seqdatadir = "tmp/ffq"
    run:
        # Prepare table
        table = {'Sample_Name': SAMPLES.keys(),
                 'Library_ID': SAMPLES.values(),
                 'Lane': [1] * 6,
                 'Colour_Chemistry': [2] * 6,
                 'SeqType': ['PE'] * 6,
                 'Organism': 'sediment',
                 'Strandedness': ['double'] * 6,
                 'UDG_Treatment': ['none'] * 6, 
                 'R1': [f'{os.getcwd()}/{params.seqdatadir}/{ers}_1.fastq.gz'
                        for ers in SAMPLES.values()],
                 'R2': [f'{os.getcwd()}/{params.seqdatadir}/{ers}_2.fastq.gz'
                        for ers in SAMPLES.values()],
                 'BAM': ['NA'] * 6}
        pd.DataFrame.from_dict(table) \
            .sort_values(['Sample_Name']) \
            .to_csv(output[0], sep="\t", index=False)

################################################################################

#### Run nf-core/eager per library type ########################################

rule run_eager:
    input:
        "04-analysis/eager/elmiron_sediments.tsv"
    output:
        touch("04-analysis/eager/elmiron_sediments.done")
    message: "Run nf-core/EAGER to process the sequencing data of the El Miron sediments"
    params:
        outdir = "04-analysis/eager/elmiron_sediments"
    shell:
        """
        nextflow run nf-core/eager -r 2.4.3 \
            -profile eva,archgen \
            --input "{input}" \
            --fasta '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa' \
            --fasta_index '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa.fai' \
            --bwa_index '/mnt/archgen/Reference_Genomes/Human/hs37d5' \
            --seq_dict '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.dict' \
            --skip_deduplication --skip_damage_calculation \
            --complexity_filter_poly_g \
            --outdir "{params.outdir}"
        """

################################################################################

#### Extract the unaligned reads in FastQ format ###############################

rule samtools_sort_by_name:
    input:
        "04-analysis/eager/elmiron_sediments.done"
    output:
        pipe("tmp/eager_extract_unmapped/{sample}.nsorted.bam")
    message: "Sort the BAM file by name: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 12,
        cores = 4
    params:
        bam = "04-analysis/eager/elmiron_sediments/mapping/bwa/{sample}_PE.mapped.bam"
    threads: 4
    shell:
        """
        samtools sort -n -o {output} {params.bam}
        """

rule samtools_fixmate:
    input:
        "tmp/eager_extract_unmapped/{sample}.nsorted.bam"
    output:
        pipe("tmp/eager_extract_unmapped/{sample}.fixmate.bam")
    message: "Fix mate flags: {wildcards.sample}"
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
        "tmp/eager_extract_unmapped/{sample}.fixmate.bam"
    output:
        pe1 = "03-data/eager_elmironsediments/{sample}_1.fastq.gz",
        pe2 = "03-data/eager_elmironsediments/{sample}_2.fastq.gz",
        pe0 = "03-data/eager_elmironsediments/{sample}_0.fastq.gz"
    message: "Extract all reads for which are not aligned in a proper pair and convert to fastq: {wildcards.sample}"
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

rule count_reads:
    input:
        pe1 = "03-data/eager_elmironsediments/{sample}_1.fastq.gz",
        pe2 = "03-data/eager_elmironsediments/{sample}_2.fastq.gz",
        pe0 = "03-data/eager_elmironsediments/{sample}_0.fastq.gz"
    output:
        temp("03-data/eager_elmironsediments/{sample}.n")
    message: "Count the number of reads: {wildcards.sample}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    shell:
        """
        reads_PE0=$(bioawk -c fastx 'END{{print NR}}' {input.pe0})
        echo -e "{wildcards.sample}\t${{reads_PE0}}" > {output}
        """

rule summarise_count_reads:
    input:
        expand("03-data/eager_elmironsediments/{sample}.n", sample=SAMPLES.keys())
    output:
        "05-results/PREP_Nextflow_EAGER_noReads_ElMiron_sediments.tsv"
    message: "Summarise the number of reads per sample"
    run:
        pd.concat([pd.read_csv(fn, sep="\t", header=None, names=['sample', 'R0'])
                   for fn in input]) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)
