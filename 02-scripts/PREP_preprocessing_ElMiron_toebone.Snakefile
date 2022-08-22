################################################################################
# Project: Natural products from the Palaeolithic
# Part: Data preparation
# Step: Prepare the metagenomic sequencing data of the toe bone sample
#       collected at El Miron (Fu et al. (2016))
#
# Download the sequencing data of the El Miron toe bone sample from ENA,
# pre-process them using nf-core/eager before extracting the non-human
# sequencing data.
#
# Alex Huebner, 22/08/22
################################################################################

import json

import pandas as pd

#### SAMPLES ###################################################################
SAMPLES = {'ElMiron': ["A5268", "A5279", "A5301"]}  # six IDs
################################################################################

#### Auxilliary functions ######################################################

def return_url(wildcards):
    with open(checkpoints.fetch_ena_info.get(**wildcards).output[0], "rt") as jsonfile:
        ena_rep = json.load(jsonfile)
    return ena_rep[int(wildcards.i) - 1]['url']

################################################################################

rule all:
    input:
        "05-results/PREP_Nextflow_EAGER_noReads_ElMiron_toebone.tsv"

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
        "04-analysis/eager/elmiron_toebone.tsv"
    message: "Generate the EAGER input TSV"
    params:
        seqdatadir = "tmp/ffq"
    run:
        # Prepare table
        table = {'Sample_Name': "ElMiron",
                 'Library_ID': SAMPLES["ElMiron"],
                 'Lane': [1] * 6,
                 'Colour_Chemistry': [2] * 6,
                 'SeqType': ['PE'] * 6,
                 'Organism': 'Homo sapiens',
                 'Strandedness': ['single'] * 6,
                 'UDG_Treatment': ['full'] * 6, 
                 'R1': [f'{os.getcwd()}/{params.seqdatadir}/{ers}_1.fastq.gz'
                        for ers in SAMPLES["ElMiron"]],
                 'R2': [f'{os.getcwd()}/{params.seqdatadir}/{ers}_2.fastq.gz'
                        for ers in SAMPLES["ElMiron"]],
                 'BAM': ['NA'] * 6}
        pd.DataFrame.from_dict(table) \
            .sort_values(['Library_ID']) \
            .to_csv(output[0], sep="\t", index=False)

################################################################################

#### Run nf-core/eager per library type ########################################

rule run_eager:
    input:
        "04-analysis/eager/elmiron_toebone.tsv"
    output:
        touch("04-analysis/eager/elmiron_toebone.done")
    message: "Run nf-core/EAGER to process the sequencing data of the El Miron toe bone"
    params:
        outdir = "04-analysis/eager/elmiron_toebone"
    shell:
        """
        nextflow run nf-core/eager -r 2.4.5 \
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
        "04-analysis/eager/elmiron_toebone.done"
    output:
        pipe("tmp/eager_extract_unmapped/{sample}.nsorted.bam")
    message: "Sort the BAM file by name: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 12,
        cores = 4
    params:
        bam = "04-analysis/eager/elmiron_toebone/mapping/bwa/{sample}_PE.mapped.bam"
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
        pe1 = "tmp/eager_extract_unmapped/{sample}_1.fastq.gz",
        pe2 = "tmp/eager_extract_unmapped/{sample}_2.fastq.gz",
        pe0 = "tmp/eager_extract_unmapped/{sample}_0.fastq.gz"
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

rule concat_fastqs:
    input:
        lambda wildcards: [f"tmp/eager_extract_unmapped/{lib}_{i}.fastq.gz" for lib in LIBRARIES for i in range(3)]
    output:
        pe0 = "03-data/eager_elmirontoebone/ElMiron_0.fastq.gz",
        pe1 = "03-data/eager_elmirontoebone/ElMiron_1.fastq.gz",
        pe2 = "03-data/eager_elmirontoebone/ElMiron_2.fastq.gz"
    message: "Concatenate all libraries"
    params:
        tmpdir = "tmp/eager_extract_unmapped",
        outdir = "03-data/eager_elmirontoebone"
    run:
        for i in range(3):
            with open(f"{params.outdir}/ElMiron_{i}.fastq.gz", "wb") as outfile:
                for lib in LIBRARIES:
                    with open(f"{params.tmpdir}/{lib}_{i}.fastq.gz", "rb") as libfile:
                        for line in libfile:
                            outfile.write(line)

rule count_reads:
    input:
        pe1 = "03-data/eager_elmirontoebone/ElMiron_1.fastq.gz",
        pe2 = "03-data/eager_elmirontoebone/ElMiron_2.fastq.gz",
        pe0 = "03-data/eager_elmirontoebone/ElMiron_0.fastq.gz"
    output:
        "05-results/PREP_Nextflow_EAGER_noReads_ElMiron_toebone.tsv"
    message: "Count the number of reads"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    shell:
        """
        reads_PE1=$(bioawk -c fastx 'END{{print NR}}' {input.pe1})
        reads_PE2=$(bioawk -c fastx 'END{{print NR}}' {input.pe2})
        echo -e "sample\tR1\tR2\nElMiron\t${{reads_PE1}}\t${{reads_PE2}}\n" > {output}
        """
