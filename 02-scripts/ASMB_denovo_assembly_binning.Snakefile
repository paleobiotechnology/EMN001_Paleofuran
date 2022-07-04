################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: Preparation of the config files for the de novo assembly pipeline
#
# Prepares the config files to run the Snakemake pipeline
# https://github.com/alexhbnr/ancient_metagenome_assembly that combines the
# features of a classic de novo assembly pipeline for metagenomic samples with
# some additional steps that are customised for ancient DNA samples.
#
# Due to a different choice in the assembler, MEGAHIT for the ancient dental
# calculus samples, metaSPAdes for the modern ones, we will run the pipeline on
# these subsets separately.
#
# Dependent on:
#   - PREP_preprocessing_dentalcalculus_sequencing_data.Snakefile
#
# Alex Huebner, 12/06/22
################################################################################

from glob import glob
import os
import yaml

import pandas as pd

#### SAMPLES ###################################################################
SAMPLES = pd.read_csv("01-resources/overview_sequencingdata.tsv", sep="\t")['individualId'].tolist()
################################################################################

rule all:
    input:
        "05-results/ASMB_assemblystats_calN50_metaQUAST.tsv",
        "05-results/ASMB_MAGS_metaWRAP_prefilter.tsv"

#### Prepare input files for running pipeline ##################################

rule prepare_input_table:
    output:
        aDNA_human = "04-analysis/ancient_metagenome_assembly-aDNA_samples.tsv",
        mDNA = "04-analysis/ancient_metagenome_assembly-mDNA_samples.tsv"
    message: "Write the sample tables for the different assembly batches"
    params:
        samples = "01-resources/overview_sequencingdata.tsv",
        dir = "03-data/eager_fastqs"
    run:
        lib_overview = pd.read_csv(params.samples, sep="\t")

        sampletsv = []
        for sample in lib_overview.itertuples():
            alldata_fns = sorted(glob(f"{os.getcwd()}/{params.dir}/{sample.individualId}-alldata_*.fastq.gz"))
            nonudg_fns = sorted(glob(f"{os.getcwd()}/{params.dir}/{sample.individualId}-nonUDG_*.fastq.gz"))
            if len(alldata_fns) == 2 and len(nonudg_fns) == 2:
                sampletsv.append((sample.individualId, alldata_fns[0], alldata_fns[1], "",
                                  nonudg_fns[0], nonudg_fns[1], ""))
            elif len(alldata_fns) == 2 and len(nonudg_fns) == 0:
                sampletsv.append((sample.individualId, alldata_fns[0], alldata_fns[1], "",
                                  alldata_fns[0], alldata_fns[1], ""))
            elif len(alldata_fns) == 3 and len(nonudg_fns) == 3:
                sampletsv.append((sample.individualId, alldata_fns[1], alldata_fns[2], alldata_fns[0],
                                  nonudg_fns[1], nonudg_fns[2], nonudg_fns[0]))
            elif len(alldata_fns) == 3 and len(nonudg_fns) == 2:
                sampletsv.append((sample.individualId, alldata_fns[1], alldata_fns[2], alldata_fns[0],
                                  nonudg_fns[0], nonudg_fns[1], ""))

        sampletsv = pd.DataFrame(sampletsv, columns=['sample', 'R1', 'R2', 'R0', 'nonUDG_R1', 'nonUDG_R2', 'nonUDG_R0'])
        sampletsv = sampletsv.drop_duplicates()
        # Write to file
        ## aDNA human
        adna_humans = lib_overview.loc[(~lib_overview['individualId'].str.startswith("VLC")) &
                                       (~lib_overview['individualId'].str.startswith("JAE")), 'individualId'].tolist()
        sampletsv.loc[sampletsv['sample'].isin(adna_humans)] \
            .to_csv(output.aDNA_human, sep="\t", index=False)
        ## mDNA
        mdna_humans = lib_overview.loc[(lib_overview['individualId'].str.startswith("VLC")) |
                                       (lib_overview['individualId'].str.startswith("JAE")), 'individualId'].tolist()
        sampletsv.loc[sampletsv['sample'].isin(mdna_humans)] \
            .to_csv(output.mDNA, sep="\t", index=False)

rule write_config:
    output:
        aDNA_human = "04-analysis/ancient_metagenome_assembly-aDNA_samples.config",
        mDNA = "04-analysis/ancient_metagenome_assembly-mDNA_samples.config"
    message: "Prepare the YAML config files"
    params:
        yaml = "/mnt/archgen/users/huebner/ancient_metagenome_assembly/config/config.yaml",
        pipelinedir = "/mnt/archgen/users/huebner/ancient_metagenome_assembly"
    run:
        with open(f"{params.pipelinedir}/config/config.yaml", "r") as f:
            config = yaml.safe_load(f)
        config['tmpdir'] = f"{os.getcwd()}/tmp/ancient_metagenome_assembly"
        config['resultdir'] = f"{os.getcwd()}/04-analysis/ancient_metagenome_assembly"
        # Write aDNA configs
        config['sampletsv'] = f"{os.getcwd()}/04-analysis/ancient_metagenome_assembly-aDNA_samples.tsv"
        with open("04-analysis/ancient_metagenome_assembly-aDNA_samples.config", "w", encoding="utf8") as outfile:
            yaml.dump(config, outfile, default_flow_style=False, allow_unicode=True)
        # Write mDNA config
        config['sampletsv'] = f"{os.getcwd()}/04-analysis/ancient_metagenome_assembly-mDNA_samples.tsv"
        config['assembler'] = "metaspades"
        config['readcorrection'] = True
        config['bowtie2_seed_mismatches'] = 0
        with open(f"04-analysis/ancient_metagenome_assembly-mDNA_samples.config", "w", encoding="utf8") as outfile:
            yaml.dump(config, outfile, default_flow_style=False, allow_unicode=True)

################################################################################

#### Run assembly pipeline #####################################################

rule run_assembly_pipeline:
    input:
        tsv = "04-analysis/ancient_metagenome_assembly-{sampletype}_samples.tsv",
        config = "04-analysis/ancient_metagenome_assembly-{sampletype}_samples.config"
    output:
        "04-analysis/ancient_metagenome_assembly/summary_{sampletype}_samples.tsv"
    message: "Run the assembly pipeline on {wildcards.sampletype} samples"
    params:
        path_to_pipeline = "../../ancient_metagenome_assembly"
    shell:
        """
        cd {params.path_to_pipeline}
        snakemake --configfile {input.config} \
            --use-conda --conda-prefix conda --profile sge -j 10 -k
        """

################################################################################

#### Evaluate performance ######################################################

rule no_bins_per_binner:
    input:
        expand("04-analysis/ancient_metagenome_assembly/summary_{sampletype}_samples.tsv", sampletype=['aDNA', 'mDNA'])
    output:
        temp("tmp/binner_stats/{sample}.nbins.txt")
    message: "Count the number of bins per binner: {wildcards.sample}"
    params:
        dir = "04-analysis/ancient_metagenome_assembly/binning/metawrap/INITIAL_BINNING"
    run:
        with open(output[0], "wt") as outfile:
            s = wildcards.sample
            if not wildcards.sample.startswith("JAE") and not wildcards.sample.startswith("VLC"):
                assembler = "megahit"
            else:
                assembler = "metaspades"
            for binner in ['metabat2', 'maxbin2', 'concoct']:
                n = sum([fn.split(".")[1].isnumeric()
                         for fn in glob(f"{params.dir}/{wildcards.sample}-{assembler}/{binner}_bins/bin.*.fa")])
                s += "\t" + str(n)
            outfile.write(s + "\n")

rule summarise_nbins:
    input:
        expand("tmp/binner_stats/{sample}.nbins.txt", sample=SAMPLES)
    output:
        "05-results/ASMB_nBins_initialbinning.tsv"
    message: "Combine the number of bin results into a single table"
    run:
        pd.concat([pd.read_csv(fn, sep='\t', header=None,
                               names=['sample', 'metabat2', 'maxbin2', 'concoct'])
                   for fn in input]) \
            .sort_values(['sample']) \
            .drop_duplicates() \
            .to_csv(output[0], sep="\t", index=False)

################################################################################

#### Summary ###################################################################

rule summarise_assembly_stats:
    input:
        expand("04-analysis/ancient_metagenome_assembly/summary_{sampletype}_samples.tsv", sampletype=['aDNA', 'mDNA'])
    output:
        "05-results/ASMB_assemblystats_calN50_metaQUAST.tsv"
    message: "Summarise the assembly stats generated by calN50 and metaQUAST"
    params:
        nreads = "05-results/PREP_Nextflow_EAGER_noReads_per_sample.tsv"
    run:
        nreads = pd.read_csv(params.nreads, sep="\t") \
            .query("sample.str.contains('alldata')")
        nreads['total'] = nreads['R1'] + nreads['R0']
        excluded_samples = nreads.loc[nreads['total'] < 5e6]['sample'] \
            .str.split("-").str[0].tolist()
        pd.concat([pd.read_csv(fn, sep='\t')
                             for fn in input]) \
            .query("~sample.isin(@excluded_samples)") \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False, float_format="%0.3f")

rule summarise_binning_stats:
    input:
        nbins = "05-results/ASMB_nBins_initialbinning.tsv",
        tsv = expand("04-analysis/ancient_metagenome_assembly/summary_{sampletype}_samples.tsv", sampletype=['aDNA', 'mDNA'])
    output:
        "05-results/ASMB_MAGS_metaWRAP_prefilter.tsv"
    message: "Summarise the list of bins after the refinement using MetaWRAP"
    params:
        qual = "04-analysis/ancient_metagenome_assembly/stats/bin_quality",
    run:
        # Summarise the overview of the bin quality (checkM, GUNC, BUSCO)
        binqual = pd.concat([pd.read_csv(fn, sep="\t")
                             for fn in glob(f"{params.qual}/*.tsv")
                             if os.stat(fn).st_size > 28])

        binqual.to_csv(output[0], sep="\t", index=False, float_format="%.2f")

################################################################################
