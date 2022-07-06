################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: Write the config files for running the automatic refinement pipeline
#       based on the taxonomic assignment of contigs against the GTDB using 
#       MMSeqs2
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 06/07/22
################################################################################

from glob import glob
import os
import yaml

import pandas as pd

#### SAMPLES ###################################################################
# ADNA samples
ASAMPLES = {os.path.basename(os.path.dirname(fn)): f"{os.getcwd()}/{fn}"
            for fn in glob("04-analysis/ancient_metagenome_assembly/binning/metawrap/BIN_REFINEMENT/**/metawrap_50_10_bins.stats")
            if not os.path.basename(os.path.dirname(fn)).startswith("JAE") and not os.path.basename(os.path.dirname(fn)).startswith("VLC")}
# MDNA samples
MSAMPLES = {os.path.basename(os.path.dirname(fn)): f"{os.getcwd()}/{fn}"
            for fn in glob("04-analysis/ancient_metagenome_assembly/binning/metawrap/BIN_REFINEMENT/**/metawrap_50_10_bins.stats")
            if os.path.basename(os.path.dirname(fn)).startswith("JAE") or os.path.basename(os.path.dirname(fn)).startswith("VLC")}
SAMPLETYPE = {'aDNA': ASAMPLES, 'mDNA': MSAMPLES}
################################################################################

rule all:
    input:
        "05-results/ASMB_MAGS_metaWRAP_postfilter.tsv"

#### Prepare input files for running pipeline ##################################

rule write_sampletsv:
    output:
        "04-analysis/{sampletype}_samples_human-autoFiltering.sample.tsv"
    message: "Prepare sample TSV for the human {wildcards.sampletype} samples for the auto-filtering pipeline"
    params:
        fastadir = "04-analysis/ancient_metagenome_assembly/alignment",
        bamdir = "04-analysis/ancient_metagenome_assembly/alignment"
    run:
        sampletsv = pd.DataFrame.from_dict(SAMPLETYPE[wildcards.sampletype], orient="index", columns=['metawrapreport'])
        sampletsv.index.name = "sample"
        sampletsv['fastafn'] = [f"{os.getcwd()}/{params.fastadir}/{sample.split('-')[1]}/{sample}.fasta.gz"
                                for sample in SAMPLETYPE[wildcards.sampletype]]
        sampletsv['bamfn'] = [f"{os.getcwd()}/{params.bamdir}/{sample.split('-')[1]}/{sample.split('-')[0]}.sorted.dedup.bam"
                                for sample in SAMPLETYPE[wildcards.sampletype]]

        sampletsv.reset_index()[['sample', 'fastafn', 'bamfn', 'metawrapreport']] \
            .to_csv(output[0], sep="\t", index=False)

rule write_config:
    input:
        "04-analysis/{sampletype}_samples_human-autoFiltering.sample.tsv"
    output:
        "04-analysis/{sampletype}_samples_human-autoFiltering.sample.config"
    message: "Prepare config for the human {wildcards.sampletype} samples for the auto-filtering pipeline"
    params:
        config = "/mnt/archgen/users/huebner/automatic_MAG_refinement/config/config.yaml",
        resultdir = "04-analysis/automatic_MAG_refinement/{sampletype}_samples_human",
        tmpdir = "tmp/automatic_MAG_refinement-{sampletype}_samples_human_tmpdir",
    run:
        with open(params.config, "r") as f:
            configfile = yaml.safe_load(f)

        configfile['sampletsv'] = f"{os.getcwd()}/{input[0]}"
        configfile['resultdir'] = f"{os.getcwd()}/{params.resultdir}"
        configfile['tmpdir'] = f"{os.getcwd()}/{params.tmpdir}"

        with open(output[0], "w", encoding="utf8") as outfile:
            yaml.dump(configfile, outfile, default_flow_style=False, allow_unicode=True, sort_keys=False)

################################################################################

#### Run refinment pipeline ####################################################

rule run_refinement_pipeline:
    input:
        tsv = "04-analysis/{sampletype}_samples_human-autoFiltering.sample.tsv",
        config = "04-analysis/{sampletype}_samples_human-autoFiltering.sample.config"
    output:
        "04-analysis/automatic_MAG_refinement/{sampletype}_samples_human/MAG_automaticrefinement_summary.tsv"
    message: "Run the automatic refinement pipeline on {wildcards.sampletype} samples"
    params:
        path_to_pipeline = "../../automatic_MAG_refinement"
    shell:
        """
        cd {params.path_to_pipeline}
        snakemake --configfile {input.config} \
            --use-conda --conda-prefix conda --profile sge -j 10 -k
        """

################################################################################

#### Summary ###################################################################

rule summary:
    input:
        expand("04-analysis/automatic_MAG_refinement/{sampletype}_samples_human/MAG_automaticrefinement_summary.tsv", sampletype=SAMPLETYPE)
    output:
        "05-results/ASMB_MAGS_metaWRAP_postfilter.tsv"
    message: "Summarise the MAG quality post filtering"
    params:
        dir = "04-analysis/automatic_MAG_refinement"
    run:
        # Concatenate the reports
        qual = pd.concat([pd.read_csv(fn, sep='\t')
                          for fn in input])

        # Extract the GUNC columns misisng in the report
        qual['sampletype'] = ["aDNA" if "megahit" in sample else "mDNA"
                              for sample in qual['sample'].tolist()]
        gunc = pd.concat([pd.read_csv(f"{params.dir}/{s.sampletype}_samples_human/{s.sample}/{s.sample}_GUNC_checkM.merged.tsv", sep='\t')
                          for s in qual[['sample', 'sampletype']]
                            .drop_duplicates()
                            .itertuples()])[['genome', 'GUNC.contamination_portion', 'GUNC.n_effective_surplus_clades']] \
                            .rename({'genome': 'sample_binID'}, axis=1)
        qual = qual.merge(gunc, how="left", on="sample_binID")


        # Clear bin and sample IDs
        qual['sample_binID'] = qual['sample_binID'].str.replace(r'-(megahit|metaspades)', '', regex=True)
        qual['sample'] = qual['sample'].str.replace(r'-(megahit|metaspades)', '', regex=True)


        # Select columns
        qual[['sample_binID', 'sample', 'bin',
              'pass.MIMAG_medium', 'pass.MIMAG_high', 'pass.GUNC', 'pass.Polyrate',
              'GUNC.n_contigs', 'GUNC.divergence_level', 'GUNC.contamination_portion',
              'GUNC.n_effective_surplus_clades', 'checkM.genome_size', 'checkM.GC',
              'checkM.coding_density', 'checkM.N50_contigs', 'checkM.lineage',
              'checkM.completeness', 'checkM.contamination', 'checkM.strain_heterogeneity',
              'meanCov', 'breadth', 'polyrate']] \
            .rename({'sample_binID': 'binID',
                     'bin': 'metaWRAP_binID'}, axis=1) \
            .sort_values(['binID']) \
            .to_csv(output[0], sep="\t", index=False)
################################################################################
