################################################################################
# Project: Natural products from the Palaeolithic
# Part: Quality analyses
# Step: Summarise the statistics that pyDamage inferred from the de novo
#       assembled contigs with respect to the presence of ancient DNA damage
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 19/06/22
################################################################################

import os

import numpy as np
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("04-analysis/ancient_metagenome_assembly/pydamage/{sample}-megahit.pydamage.csv.gz")
SAMPLES = [s for s in SAMPLES if s not in ['DLV001', 'DLV002', 'SPM001', 'SPM002', 'TAF016']]
################################################################################

rule all:
    input:
        "05-results/QUAL_pyDamage_summary_qvalue_predaccuracy.tsv"

rule summarise_pydamage:
    output:
        temp("tmp/pydamage/{sample}.pydamage.tsv")
    message: "Summarise the pyDamage results: {wildcards.sample}"
    params:
        pydamage = "04-analysis/ancient_metagenome_assembly/pydamage/{sample}-megahit.pydamage.csv.gz"
    threads: 1
    run:
        pydamage = pd.read_csv(params.pydamage, sep=",")

        total_contigs = pydamage.shape[0]
        pydamage = pydamage \
            .query("coverage >= 5")

        # Bin the prediction accuracy in 1% bins
        df = pd.DataFrame.from_dict({"bin": range(21),
                                     "count": np.bincount(np.digitize(pydamage['predicted_accuracy'].values,
                                                                      np.arange(0, 1, 0.05)))}) \
            .set_index(['bin'])
        df['freq'] = df['count'] / pydamage.shape[0]
        df = df.drop(['count'], axis=1).transpose() \
            .assign(sample=wildcards.sample,
                    totalContigs=total_contigs,
                    evalContigs=pydamage.shape[0],
                    ancContigs=(pydamage.qvalue < 0.05).sum())

        df.iloc[:, list(range(21, 25)) + list(range(21))] \
            .to_csv(output[0], sep="\t", index=False, float_format="%.5f")

rule concat_summaries:
    input:
        expand("tmp/pydamage/{sample}.pydamage.tsv", sample=SAMPLES)
    output:
        "05-results/QUAL_pyDamage_summary_qvalue_predaccuracy.tsv"
    message: "Summarise the number of evaluated contigs (cov >= 5), the number of contigs with a q-value <= 0.05, and the predicted accuracy"
    run:
        pd.concat([pd.read_csv(fn, sep="\t")
                   for fn in input]) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.5f")
