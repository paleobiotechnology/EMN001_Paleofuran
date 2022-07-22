################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: Summmarise the pyDamage estimates for the presence of aDNA damage per
#       MAG
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#   - ASMB_automaticRefinement.Snakefile
#
# Alex Huebner, 22/07/22
################################################################################

import os

import pandas as pd
import pyfastx

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
MAGS_OVERVIEWS = pd.read_csv("05-results/ASMB_MAGS_metaWRAP_postfilter.tsv", "\t")
MAGS = [i for i in MAGS_OVERVIEWS['binID'].tolist()
        if not i.startswith("JAE") and not i.startswith("VLC")]
SAMPLES = [s for s in MAGS_OVERVIEWS['sample'].unique().tolist()
           if not s.startswith("JAE") and not s.startswith("VLC")]
################################################################################

#### Auxilliary functions ######################################################
def return_fasta_fn(binID):
    sample = f"{binID.split('_')[0]}-megahit"
    return f"04-analysis/automatic_MAG_refinement/aDNA_samples_human/{sample}/bins/{sample}_{binID.split('_')[1]}.fasta.gz"
################################################################################

rule all:
    output:
        "05-results/ASMB_pyDamage_MAGs.tsv"
    message: "Summarise the pyDamage results on MAG-level"
    params: 
        pydamage_dir = "04-analysis/ancient_metagenome_assembly/pydamage"
    run:
        # Load pyDamage results of all samples
        pydamage = pd.concat([pd.read_csv(f"{params.pydamage_dir}/{sample}-megahit.pydamage.csv.gz", sep=",")
                              for sample in SAMPLES])

        # Calculate the median qvalue, pred. accuracy, anc CtoT freq for the
        # first 5 positions
        pydamage_summary = []
        for binID in MAGS:
            contigs = [name for name, _ in pyfastx.Fasta(return_fasta_fn(binID), build_index=False)]
            pyres_bin = pydamage.loc[pydamage['reference'].isin(contigs)]
            pydamage_summary.append(
                pyres_bin[['qvalue', 'predicted_accuracy', 'CtoT-0', 'CtoT-1',
                           'CtoT-2', 'CtoT-3', 'CtoT-4']].median() \
                .to_frame() \
                .transpose() \
                .assign(binID=binID))

        pd.concat(pydamage_summary) \
            .iloc[:, [7] + list(range(7))] \
            .to_csv(output[0], sep="\t", index=False, float_format="%.4f")
