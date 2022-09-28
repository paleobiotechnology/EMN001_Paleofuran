################################################################################
# Project: Natural products from the Palaeolithic
# Part: De novo assembly of metagenome-assembled genomes
# Step: Concatenate the pyDamage results for the individual contigs into a
#       single table
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 22/07/22
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

rule all:
    output:
        "05-results/ASMB_pyDamage_individualContigs.tsv.gz"
    message: "Concatenate the pyDamage results of the individual contigs into a single file"
    params:
        dir = "04-analysis/ancient_metagenome_assembly/pydamage"
    run:
        pydamage_res = [pd.read_csv(fn, sep=",")
                        .assign(sample=os.path.basename(fn)[:6])
                        for fn in glob(f"{params.dir}/*.csv.gz")
                        if os.stat(fn).st_size > 50]

        pd.concat(pydamage_res) \
            .iloc[:, [51] + list(range(51))] \
            .sort_values(['sample', 'reflen'], ascending =[True, False]) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.4f", compression="gzip")
