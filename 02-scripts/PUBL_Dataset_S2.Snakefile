################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S2
#
# Dependent on:
#   - QUAL_freeBayes_nomismatches_MEGAHIT.Snakefile
#
# Alex Huebner, 15/06/22
################################################################################

import os

import pandas as pd


rule all:
    output:
        "06-figures_tables/Dataset_S2.xlsx"
    message: "Generate the XLSX file for the Dataset S2"
    params: 
        fb_nsubst = "05-results/QUAL_freeBayes_nSubsts.tsv",
        fb_mafs = "05-results/QUAL_freeBayes_distributionMAF.tsv",
        fb_ncontigs = "05-results/QUAL_freeBayes_nContigs.tsv"
    threads: 1
    run:
        # Auxilliary functions
        ## Add formatting
        def determine_col_width(col, colname):
            if col is not None:
                text_width = col.str.len().max()
            else:
                text_width = 0
            col_width = len(colname)
            if text_width >= col_width:
                return text_width
            else:
                return col_width

        writer = pd.ExcelWriter(output[0], engine='xlsxwriter')
        workbook = writer.book

        # Dataset S2b: genotype calling using freeBayes vs. MEGAHIT
        ## Number of substitutions
        fb_subst = pd.read_csv(params.fb_nsubst, sep="\t")
        fb_subst.iloc[:,1:] = fb_subst.iloc[:,1:].astype(int)
        fb_subst = fb_subst[['sample', 'AG', 'CT', 'GA', 'TC',
                             'AC', 'AT', 'CA', 'CG',
                             'GC', 'GT', 'TA', 'TG']]
        fb_subst['# of differences'] = fb_subst.iloc[:, 1:].sum(axis=1)
        fb_subst['transitions'] = fb_subst.iloc[:, 1:5].sum(axis=1) / fb_subst['# of differences']
        fb_subst['transversions'] = fb_subst.iloc[:, 5:13].sum(axis=1) / fb_subst['# of differences']
        fb_subst = fb_subst[['sample', '# of differences', 'transitions', 'transversions'] +
                            fb_subst.columns.tolist()[1:13]]

        ## Minor allele frequency distribution
        fb_maf = pd.read_csv(params.fb_mafs, sep="\t")
        fb_maf.iloc[:,1:] = fb_maf.iloc[:,1:].astype(int)

        fb_maf_long = fb_maf.set_index(['sample']) \
            .stack() \
            .reset_index(level=1)
        fb_maf_long['level_1'] = fb_maf_long['level_1'].astype(int) + 1
        fb_maf_long['10p_bins'] = fb_maf_long['level_1'].floordiv(2)
        fb_maf = fb_maf_long.groupby([fb_maf_long.index, '10p_bins'])[0].sum() \
            .reset_index() \
            .query('`10p_bins` >= 4') \
            .set_index(['sample', '10p_bins']) \
            .unstack() \
            .reset_index()
        fb_maf.columns = ['sample', "30% < AF <= 40%",
                          "40% < AF <= 50%",
                          "50% < AF <= 60%",
                          "60% < AF <= 70%",
                          "70% < AF <= 80%",
                          "80% < AF <= 90%",
                          "90% < AF <= 100%"]
        ## Number of affected contigs
        fb_ncontigs = pd.read_csv(params.fb_ncontigs, sep="\t") \
            .rename({'nContigs': 'number of contigs'}, axis=1)
        dataset_s2b = fb_ncontigs.merge(fb_subst, how="left", on="sample") \
            .merge(fb_maf, how="left", on="sample") \
            .rename({'sample': 'individualId'}, axis=1)
        dataset_s2b = dataset_s2b.loc[~dataset_s2b['individualId'].str.startswith("DLV")]

        dataset_s2b.to_excel(writer, sheet_name="S2b - consensus correction", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s2b_sheet = writer.sheets["S2b - consensus correction"]
        s2b_sheet.write(0, 0, "Table S2b: Overview of the number of differences "
                              "between the contig sequences that were exported by "
                              "MEGAHIT compared to the genotypes obtained from "
                              "freeBayes after aligning the short-read sequencing "
                              "data back against the assembled contigs. The columns "
                              "'AG' to 'TG' indicate how often each substitution was "
                              "observed. The columns 'frequency < AF <= frequency' "
                              "show the distrubtion of the frequency of the allele "
                              "(AF) of the freeBayes genotype binned into 10% bins.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(dataset_s2b.columns.values):
            s2b_sheet.write(2, ci, cname, header_format)
        s2b_sheet.set_column(0, 0, determine_col_width(dataset_s2b.iloc[:, 0],
                                                       dataset_s2b.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in [1, 2] + list(range(17, 24)):  # integer columns
            s2b_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     dataset_s2b.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "0"}))
        for i in range(5, 17):  # integer columns
            s2b_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     dataset_s2b.columns[i]) + 4,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "0"}))
        for i in [3, 4]:  # float columns
            s2b_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     dataset_s2b.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "0.00"}))

        # Save XLSX file
        writer.save()
