################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S2
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#   - QUAL_freeBayes_nomismatches_MEGAHIT.Snakefile
#   - QUAL_pyDamage_summary.Snakefile
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
        asmb_stats = "05-results/ASMB_assemblystats_calN50_metaQUAST.tsv",
        fb_nsubst = "05-results/QUAL_freeBayes_nSubsts.tsv",
        fb_mafs = "05-results/QUAL_freeBayes_distributionMAF.tsv",
        fb_ncontigs = "05-results/QUAL_freeBayes_nContigs.tsv",
        pydamage = "05-results/QUAL_pyDamage_summary_qvalue_predaccuracy.tsv"
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

        # Dataset S2a: summary of the assembly stats
        asmb_stats = pd.read_csv(params.asmb_stats, sep="\t")
        asmb_stats = asmb_stats.iloc[:,:19]
        asmb_stats.columns = ['sample', 'total length [bp]', '# of contigs',
                              '# of contigs >= 1,000 bp', '# of contigs >= 5,000 bp',
                              '# of contigs >= 10,000 bp', '# of contigs >= 25,000 bp',
                              '# of contigs >= 50,000 bp', 'N0', 'N10', 'N20', 'N30',
                              'N40', 'N50', 'N60', 'N70', 'N80', 'N90', 'N100']
        asmb_stats.to_excel(writer, sheet_name="S2a - de novo assembly stats", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s2a_sheet = writer.sheets["S2a - de novo assembly stats"]
        s2a_sheet.write(0, 0, "Table S2a: Overview of performance of the de novo "
                        "assembly of the metagenomic sequencing data. Only contigs "
                        "with a minimal length of 500 bp were considered. The columns N0 "
                        "to N100 represent the maximal and minimal contig length plus "
                        "the nine deciles. The number of coding sequences was determined "
                        "prodigal and the number of genes using Prokka.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(asmb_stats.columns.values):
            s2a_sheet.write(2, ci, cname, header_format)
        s2a_sheet.set_column(0, 0, determine_col_width(asmb_stats.iloc[:, 0],
                                                       asmb_stats.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in range(1, 8):
            s2a_sheet.set_column(i, i,
                                 determine_col_width(asmb_stats.iloc[:, i].astype(str),
                                                     asmb_stats.columns[i]),
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))
        for i in range(8, asmb_stats.shape[1]):
            s2a_sheet.set_column(i, i,
                                 determine_col_width(asmb_stats.iloc[:, i].astype(str),
                                                     asmb_stats.columns[i]) + 2,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))

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

        # Dataset S2c: PyDamage results
        pydamage = pd.read_csv("05-results/QUAL_pyDamage_summary_qvalue_predaccuracy.tsv", sep="\t")
        pydamage['ancContigs'] = (pydamage['ancContigs'] * 100 / pydamage['evalContigs']).round(decimals=1)
        pydamage['evalContigs'] = (pydamage['evalContigs'] * 100 / pydamage['totalContigs']).round(decimals=1)
        pydamage.columns = (['sample', 'total # of contigs', '% of contigs with coverage >= 5x',
                             '% of contigs with q-value < 0.05', 'PA == 0%'] + 
                            [f"{i - 5}% <= PA < {i}%" for i in range(5, 105, 5)])
        pydamage = pydamage.iloc[:, list(range(4)) + list(range(8, 25))]
        pydamage.to_excel(writer, sheet_name="S2c - pyDamage analysis", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s2c_sheet = writer.sheets["S2c - pyDamage analysis"]
        s2c_sheet.write(0, 0, "Table S2c: Overview of the evaluation of the "
                              "presence of ancient DNA damage using pyDamage. "
                              "We only considered contigs that had a minimal coverage "
                              "of 5-fold to avoid wrong inference due to the lack of "
                              "sufficient data. From these contigs, we considered "
                              "the contigs 'ancient', when the q-value returned by "
                              "pyDamage was < 0.05. The columns '15% <= PA < 20%' to "
                              "'95% <= PA < 100%' summarise the number of contigs that "
                              "had a predicted accuracy (PA) in this range.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(pydamage.columns.values):
            s2c_sheet.write(2, ci, cname, header_format)
        s2c_sheet.set_column(0, 0, determine_col_width(pydamage.iloc[:, 0],
                                                       pydamage.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        s2c_sheet.set_column(1, 1,
                             determine_col_width(pydamage.iloc[:, 1].astype(str),
                                                 pydamage.columns[i]) + 1,
                             workbook.add_format({'align': 'center',
                                                 'num_format': "0"}))
        for i in [2, 3]:  # float columns with a single digit
            s2c_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     pydamage.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "0.0"}))
        for i in range(4, 21):  # float columns with five digits
            s2c_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     pydamage.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "0.00000"}))

        # Save XLSX file
        writer.save()
