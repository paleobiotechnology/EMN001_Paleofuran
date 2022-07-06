################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S2
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 04/07/22
################################################################################

import os
from math import ceil

import pandas as pd


rule all:
    output:
        "06-figures_tables/Dataset_S3.xlsx"
    message: "Generate the XLSX file for the Dataset S3"
    params:
        nbins = "05-results/ASMB_nBins_initialbinning.tsv",
        metawrap = "05-results/ASMB_MAGS_metaWRAP_prefilter.tsv",
        mag_preflt = "05-results/ASMB_MAGS_metaWRAP_prefilter.tsv",
        mag_postflt = "05-results/ASMB_MAGS_metaWRAP_postfilter.tsv"
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

        # Dataset S3a: number of bins
        nbins = pd.read_csv(params.nbins, sep="\t")
        nbins.columns = ['sample', 'metaBAT2', 'MaxBin2', 'CONCOCT']
        nrefined = pd.read_csv(params.metawrap, sep="\t") \
            .groupby(['sample'])['pass.MIMAG_medium'].agg(sum).to_frame() \
            .reset_index() \
            .rename({'pass.MIMAG_medium': 'metaWRAP bin refinement'}, axis=1)
        nbins = nbins.merge(nrefined, how="left", on="sample").fillna(0)
        nbins['metaWRAP bin refinement'] = nbins['metaWRAP bin refinement'].astype(int)
        nbins.to_excel(writer, sheet_name="S3a - number of bins", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s3a_sheet = writer.sheets["S3a - number of bins"]
        s3a_sheet.write(0, 0, "Table S3a: Overview of the number of bins "
                        "identified by three non-reference metaBAT2, MaxBin2, "
                        "and CONCOCT and the number of bins retained after "
                        "refining their results with metaWRAP's bin refinement module.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(nbins.columns.values):
            s3a_sheet.write(2, ci, cname, header_format)
        s3a_sheet.set_column(0, 0, determine_col_width(nbins.iloc[:, 0],
                                                       nbins.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in range(1, 5):
            s3a_sheet.set_column(i, i,
                                 determine_col_width(nbins.iloc[:, i].astype(str),
                                                     nbins.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))

        # Dataset S3b: overview of the MAG quality pre-filtering
        mag_preflt = pd.read_csv(params.mag_preflt, sep="\t")
        # Introduce the naming scheme for bins
        mag_preflt['binID'] = [f"{sample}_{int(i):03d}"
                               for sample, i in zip(mag_preflt['sample'],
                                                    mag_preflt['binID'].str.extract(r'bin.([0-9]+)')[0].tolist())]
        # Change format of columns
        mag_preflt['GUNC.contamination_portion'] *= 100
        mag_preflt['GUNC.contamination_portion'] = mag_preflt['GUNC.contamination_portion'].astype(int)
        mag_preflt['checkM.genome_size'] /= 1000000
        # Generate BUSCO
        mag_preflt['BUSCO_generic.C'] = mag_preflt['BUSCO_generic.S'] + mag_preflt['BUSCO_generic.D']
        for i in list(range(25, 29)) + [36]:
            mag_preflt.iloc[:, i] /=  mag_preflt['BUSCO_generic.T']
            mag_preflt.iloc[:, i] *= 100
        mag_preflt['BUSCO_specific.C'] = mag_preflt['BUSCO_specific.S'] + mag_preflt['BUSCO_specific.D']
        for i in list(range(31, 35)) + [37]:
            mag_preflt.iloc[:, i] /=  mag_preflt['BUSCO_specific.T']
            mag_preflt.iloc[:, i] *= 100
        mag_preflt.columns = (list(mag_preflt.columns[:25]) +
                              [f"BUSCO generic {x} [%]" for x in ['duplicated', 'single-copy', 'fragmented', 'missing']] +
                              list(mag_preflt.columns[29:31]) +
                              [f"BUSCO specific {x} [%]" for x in ['duplicated', 'single-copy', 'fragmented', 'missing']] +
                              [mag_preflt.columns[35]] + ['BUSCO generic complete [%]', 'BUSCO specific complete [%]'])
         
        # Select columns and reorder
        dataset_s3b = mag_preflt[['binID', 'checkM.genome_size', 'GUNC.n_contigs',
                                  'checkM.N50', 'checkM.longest_contig', 'checkM.GC',
                                  'checkM.completeness', 'checkM.contamination',
                                  'checkM.strain_heterogeneity',
                                  'GUNC.divergence_level',
                                  'GUNC.contamination_portion',
                                  'GUNC.n_effective_surplus_clades', 'GUNC.CSS',
                                  'BUSCO generic complete [%]', 'BUSCO generic fragmented [%]',
                                  'BUSCO generic missing [%]', 'BUSCO_specific.lineage',
                                  'BUSCO specific complete [%]', 'BUSCO specific fragmented [%]',
                                  'BUSCO specific missing [%]']] \
            .rename({'checkM.genome_size': 'genome size [Mb]',
                     'GUNC.n_contigs': 'no. of contigs',
                     'checkM.N50': 'N50 [bp]',
                     'checkM.longest_contig': 'longest contig [bp]',
                     'checkM.GC': 'GC [%]',
                     'checkM.completeness': 'checkM completeness [%]',
                     'checkM.contamination': 'checkM contamination [%]',
                     'checkM.strain_heterogeneity': 'checkM strain heterogeneity [%]',
                     'GUNC.divergence_level': 'GUNC divergence level',
                     'GUNC.contamination_portion': 'GUNC contamination [%]',
                     'GUNC.n_effective_surplus_clades': 'GUNC effective no. of surplus clades',
                     'GUNC.CSS': 'GUNC clade separation score',
                     'BUSCO_specific.lineage': 'BUSCO specific lineage'}, axis=1) \
            .sort_values(['binID'])

        dataset_s3b.to_excel(writer, sheet_name="S3b - MAG results prefilter", index=False,
                             header=False, startrow=3)
        ## Sheet: MAG overview pre-filtering
        s3b_sheet = writer.sheets["S3b - MAG results prefilter"]
        s3b_sheet.write(0, 0, "Table S3b: Overview of the quality of the metagenome-"
                        "assembled genomes prior to filtering. The quality was "
                        "estimated using checkM, GUNC, and BUSCO. For BUSCO, the "
                        "estimates based on either the set of 124 bacterial genes "
                        "(generic) or for a lineage specific set of genes (specific) "
                        "were reported.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0,
            'text_wrap': True
        })
        for ci, cname in enumerate(dataset_s3b.columns.values):
            s3b_sheet.write(2, ci, cname, header_format)
        for i in [0, 9, 16]:
            s3b_sheet.set_column(i, i, determine_col_width(dataset_s3b.iloc[:, i],
                                                           dataset_s3b.columns[i]) + 1,
                                workbook.add_format({'align': 'center'}))
        for i in list(range(2, 5)) + [10]:  # integer columns
            s3b_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     dataset_s3b.columns[i]),
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0"}))
        s3b_sheet.set_column(5, 5,
                             determine_col_width(None,
                                                 dataset_s3b.columns[5]) + 1,
                             workbook.add_format({'align': 'center',
                                                  'num_format': "0.00"}))
        for i in [1] + list(range(6, 9)) + list(range(11, 16)) + list(range(17, 20)):  # float columns
            s3b_sheet.set_column(i, i,
                                 ceil(determine_col_width(None,
                                                     dataset_s3b.columns[i]) / 2),
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0.00"}))

        # Dataset S3c: overview of the MAG quality post-filtering
        mag_postflt = pd.read_csv(params.mag_postflt, sep="\t")
        mag_postflt = pd.read_csv("05-results/ASMB_MAGS_metaWRAP_postfilter.tsv", sep="\t")

        # Change format of columns
        mag_postflt['GUNC.contamination_portion'] *= 100
        mag_postflt['GUNC.contamination_portion'] = mag_postflt['GUNC.contamination_portion'].astype(int)
        mag_postflt['checkM.genome_size'] /= 1000000
        mag_postflt['polyrate'] *= 100
         
        # Select columns and reorder
        dataset_s3c = mag_postflt[['binID', 'checkM.genome_size', 'GUNC.n_contigs',
                                  'checkM.N50_contigs', 'checkM.GC', 'meanCov',
                                  'checkM.completeness', 'checkM.contamination',
                                  'checkM.strain_heterogeneity',
                                  'GUNC.divergence_level',
                                  'GUNC.contamination_portion',
                                  'GUNC.n_effective_surplus_clades', 'GUNC.CSS', 'polyrate']] \
            .rename({'checkM.genome_size': 'genome size [Mb]',
                     'GUNC.n_contigs': 'no. of contigs',
                     'checkM.N50_contigs': 'N50 [bp]',
                     'checkM.longest_contig': 'longest contig [bp]',
                     'checkM.GC': 'GC [%]',
                     'meanCov': 'mean coverage',
                     'checkM.completeness': 'checkM completeness [%]',
                     'checkM.contamination': 'checkM contamination [%]',
                     'checkM.strain_heterogeneity': 'checkM strain heterogeneity [%]',
                     'GUNC.divergence_level': 'GUNC divergence level',
                     'GUNC.contamination_portion': 'GUNC contamination [%]',
                     'GUNC.n_effective_surplus_clades': 'GUNC effective no. of surplus clades',
                     'GUNC.CSS': 'GUNC clade separation score',
                     'polyrate': 'ratio non-syn. to syn. minor alleles [%]'}, axis=1) \
            .sort_values(['binID'])

        # Evaluate qualities
        MIMAG = {0: 'low', 1: 'medium', 2: 'high'}
        dataset_s3c['MIMAG'] = [MIMAG[i] for i in (mag_postflt['pass.MIMAG_medium'].astype(int) +
                                                   mag_postflt['pass.MIMAG_high'].astype(int))]


        dataset_s3c.to_excel(writer, sheet_name="S3c - MAG results postfilter", index=False,
                             header=False, startrow=3)
        ## Sheet: MAG overview post-filtering
        s3c_sheet = writer.sheets["S3c - MAG results postfilter"]
        s3c_sheet.write(0, 0, "Table S3c: Overview of the quality of the metagenome-"
                        "assembled genomes prior to filtering. The quality was "
                        "estimated using checkM, GUNC, and the script polymut.py from "
                        "the GitHub repository https://github.com/SegataLab/cmseq. This "
                        "script evaluates all sites in coding sequences with minor "
                        "alleles that have a frequency of at least than 20% and "
                        "calculates the ratio between the ones that would cause a "
                        "non-synonymous substitution compared to the major allele "
                        "and the ones that would cause a synonymous substitution.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0,
            'text_wrap': True
        })
        for ci, cname in enumerate(dataset_s3c.columns.values):
            s3c_sheet.write(2, ci, cname, header_format)
        for i in [0, 9, 14]:
            s3c_sheet.set_column(i, i, determine_col_width(dataset_s3c.iloc[:, i],
                                                           dataset_s3c.columns[i]) + 1,
                                workbook.add_format({'align': 'center'}))
        for i in list(range(2, 4)) + [10]:  # integer columns
            s3c_sheet.set_column(i, i,
                                 determine_col_width(None,
                                                     dataset_s3c.columns[i]),
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0"}))
        s3c_sheet.set_column(4, 4,
                             determine_col_width(None,
                                                 dataset_s3c.columns[5]) + 1,
                             workbook.add_format({'align': 'center',
                                                  'num_format': "0.00"}))
        for i in [1] + list(range(5, 9)) + list(range(11, 14)):  # float columns
            s3c_sheet.set_column(i, i,
                                 ceil(determine_col_width(None,
                                                          dataset_s3c.columns[i]) / 2),
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0.00"}))

        # Save XLSX file
        writer.save()
