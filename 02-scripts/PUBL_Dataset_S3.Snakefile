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
        mag_postflt = "05-results/ASMB_MAGS_metaWRAP_postfilter.tsv",
        rnas = "05-results/ASMB_rRNA_tRNA_presence.tsv",
        repr_mags = "05-results/ASMB_representativeMAGs.tsv",
        taxprof = "05-results/ASMB_taxonomic_profiling_MAGs.tsv",
        oraltaxon = "05-results/ASMB_oraltaxon_classification.tsv"
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

        # Dataset S3d: the number of tRNAs and rRNAs
        rnas = pd.read_csv(params.rnas, sep="\t")
        rnas = rnas.fillna("")
        rnas.to_excel(writer, sheet_name="S3d - tRNAs and rRNAs", index=False,
                      header=False, startrow=3)
        ## Sheet: Sample overview
        s3d_sheet = writer.sheets["S3d - tRNAs and rRNAs"]
        s3d_sheet.write(0, 0, "Table S3d: Overview of the number of tRNAs and "
                        "rRNAs that were identified by tRNAscan-SE and INFERNAL, "
                        "respectively, for all MAGs that fulfilled the MIMAG "
                        "criteria for a high-quality MAG with respect to "
                        "completeness and contamination.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(rnas.columns.values):
            s3d_sheet.write(2, ci, cname, header_format)
        s3d_sheet.set_column(0, 0, determine_col_width(rnas.iloc[:, 0],
                                                       rnas.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in range(1, 26):
            s3d_sheet.set_column(i, i,
                                 determine_col_width(rnas.iloc[:, i].astype(str),
                                                     rnas.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))
        s3d_sheet.set_column(26, 26, determine_col_width(rnas.iloc[:, 26],
                                                         rnas.columns[26]) + 1,
                             workbook.add_format({'align': 'center'}))

        # Dataset S3e: representative MAGs with taxonomic assignments
        ## Load data
        repr_mags = pd.read_csv(params.repr_mags, sep="\t")
        taxprof = pd.read_csv(params.taxprof, sep="\t")
        oraltaxon = pd.read_csv(params.oraltaxon, sep="\t")
        ## Adjust columns
        repr_mags['primary cluster'] = repr_mags['cluster'].str.split("_").str[0].astype(int)
        repr_mags['secondary cluster'] = repr_mags['cluster'].str.split("_").str[1].astype(int)

        dataset_s3e = repr_mags[['binID', 'primary cluster', 'secondary cluster',
                                 'genome size [Mb]', 'N50', 'completeness [%]',
                                 'contamination [%]', 'cluster size', 'members of cluster']] \
            .merge(taxprof, how="left", on="binID") \
            .merge(oraltaxon.iloc[:, [0, 1, 3, 2]], how="left", on="binID")
        dataset_s3e.to_excel(writer, sheet_name="S3e - representative MAGs", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s3e_sheet = writer.sheets["S3e - representative MAGs"]
        s3e_sheet.write(0, 0, "Table S3e: Overview of the representative MAGs obtained "
                        "from all dental calculus samples. The primary and secondary "
                        "cluster indicate whether the genomes had a similarity of >= 90% "
                        "or >= 95%, respectively. The completeness and contamination "
                        "estimates were calculated using checkM. The GTDB and "
                        "single-genome bin (SGB) classification shows the taxonomic "
                        "assignment using either GTDBTK or PhyloPhlAn3's metagenomic "
                        "module. The SGB classification-based metrics consider whether "
                        "this MAG has a known (kSGB) or unknown (uSGB) closely related "
                        "genome. The column 'oral taxon' indicates whether there was a "
                        "reference genome present in the HOMD that shared an ANI of >= "
                        "80% (primary) or >= 95% (secondary).",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(dataset_s3e.columns.values):
            s3e_sheet.write(2, ci, cname, header_format)
        for i in [0, 11, 12, 15, 16, 17, 18]:  # string columns
            s3e_sheet.set_column(i, i,
                                 determine_col_width(dataset_s3e.iloc[:, i],
                                                     dataset_s3e.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                      'valign': 'vcenter'}))
        s3e_sheet.set_column(8, 8,
                             determine_col_width(None,
                                                 dataset_s3e.columns[i]) + 24,
                             workbook.add_format({'align': 'left',
                                                  'valign': 'vcenter',
                                                  'text_wrap': True}))
        for i in [9, 10]:  # string columns with left alignment
            s3e_sheet.set_column(i, i,
                                 determine_col_width(dataset_s3e.iloc[:, i],
                                                     dataset_s3e.columns[i]) + 1,
                                 workbook.add_format({'align': 'left',
                                                      'valign': 'vcenter'}))
        for i in [1, 2, 4, 7]:  # integers
            s3e_sheet.set_column(i, i,
                                 determine_col_width(dataset_s3e.iloc[:, i].astype(str),
                                                     dataset_s3e.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'valign': 'vcenter',
                                                     'num_format': "#,##0"}))
        for i in [3, 5, 6, 13, 14, 19]:  # integers
            s3e_sheet.set_column(i, i,
                                 determine_col_width(dataset_s3e.iloc[:, i].astype(str),
                                                     dataset_s3e.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                      'valign': 'vcenter',
                                                      'num_format': "0.00"}))

        # Dataset S3f: number of modern and ancient MAGs per representative MAG 
        all_mags = []
        for mag in repr_mags[['binID', 'members of cluster']].itertuples():
            if type(mag._2) == float:
                all_mags.append((mag.binID, mag.binID))
            else:
                for redmag in mag._2.split(", "):
                    all_mags.append((mag.binID, redmag))
                all_mags.append((mag.binID, mag.binID))
        all_mags_df = pd.DataFrame(all_mags, columns=['binID', 'redundant_MAG']) \
            .merge(mag_postflt[['binID', 'pass.MIMAG_high']], how="left",
                   left_on="redundant_MAG", right_on="binID") \
            .rename({'binID_x': 'representative MAG',
                     'redundant_MAG': 'redundant MAG',
                     'pass.MIMAG_high': 'MIMAG'}, axis=1) \
            .drop(['binID_y'], axis=1)
        all_mags_df['sampletype'] = ((all_mags_df['redundant MAG'].str.startswith("JAE")) |
                                     (all_mags_df['redundant MAG'].str.startswith("VLC")))
        nonred_mag_clustercount = all_mags_df.groupby(['representative MAG', 'sampletype'])['MIMAG'].value_counts() \
            .to_frame() \
            .rename({'MIMAG': 'n'}, axis=1) \
            .reset_index()
        nonred_mag_clustercount['sampletype'] = ["modern" if s else "ancient"
                                                 for s in nonred_mag_clustercount['sampletype'].tolist()]
        nonred_mag_clustercount['MIMAG'] = ["HQ" if m else "MQ"
                                            for m in nonred_mag_clustercount['MIMAG'].tolist()]
        nonred_mag_clustercount['columnname'] = nonred_mag_clustercount['sampletype'] + " - " + nonred_mag_clustercount['MIMAG']
        nonred_mag_clustercount = nonred_mag_clustercount[['representative MAG', 'columnname', 'n']] \
            .set_index(['representative MAG', 'columnname']) \
            .unstack() \
            .fillna(0) \
            .astype(int) \
            .reset_index()
        nonred_mag_clustercount.columns = [c[0] if i == 0 else c[1]
                                           for i, c in enumerate(nonred_mag_clustercount.columns)]
        nonred_mag_clustercount.to_excel(writer, sheet_name="S3f - overview of MAG clusters", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s3f_sheet = writer.sheets["S3f - overview of MAG clusters"]
        s3f_sheet.write(0, 0, "Table S3f: Overview of the number of MAGs "
                        "that were clustered with the representative MAGs "
                        "with respect to their sample type (ancient or sample) "
                        "and their quality following MIMAG criteria.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(nonred_mag_clustercount.columns.values):
            s3f_sheet.write(2, ci, cname, header_format)
        s3f_sheet.set_column(0, 0, determine_col_width(nonred_mag_clustercount.iloc[:, 0],
                                                       nonred_mag_clustercount.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in range(1, 5):
            s3f_sheet.set_column(i, i,
                                 determine_col_width(nonred_mag_clustercount.iloc[:, i].astype(str),
                                                     nonred_mag_clustercount.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))

        # Dataset S3g: aDNA damage against 
        ## Load data
        # repr_mags = pd.read_csv(params.repr_mags, sep="\t")
        # taxprof = pd.read_csv(params.taxprof, sep="\t")
        # oraltaxon = pd.read_csv(params.oraltaxon, sep="\t")
        # ## Adjust columns
        # repr_mags['primary cluster'] = repr_mags['cluster'].str.split("_").str[0].astype(int)
        # repr_mags['secondary cluster'] = repr_mags['cluster'].str.split("_").str[1].astype(int)

        # dataset_s3e = repr_mags[['binID', 'primary cluster', 'secondary cluster',
                                 # 'genome size [Mb]', 'N50', 'completeness [%]',
                                 # 'contamination [%]', 'cluster size', 'members of cluster']] \
            # .merge(taxprof, how="left", on="binID") \
            # .merge(oraltaxon.iloc[:, [0, 1, 3, 2]], how="left", on="binID")
        # dataset_s3e.to_excel(writer, sheet_name="S3e - representative MAGs", index=False,
                             # header=False, startrow=3)
        # ## Sheet: Sample overview
        # s3e_sheet = writer.sheets["S3e - representative MAGs"]
        # s3e_sheet.write(0, 0, "Table S3e: Overview of the representative MAGs obtained "
                        # "from all dental calculus samples. The primary and secondary "
                        # "cluster indicate whether the genomes had a similarity of >= 90% "
                        # "or >= 95%, respectively. The completeness and contamination "
                        # "estimates were calculated using checkM. The GTDB and "
                        # "single-genome bin (SGB) classification shows the taxonomic "
                        # "assignment using either GTDBTK or PhyloPhlAn3's metagenomic "
                        # "module. The SGB classification-based metrics consider whether "
                        # "this MAG has a known (kSGB) or unknown (uSGB) closely related "
                        # "genome. The column 'oral taxon' indicates whether there was a "
                        # "reference genome present in the HOMD that shared an ANI of >= "
                        # "80% (primary) or >= 95% (secondary).",
                        # workbook.add_format({'bold': True, 'align': 'left'}))
        # header_format = workbook.add_format({
            # 'bold': True,
            # 'align': 'center',
            # 'valign': 'vcenter',
            # 'border': 0
        # })
        # for ci, cname in enumerate(dataset_s3e.columns.values):
            # s3e_sheet.write(2, ci, cname, header_format)
        # for i in [0, 11, 12, 15, 16, 17, 18]:  # string columns
            # s3e_sheet.set_column(i, i,
                                 # determine_col_width(dataset_s3e.iloc[:, i],
                                                     # dataset_s3e.columns[i]) + 1,
                                 # workbook.add_format({'align': 'center',
                                                      # 'valign': 'vcenter'}))
        # s3e_sheet.set_column(8, 8,
                             # determine_col_width(None,
                                                 # dataset_s3e.columns[i]) + 24,
                             # workbook.add_format({'align': 'left',
                                                  # 'valign': 'vcenter',
                                                  # 'text_wrap': True}))
        # for i in [9, 10]:  # string columns with left alignment
            # s3e_sheet.set_column(i, i,
                                 # determine_col_width(dataset_s3e.iloc[:, i],
                                                     # dataset_s3e.columns[i]) + 1,
                                 # workbook.add_format({'align': 'left',
                                                      # 'valign': 'vcenter'}))
        # for i in [1, 2, 4, 7]:  # integers
            # s3e_sheet.set_column(i, i,
                                 # determine_col_width(dataset_s3e.iloc[:, i].astype(str),
                                                     # dataset_s3e.columns[i]) + 1,
                                 # workbook.add_format({'align': 'center',
                                                     # 'valign': 'vcenter',
                                                     # 'num_format': "#,##0"}))
        # for i in [3, 5, 6, 13, 14, 19]:  # integers
            # s3e_sheet.set_column(i, i,
                                 # determine_col_width(dataset_s3e.iloc[:, i].astype(str),
                                                     # dataset_s3e.columns[i]) + 1,
                                 # workbook.add_format({'align': 'center',
                                                      # 'valign': 'vcenter',
                                                      # 'num_format': "0.00"}))


        # Save XLSX file
        writer.save()
