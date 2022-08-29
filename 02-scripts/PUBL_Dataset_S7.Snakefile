################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S7
#
# Dependent on:
#   - PREP_preprocessing_ElMiron_sediments.Snakefile
#   - PREP_preprocessing_ElMiron_toebone.Snakefile
#   - PREP_preprocessing_lab_negcontrols.Snakefile
#   - QUAL_refalignment_Chlorobiaceae_controls.Snakefile
#
# Alex Huebner, 28/08/22
################################################################################

import pandas as pd

rule prepare_dataset_s7:
    output:
        "06-figures_tables/Dataset_S7.xlsx"
    message: "Prepare the overview table of the sediment and control samples plus the alignment against the Chlorobium MAG of EMN001"
    params:
        sediments = "01-resources/overview_sediments.tsv",
        sediments_nreads = "05-results/PREP_Nextflow_EAGER_noReads_ElMiron_sediments.tsv",
        toebone = "01-resources/overview_ElMiron_humanremains.tsv",
        toebone_nreads = "05-results/PREP_Nextflow_EAGER_noReads_ElMiron_toebone.tsv",
        labcontrols = "01-resources/overview_labcontrols.tsv",
        labcontrols_nreads = "05-results/PREP_Nextflow_EAGER_noReads_labcontrols.tsv",
        ref_alignment_calc = "05-results/QUAL_controls_Chlorobiaceae_refalignment.tsv"
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

        # Dataset S7a: overview of the sequencing data from El Miron sediment samples
        sediments = pd.read_csv(params.sediments, sep="\t") \
            .drop(['latitude', 'longitude', 'country', 'locality'], axis=1)
        sediments['description'] = "sediment from " + sediments['description']
        sediments_nreads = pd.read_csv(params.sediments_nreads, sep="\t")
        sediments_nreads['number of reads'] = sediments_nreads['R1'] + sediments_nreads['R2']
        sediments = sediments.merge(sediments_nreads[['sample', 'number of reads']],
                                    how="left", left_on="sampleId", right_on="sample") \
            .drop(['sample'], axis=1)
        toebone = pd.read_csv(params.toebone, sep="\t") \
            .assign(description="The Red Lady of El Miron (toebone)")
        toebone_nreads = pd.read_csv(params.toebone_nreads, sep="\t") \
            .rename({'R0': 'number of reads'}, axis=1)
        toebone = toebone.merge(toebone_nreads,
                                how="left", left_on="sampleId", right_on="sample") \
            .drop(['sample'], axis=1)
        elmiron = pd.concat([sediments, toebone])

        elmiron.to_excel(writer, sheet_name="S7a - El Miron controls", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s7a_sheet = writer.sheets["S7a - El Miron controls"]
        s7a_sheet.write(0, 0, "Table S7a: Overview of the shotgun sequencing "
                        "data of the sediment samples that were collected "
                        "at the archaeological layers of the site of El Miron "
                        "(Cantabria, Spain) and the El Miron toe bone sample."
                        "The data were processed using nf-core/eager.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(sediments.columns.values):
            s7a_sheet.write(2, ci, cname, header_format)
        for i in [0, 1] + list(range(4, 10)):  # text columns
            s7a_sheet.set_column(i, i, determine_col_width(sediments.iloc[:, i],
                                                           sediments.columns[i]) + 1,
                                workbook.add_format({'align': 'center'}))
        s7a_sheet.set_column(10, 10,  # numerical columns
                                determine_col_width(sediments.iloc[:, 10].astype(str),
                                                    sediments.columns[10]) + 2,
                                workbook.add_format({'align': 'center',
                                                    'num_format': "#,##0"}))

        # Dataset S7b: overview of the sequencing data from the lab controls
        labcontrols = pd.read_csv(params.labcontrols, sep="\t", \
                                  dtype={'extraction_batch': str,
                                         'library_batch': str})
        labcontrols_nreads = pd.read_csv(params.labcontrols_nreads, sep="\t")
        labcontrols_nreads['number of reads'] = labcontrols_nreads.drop(['sample'], axis=1).sum(axis=1)
        labcontrols = labcontrols.merge(labcontrols_nreads[['sample', 'number of reads']],
                                        how="left", left_on="sampleId", right_on="sample") \
            .drop(['sample'], axis=1) \
            .sort_values(['sampleId'])
        labcontrols.to_excel(writer, sheet_name="S7b - laboratory controls", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s7b_sheet = writer.sheets["S7b - laboratory controls"]
        s7b_sheet.write(0, 0, "Table S7b: Overview of the shotgun sequencing "
                        "data of the extraction (EXB) and library (LIB) "
                        "negative control samples that were run alongside "
                        "the dental calculus and the El Miron sediment samples. "
                        "The data were processed using nf-core/eager.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(labcontrols.columns.values):
            s7b_sheet.write(2, ci, cname, header_format)
        for i in range(11):  # text columns
            s7b_sheet.set_column(i, i, determine_col_width(labcontrols.iloc[:, i].astype(str),
                                                           labcontrols.columns[i]) + 1,
                                workbook.add_format({'align': 'center'}))
        s7b_sheet.set_column(11, 11,
                             determine_col_width(labcontrols.iloc[:, 11].astype(str),
                                                 labcontrols.columns[11]) + 2,
                             workbook.add_format({'align': 'center',
                                                 'num_format': "#,##0"}))

        # Dataset S7c: overview of fraction of reads aligned to the contigs of
        # the Chlorobiaceae MAG of EMN001
        ref_aln_calc = pd.read_csv(params.ref_alignment_calc, \
                                   sep="\t", index_col=['sample']) \
            .rename({'alignedReads': '# of reads against the Chlorobium MAG',
                     'totalReads': 'total # of reads'}, axis=1)
        ref_aln_calc['% aligned'] = (ref_aln_calc['# of reads against the Chlorobium MAG'] * 100 /
                                     ref_aln_calc['total # of reads']).round(4)
        ref_aln_calc = pd.concat([ref_aln_calc.loc[["ElMiron_toebone"]],
                                  ref_aln_calc.drop(["ElMiron_toebone"])])
        ref_aln_calc.index.values[0] = 'ElMiron'
        ref_aln_calc = ref_aln_calc.reset_index()
        ref_aln_calc.loc[ref_aln_calc['sample'].str.startswith("EMN"), "sample"] += ".A"
        ref_aln_calc.loc[ref_aln_calc['sample'].str.startswith("EXB"), "sample"] = \
            ref_aln_calc.loc[ref_aln_calc['sample'].str.startswith("EXB"), "sample"].str.replace("_", ".", regex=False)
        ref_aln_calc.loc[ref_aln_calc['sample'].str.startswith("LIB"), "sample"] = \
            ref_aln_calc.loc[ref_aln_calc['sample'].str.startswith("LIB"), "sample"].str.replace("_", ".", regex=False)
        ref_aln_calc.to_excel(writer, sheet_name="S7c - Chlorobium in controls", index=False,
                              header=False, startrow=3)
        ## Sheet: Sample overview
        s7c_sheet = writer.sheets["S7c - Chlorobium in controls"]
        s7c_sheet.write(0, 0, "Table S7c: Overview of the number of reads that "
                        "could be aligned to the contigs of the Chlorobium MAG "
                        "of EMN001 for El Miron sediment samples, the El Miron "
                        "human remain, and the laboratory controls. The available "
                        "short-read sequencing data were aligned against all contigs "
                        "assembled from the sample EMN001 to avoid spurious "
                        "alignments in the absence of any other reference genomes.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0,
            'text_wrap': True,
        })
        for ci, cname in enumerate(ref_aln_calc.columns.values):
            s7c_sheet.write(2, ci, cname, header_format)
        s7c_sheet.set_column(0, 0, determine_col_width(ref_aln_calc.iloc[:, 0],
                                                       ref_aln_calc.columns[0]) + 1,
                            workbook.add_format({'align': 'center'}))
        s7c_sheet.set_column(1, 2,  # numerical columns
                             determine_col_width(ref_aln_calc.iloc[:, 2].astype(str),
                                                 ref_aln_calc.columns[2]),
                             workbook.add_format({'align': 'center',
                                                 'num_format': "#,##0"}))
        s7c_sheet.set_column(3, 3,  # float columns
                             determine_col_width(ref_aln_calc.iloc[:, 3].astype(str),
                                                 ref_aln_calc.columns[3]),
                             workbook.add_format({'align': 'center',
                                                 'num_format': "0.00"}))

        # Save XLSX file
        writer.save()
