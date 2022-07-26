################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S1
#
# Dependent on:
#   - PREP_preprocessing_dentalcalculus_sequencing_data.Snakefile
#   - QUAL_fragmentlengths.Snakefile
#
# Alex Huebner, 08/06/22
################################################################################

import pandas as pd

rule prepare_dataset_s1:
    output:
        "06-figures_tables/Dataset_S1.xlsx"
    message: "Prepare the overview table of the samples and sequencing data"
    params:
        sample = "01-resources/FellowsYates_PNAS2021_TableS1_SheetA.tsv",
        seqdata = "01-resources/overview_sequencingdata.tsv",
        nreads = "05-results/PREP_Nextflow_EAGER_noReads_per_sample.tsv",
        fraglength = "05-results/QUAL_fragmentlength_distribution.tsv"
    run:
        # Read the sequencing data information
        seqdata = pd.read_csv(params.seqdata, sep="\t")
        seqdata.columns = ['individual Id', 'sample Id', 'library Id', 'library type', 'library treatment',
                           'sequencing platform', 'date of library preparation', 'extraction date',
                           'sequencing data', 'sequencing setup', 'ENA secondary sample accession',
                           'ENA run accession', 'publication']

        # Read the sample info from the published table from Fellows Yates et al. (2021)
        sample_info = pd.read_csv(params.sample, sep="\t",
                                  usecols=['Analysis group', 'Common name', 'MPI-SHH Lab Sample ID',
                                           'Archaeological or Museum ID', 'Collection site',
                                           'Collection location', 'Latitude', 'Longitude',
                                           'Calculus collection tooth (FDI)', 'Specimen period',
                                           'Specimen date (A-associated, D-direct; Method, cal. BP: 2Sigma)',
                                           'Museum/Institution']) \
            .rename({'MPI-SHH Lab Sample ID': 'sampleId'}, axis=1)
        sample_info = sample_info.loc[sample_info['sampleId'].isin(seqdata['sample Id'].tolist())]
        sample_info = sample_info.iloc[:, [2, 1, 0] + list(range(3, 12))]
        sample_info['Latitude'] = sample_info['Latitude'].astype(float)
        sample_info['Longitude'] = sample_info['Longitude'].astype(float)

        # Number of DNA molecules per sample and library type
        nreads = pd.read_csv(params.nreads, sep="\t") \
            .drop(['R2'], axis=1) \
            .rename({'sample': 'ID',
                     'R1': 'PE',
                     'R0': 'SE'}, axis=1)
        nreads[['individualId', 'libraryType']] = nreads['ID'].str.split("-", 1, expand=True)
        nreads['libraryType'] = ['all data' if lt == "alldata" else "non-UDG" for lt in nreads['libraryType']]
        nreads = pd.concat([pd.pivot_table(nreads[['individualId', 'libraryType', 'PE']], index="individualId",
                            columns='libraryType', values="PE", fill_value=0, sort=False)
                            .rename({'all data': 'all data - PE',
                                     'non-UDG': 'non-UDG - PE'}, axis=1)
                            .reset_index(),
                            pd.pivot_table(nreads[['individualId', 'libraryType', 'SE']], index="individualId",
                            columns='libraryType', values="SE", fill_value=0, sort=False)
                            .rename({'all data': 'all data - SE',
                                     'non-UDG': 'non-UDG - SE'}, axis=1)
                            .reset_index()
                            .drop(['individualId'], axis=1)], axis=1)
        nreads = nreads.iloc[:, [0, 1, 3, 2, 4]]

        # Fragment length distribution
        fraglength = pd.read_csv(params.fraglength, sep="\t") \
            .drop(['no. of DNA molecules'], axis=1)

        # Export to Excel as multi-sheet document
        writer = pd.ExcelWriter(output[0], engine='xlsxwriter')
        sample_info.to_excel(writer, sheet_name="S1a - Sample overview", index=False,
                             header=False, startrow=3)
        seqdata.to_excel(writer, sheet_name="S1b - Seq. data overview", index=False,
                         header=False, startrow=3)
        nreads.to_excel(writer, sheet_name="S1c - Number of DNA mol.", index=False,
                        header=False, startrow=4)
        fraglength.to_excel(writer, sheet_name="S1d - DNA molecule length dist.",
                            index=False, header=False, startrow=3)

        # Add formatting
        def determine_col_width(col, colname):
            text_width = col.str.len().max()
            col_width = len(colname)
            if text_width >= col_width:
                return text_width
            else:
                return col_width

        workbook = writer.book
        ## Sheet: Sample overview
        samplesheet = writer.sheets['S1a - Sample overview']
        samplesheet.write(0, 0, "Table S1a: Overview of the dental calculus samples. "
                          "The analysis group refers to grouping used in Fellows Yates et al. (2021).",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(sample_info.columns.values):
            samplesheet.write(2, ci, cname, header_format)
        for i in range(sample_info.shape[1]):
            if i not in [6, 7]:
                samplesheet.set_column(i, i,
                                       determine_col_width(sample_info.iloc[:, i],
                                                           sample_info.columns[i]) + 1,
                                       workbook.add_format({'align': 'center'}))
            else:
                samplesheet.set_column(i, i,
                                       None,
                                       workbook.add_format({'align': 'center',
                                                            'num_format': "0.0000"}))
        ## Sheet: Number of reads
        seqdatasheet = writer.sheets["S1b - Seq. data overview"]
        seqdatasheet.write(0, 0, "Table S1b: Overview of the sequencing data used in this study. "
                           "The sequencing setup indicates whether the library was sequenced "
                           "paired-end (PE) or single-end (SE) and the following number indicates "
                           "the nominal sequencing length in bp. "
                           "The column publication indicates whether the sequencing data was "
                           "generated specifically for this study or has been generated previously.",
                           workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(seqdata.columns.values):
            seqdatasheet.write(2, ci, cname, header_format)
        for i in range(seqdata.shape[1]):
            seqdatasheet.set_column(i, i,
                                    determine_col_width(seqdata.iloc[:, i],
                                                        seqdata.columns[i]) + 1,
                                    workbook.add_format({'align': 'center'}))
        ## Sheet: Sequencing data overview
        nreadssheet = writer.sheets["S1c - Number of DNA mol."]
        nreadssheet.write(0, 0, "Table S1c: Overview of the number of DNA molecules per "
                          "individual that were used for assembly. All data were used for "
                          "assembly process, while the non-UDG data were used for the "
                          "verification of the presence of ancient DNA damage. "
                          "PE and SE list the number of DNA molecules that were paired-end "
                          "and single-end sequenced, respectively.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        nreadssheet.merge_range('B3:C3', 'all data', header_format)
        nreadssheet.merge_range('D3:E3', 'non-UDG data', header_format)
        nreadssheet.write(3, 0, "individual Id", header_format)
        nreadssheet.write(3, 1, "PE", header_format)
        nreadssheet.write(3, 2, "SE", header_format)
        nreadssheet.write(3, 3, "PE", header_format)
        nreadssheet.write(3, 4, "SE", header_format)
        nreadssheet.set_column(0, 0,
                               determine_col_width(nreads.iloc[:, 0],
                                                   nreads.columns[0]) + 1,
                               workbook.add_format({'align': 'center'}))
        nreadssheet.set_column(1, 4,
                               12,
                               workbook.add_format({'align': 'center',
                                                    'num_format': "#,##0"}))

        ## Sheet: fragment length distribution
        fraglensheet = writer.sheets["S1d - DNA molecule length dist."]
        fraglensheet.write(0, 0, "Table S1d: Overview of the distribution of "
                           "DNA molecule fragmenth length across samples. The "
                           "molecule length was inferred by overlapping read pairs "
                           "of sequencing data generated on 2x 75 bp paired-end "
                           "Illumina sequencing runs using fastp requiring an "
                           "overlap of at least 11 bp. The proportions of DNA "
                           "molecules with a certain length in basepair are listed. "
                           "Read pairs that could not be merged because their DNA "
                           "molecule was > 140 bp were summarised.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(fraglength.columns.values):
            fraglensheet.write(2, ci, cname, header_format)
        fraglensheet.set_column(0, 0,
                                determine_col_width(fraglength.iloc[:, 0],
                                                    fraglength.columns[0]) + 1,
                                workbook.add_format({'align': 'center'}))
        fraglensheet.set_column(1, 1,
                                17,
                                workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))
        fraglensheet.set_column(2, fraglength.shape[1],
                                6,
                                workbook.add_format({'align': 'center',
                                                     'num_format': "0.0000"}))

        # Save XLSX file
        writer.save()
