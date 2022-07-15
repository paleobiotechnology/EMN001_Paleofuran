################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S4
#
# Dependent on:
#   - PREP_preprocessing_publishedNeanderthalcalculus.Snakefile
#   - PREP_preprocessing_ElMiron_sediments.Snakefile
#
# Alex Huebner, 22/06/22
################################################################################

import pandas as pd

rule prepare_dataset_s4:
    output:
        "06-figures_tables/Dataset_S4.xlsx"
    message: "Prepare the overview table of the published Neanderthal calculus samples, the reconstruction of Chlorobiaceae genomes, and the contamination analysis"
    params:
        nea_samples = "05-results/PREP_Nextflow_EAGER_noReads_Weyrich2017_Neanderthals.tsv",
        ref_alignment_calc = "05-results/QUAL_dentalcalculus_Chlorobiaceae_refalignment.tsv",
        snpad = "05-results/REFG_Chlorobiaceae_genomes_snpAD.tsv"
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

        # Dataset S4a: overview of the sequencing data from Weyrich et al. (2017)
        nea_samples = pd.read_csv(params.nea_samples, sep="\t") \
            .rename({'secondary_sample_accession': 'ENA secondary sample accession',
                     'run_accession': 'ENA run accession'}, axis=1)
        nea_samples.to_excel(writer, sheet_name="S4a - seq. data Weyrich2017", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s4a_sheet = writer.sheets["S4a - seq. data Weyrich2017"]
        s4a_sheet.write(0, 0, "Table S4a: Overview of the available sequencing data "
                        "of four Neanderthal dental calculus samples previously "
                        "published by Weyrich et al. (2017). The data were downloaded "
                        "from ENA processed using nf-core/eager.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(nea_samples.columns.values):
            s4a_sheet.write(2, ci, cname, header_format)
        for i in range(3):  # text columns
            s4a_sheet.set_column(i, i, determine_col_width(nea_samples.iloc[:, i],
                                                           nea_samples.columns[i]) + 1,
                                workbook.add_format({'align': 'center'}))
        s4a_sheet.set_column(3, 4,  # numerical columns
                                determine_col_width(nea_samples.iloc[:, 3].astype(str),
                                                    nea_samples.columns[4]) + 2,
                                workbook.add_format({'align': 'center',
                                                    'num_format': "#,##0"}))

        # Dataset S4b: overview of fraction of reads aligned to the contigs of
        # the Chlorobiaceae MAG of EMN001
        ref_aln_calc = pd.read_csv(params.ref_alignment_calc, sep="\t") \
            .sort_values(['alignedReads'], ascending=False) \
            .rename({'alignedReads': '# of reads against the Chlorobium MAG',
                     'totalReads': 'total # of reads'}, axis=1)
        ref_aln_calc['% aligned'] = (ref_aln_calc['# of reads against the Chlorobium MAG'] * 100 /
                                     ref_aln_calc['total # of reads']).round(2)
        ref_aln_calc.to_excel(writer, sheet_name="S4b - ref. alignment ag. EMN001", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s4b_sheet = writer.sheets["S4b - ref. alignment ag. EMN001"]
        s4b_sheet.write(0, 0, "Table S4b: Overview of the number of reads that "
                        "could be aligned to the contigs of the Chlorobium MAG "
                        "of EMN001 for all dental calculus samples. The available "
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
            s4b_sheet.write(2, ci, cname, header_format)
        s4b_sheet.set_column(0, 0, determine_col_width(ref_aln_calc.iloc[:, 0],
                                                       ref_aln_calc.columns[0]) + 1,
                            workbook.add_format({'align': 'center'}))
        s4b_sheet.set_column(1, 2,  # numerical columns
                             determine_col_width(ref_aln_calc.iloc[:, 2].astype(str),
                                                 ref_aln_calc.columns[2]),
                             workbook.add_format({'align': 'center',
                                                 'num_format': "#,##0"}))
        s4b_sheet.set_column(3, 3,  # float columns
                             determine_col_width(ref_aln_calc.iloc[:, 3].astype(str),
                                                 ref_aln_calc.columns[3]),
                             workbook.add_format({'align': 'center',
                                                 'num_format': "0.00"}))

        # Dataset S4c: overview of snpAD genotyping results
        snpad = pd.read_csv(params.snpad, sep="\t") \
            .rename({'sites_covered': 'sites with coverage >= 1-fold',
                     'sites_genotyped': 'sites genotyped',
                     'HOMR': 'homozygous REF',
                     'HET': 'heterozygous',
                     'HOMA': 'homozygous ALT'}, axis=1)
        snpad = snpad.iloc[:, list(range(6)) + [11, 12, 7, 16, 6, 8, 9, 10, 13, 14, 15, 17]]
        snpad.to_excel(writer, sheet_name="S4c - genotyping snpAD", index=False,
                             header=False, startrow=3)
        ## Sheet: Sample overview
        s4c_sheet = writer.sheets["S4c - genotyping snpAD"]
        s4c_sheet.write(0, 0, "Table S4c: Overview of the snpAD genotyping results "
                        "for the samples with strong evidence of Chlorobium DNA. "
                        "The number of sites to which at least a single read could "
                        "be aligned and the sites genotypes are relative to the size "
                        "of the EMN001 Chlorobium MAG (1.88 Mb). Based on the "
                        "number of genotyped sites, the fraction of sites with a "
                        "homozygous genotype as the reference, with a heterozygous "
                        "genotype, and a homozygous alternative genotype are reported. "
                        "The remainder of the column lists the type of substitutions "
                        "that were observed at sites with a homozygous alternative genotype.",
                        workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(snpad.columns.values):
            s4c_sheet.write(2, ci, cname, header_format)
        s4c_sheet.set_column(0, 0, determine_col_width(snpad.iloc[:, 0],
                                                       snpad.columns[0]),
                            workbook.add_format({'align': 'center'}))
        for i in range(1, 6):  # float columns
            s4c_sheet.set_column(i, i, determine_col_width(snpad.iloc[:, i].astype(str),
                                                           snpad.columns[i]),
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0.00000"}))
        s4c_sheet.set_column(6, 17,  # numerical columns
                                determine_col_width(snpad.iloc[:, i].astype(str),
                                                    snpad.columns[i]) + 1,
                                workbook.add_format({'align': 'center',
                                                     'num_format': "#,##0"}))
        # Save XLSX file
        writer.save()
