################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Dataset S8
#
# Dependent on:
#   - FUNC_pangenome_Climicola.Snakefile
#   - PHYL_Flexilinea_proteinalignment.Snakefile
#   - PHYL_Yersinia_proteinalignment.Snakefile
#
# Alex Huebner, 09/09/22
################################################################################

import os

import pandas as pd

rule all:
    output:
        "06-figures_tables/Dataset_S8.xlsx"
    message: "Generate the XLSX file for the Dataset S8"
    params:
        chlorobiaceae_3kb = "05-results/PHYL_Chlorobiaeceae_fastANI_3kb.tsv",
        chlorobiaceae_1kb = "05-results/PHYL_Chlorobiaeceae_fastANI_1kb.tsv",
        roary = "05-results/FUNC_roary_pangenome.tsv",
        keggko = "05-results/FUNC_KEGG_specificKO.tsv",
        flexilinea_taxa = "05-results/PHYL_Flexilinea_proteintree_taxa.tsv",
        flexilinea_pani = "05-results/PHYL_Flexilinea_fastANI.tsv",
        ypestis_completeness = "05-results/PHYL_Yersinia_completeness.tsv",
        ypestis_pani = "05-results/PHYL_Yersinia_fastANI.tsv"
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

        # Dataset S8a: pairwise ANI for Chlorobium with contigs >= 3 kb
        ## Format sample labels
        chlorobiaceae_3kb = pd.read_csv("05-results/PHYL_Chlorobiaeceae_fastANI_3kb.tsv", sep="\t")
        chlorobiaceae_3kb.columns = [c.replace(".fna", "") if c != "43" else "sample1"
                                     for c in chlorobiaceae_3kb.columns]
        chlorobiaceae_3kb['sample1'] = chlorobiaceae_3kb['sample1'].str.replace(".fna", "", regex=False)
        chlorobiaceae_3kb.to_excel(writer, sheet_name="S8a - pANI Chlorobiaceae 3 kb", index=False,
                                 header=False, startrow=3, float_format="%.4f")
        ## Sheet: pairwise ANI Flexilinea
        s8a_sheet = writer.sheets["S8a - pANI Chlorobiaceae 3 kb"]
        s8a_sheet.write(0, 0, "Table S8a: Overview of the pairwise average "
                        "nucleotide identity (ANI) between all pairs of the "
                        "reconstructed Chlorobiaceae MAGs and published "
                        "comparative genomes of the family Chlorobiaceae. "
                        "Only contigs of a minimum length of 3 kb were considered.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(chlorobiaceae_3kb.columns.values):
            s8a_sheet.write(2, ci, cname, header_format)
        s8a_sheet.set_column(0, 0, determine_col_width(chlorobiaceae_3kb.iloc[:, 0],
                                                       chlorobiaceae_3kb.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in range(1, 43):
            s8a_sheet.set_column(i, i, determine_col_width(chlorobiaceae_3kb.iloc[:, i].astype(str),
                                                           chlorobiaceae_3kb.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0.00"}))
        

        # Dataset S8b: pairwise ANI for Chlorobium with contigs >= 1 kb
        ## Format sample labels
        chlorobiaceae_1kb = pd.read_csv("05-results/PHYL_Chlorobiaeceae_fastANI_1kb.tsv", sep="\t")
        chlorobiaceae_1kb.columns = [c.replace(".fna", "") if c != "Unnamed: 0" else "sample1"
                                     for c in chlorobiaceae_1kb.columns]
        chlorobiaceae_1kb['sample1'] = chlorobiaceae_1kb['sample1'].str.replace(".fna", "", regex=False)
        chlorobiaceae_1kb.to_excel(writer, sheet_name="S8b - pANI Chlorobiaceae 1 kb", index=False,
                                 header=False, startrow=3, float_format="%.4f")
        ## Sheet: pairwise ANI Flexilinea
        s8b_sheet = writer.sheets["S8b - pANI Chlorobiaceae 1 kb"]
        s8b_sheet.write(0, 0, "Table S8b: Overview of the pairwise average "
                        "nucleotide identity (ANI) between all pairs of the "
                        "reconstructed Chlorobiaceae MAGs and published "
                        "comparative genomes of the family Chlorobiaceae. "
                        "Only contigs of a minimum length of 1 kb were considered.",
                         workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(chlorobiaceae_1kb.columns.values):
            s8b_sheet.write(2, ci, cname, header_format)
        s8b_sheet.set_column(0, 0, determine_col_width(chlorobiaceae_1kb.iloc[:, 0],
                                                       chlorobiaceae_1kb.columns[0]) + 1,
                             workbook.add_format({'align': 'center'}))
        for i in range(1, 42):
            s8b_sheet.set_column(i, i, determine_col_width(chlorobiaceae_1kb.iloc[:, i].astype(str),
                                                           chlorobiaceae_1kb.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0.00"}))

        # Dataset S8c: pan-genome analysis with Roary
        roary = pd.read_csv(params.roary, sep="\t")
        roary.columns = [c.replace("-megahit", "") for c in roary.columns]
        roary = roary.iloc[:, list(range(4)) + list(range(14, 35))]
        roary.iloc[:, 4:14] = roary.iloc[:, 4:14].fillna(0).astype(int).astype(str)
        roary.iloc[:, 4:14] = roary.iloc[:, 4:14].apply(lambda c: c.str.replace(r'^0$', '', regex=True), axis=1)
        roary.to_excel(writer, sheet_name="S8c - pan-genome", index=False,
                       header=False, startrow=3, float_format="%.4f")
        ## Sheet: pan-genome analysis
        s8c_sheet = writer.sheets["S8c - pan-genome"]
        s8c_sheet.write(0, 0, "Table S8c: Overview of pan-genome analysis "
                        "resuls obtained from Roary. All annotated genes of "
                        "four published Chlorobium limicola and the six ancient "
                        "Chlorobium MAGs were used as input. We identified a "
                        "gene as specific for either the published genomes or "
                        "the ancient MAGs when it was present in 67% of the "
                        "members of one group and absent in 67% of the other. "
                        "We annotated specific genes using the eggNOG database "
                        "based on the alignments obtained from DIAMOND. The "
                        "gene identifiers for the samples are the ones inferred "
                        "by Prokka.",
                         workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(roary.columns.values):
            s8c_sheet.write(2, ci, cname, header_format)
        for i in list(range(2)) + list(range(4, 25)):
            s8c_sheet.set_column(0, 0, determine_col_width(roary.iloc[:, i],
                                                           roary.columns[i]) + 1,
                                 workbook.add_format({'align': 'center'}))
        s8c_sheet.set_column(3, 3,
                             determine_col_width(roary.iloc[:, 3].astype(str),
                                                 roary.columns[3]) + 1,
                             workbook.add_format({'align': 'center',
                                                 'num_format': "#,##0"}))

        # Dataset S8d: KEGG pathways for clade-specific KEGG kos
        keggko = pd.read_csv(params.keggko, sep="\t")
        keggko = pd.read_csv("05-results/FUNC_KEGG_specificKO.tsv", sep="\t")
        keggko.columns = ["KEGG KO", "pathway 1st level", "pathway 2nd level",
                          "pathway 3rd level", "description", "specificity"]
        keggko['pathway 1st level'] = keggko['pathway 1st level'].str.replace(r"^ +", "", regex=True)
        keggko.to_excel(writer, sheet_name="S8d - KEGG KO and pathways", index=False,
                        header=False, startrow=3, float_format="%.4f")
        ## Sheet: KEGG KO analysis
        s8d_sheet = writer.sheets["S8d - KEGG KO and pathways"]
        s8d_sheet.write(0, 0, "Table S8d: Overview of the KEGG orthologs (KO) "
                        "and their associated pathways for orthologs that were "
                        "found exclusively in either the C. limicola genomes or "
                        "in the ancient Chlorobium MAGs.",
                         workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(keggko.columns.values):
            s8d_sheet.write(2, ci, cname, header_format)
        for i in range(6):
            s8d_sheet.set_column(i, i, determine_col_width(keggko.iloc[:, i],
                                                           keggko.columns[i]) + 1,
                                 workbook.add_format({'align': 'center'}))

        # Dataset S8e: pairwise ANI for Flexilinea with contigs >= 3 kb
        ## Prepare the taxa information for the reference genomes
        flexilinea_taxa = pd.read_csv(params.flexilinea_taxa, sep="\t")
        flexilinea_taxa = flexilinea_taxa.loc[(flexilinea_taxa['NCBI taxonomy Id'] \
            .isin([1946754, 1946750, 2699749, 2017390,
                   1678840, 1946749, 1889813, ])) |
            (flexilinea_taxa['accession Id'] == "GCA_004294845.1")]
        ## Calculate the mean ANI values because the reported ANI values aren't symmetrical
        flexilinea_pani = pd.read_csv(params.flexilinea_pani, sep="\t",
                                      header=None, names=['sample1', 'sample2', 'ANI'], usecols=[0, 1, 2])
        flexilinea_pani['sample1'] = [os.path.basename(fn).replace(".fa", "")
                                      for fn in flexilinea_pani['sample1'].values]
        flexilinea_pani['sample2'] = [os.path.basename(fn).replace(".fa", "")
                                      for fn in flexilinea_pani['sample2'].values]
        flexilinea_pani['pair'] = flexilinea_pani.apply(lambda r: r['sample1'] + "-" + r['sample2']
                                                        if r['sample1'] < r['sample2']
                                                        else r['sample2'] + "-" + r['sample1'],
                                                        axis=1)
        flexilinea_pani = flexilinea_pani.drop(['ANI'], axis=1) \
            .merge(flexilinea_pani.groupby(['pair']).agg({'ANI': 'mean'}),
                   how="left", on="pair") \
            .drop(['pair'], axis=1)
        ## Merge 
        flexilinea_tbl = flexilinea_pani.merge(flexilinea_taxa, how="left",
                                               left_on="sample1", right_on="accession Id")
        flexilinea_tbl[['accession Id', 'name']] = flexilinea_tbl[['accession Id', 'name']].fillna("")
        flexilinea_tbl['NCBI taxonomy Id'] = flexilinea_tbl['NCBI taxonomy Id'].fillna(0).astype(int).astype(str)
        flexilinea_tbl['NCBI taxonomy Id'] = [i if i != "0" else "" for i in flexilinea_tbl['NCBI taxonomy Id'].values]
        flexilinea_wtbl = flexilinea_tbl.drop(['accession Id'], axis=1) \
            .set_index(['sample1', 'name', 'NCBI taxonomy Id', 'sample2']) \
            .unstack(level=-1) \
            .reset_index()
        flexilinea_wtbl.columns = [c[0] if i < 3 else c[1]
                                   for i, c in enumerate(flexilinea_wtbl.columns)]

        flexilinea_wtbl.to_excel(writer, sheet_name="S8e - pANI Flexilinea", index=False,
                                 header=False, startrow=3, float_format="%.4f")
        ## Sheet: pairwise ANI Flexilinea
        s8e_sheet = writer.sheets["S8e - pANI Flexilinea"]
        s8e_sheet.write(0, 0, "Table S8e: Overview of the pairwise average "
                        "nucleotide identity (ANI) between all pairs of the "
                        "reconstructed Flexilinea MAGs and published comparative "
                        "genomes. For these comparative genomes, the common name "
                        "and the NCBI taxonony Id are provided.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(flexilinea_wtbl.columns.values):
            s8e_sheet.write(2, ci, cname, header_format)
        for i in range(0, 3):
            s8e_sheet.set_column(i, i, determine_col_width(flexilinea_wtbl.iloc[:, i],
                                                           flexilinea_wtbl.columns[i]) + 1,
                                 workbook.add_format({'align': 'center'}))
        for i in range(4, 44):
            s8e_sheet.set_column(i, i, determine_col_width(flexilinea_wtbl.iloc[:, i].astype(str),
                                                           flexilinea_wtbl.columns[i]) + 1,
                                 workbook.add_format({'align': 'center',
                                                      'num_format': "0.00"}))

        # Dataset S8f: overview of the MAG quality pre-filtering
        ypestis_completeness = pd.read_csv(params.ypestis_completeness, sep="\t")
        ypestis_hq_samples = ypestis_completeness.query("fracNs < 0.1")['genome'].tolist()
        ## Calculate the mean ANI values because the reported ANI values aren't symmetrical
        ypestis_pani = pd.read_csv(params.ypestis_pani, sep="\t",
                                      header=None, names=['sample1', 'sample2', 'ANI'], usecols=[0, 1, 2])
        ypestis_pani['sample1'] = [os.path.basename(fn).replace(".fa", "")
                                      for fn in ypestis_pani['sample1'].values]
        ypestis_pani['sample2'] = [os.path.basename(fn).replace(".fa", "")
                                      for fn in ypestis_pani['sample2'].values]
        ypestis_pani = ypestis_pani.loc[(ypestis_pani['sample1'].isin(ypestis_hq_samples)) &
                                        (ypestis_pani['sample2'].isin(ypestis_hq_samples))]
        ypestis_pani['pair'] = ypestis_pani.apply(lambda r: r['sample1'] + "-" + r['sample2']
                                                        if r['sample1'] < r['sample2']
                                                        else r['sample2'] + "-" + r['sample1'],
                                                        axis=1)
        ypestis_pani = ypestis_pani.drop(['ANI'], axis=1) \
            .merge(ypestis_pani.groupby(['pair']).agg({'ANI': 'mean'}),
                   how="left", on="pair") \
            .drop(['pair'], axis=1)
        ypestis_pani.to_excel(writer, sheet_name="S8f - pANI Yersinia pestis", index=False,
                              header=False, startrow=3, float_format="%.4f")
        ## Sheet: pairwise ANI Y. pestis
        s8f_sheet = writer.sheets["S8f - pANI Yersinia pestis"]
        s8f_sheet.write(0, 0, "Table S8f: Overview of the pairwise average "
                        "nucleotide identity (ANI) between all pairs of the "
                        "Yersinia pestis genomes that were reported by "
                        "Andrades Valtueña et al. (2022). Genomes were only "
                        "considered when > 90% of the 4.65 Mb were genotyped, "
                        "i.e. had a base other than 'N'.",
                          workbook.add_format({'bold': True, 'align': 'left'}))
        header_format = workbook.add_format({
            'bold': True,
            'align': 'center',
            'valign': 'vcenter',
            'border': 0
        })
        for ci, cname in enumerate(ypestis_pani.columns.values):
            s8f_sheet.write(2, ci, cname, header_format)
        for i in range(0, 2):
            s8f_sheet.set_column(i, i, determine_col_width(ypestis_pani.iloc[:, i],
                                                        ypestis_pani.columns[i]) + 1,
                                 workbook.add_format({'align': 'center'}))
        s8f_sheet.set_column(2, 2,
                                determine_col_width(ypestis_pani.iloc[:, 2].astype(str),
                                                    ypestis_pani.columns[2]) + 1,
                                workbook.add_format({'align': 'center',
                                                    'num_format': "0.00"}))

        # Save XLSX file
        writer.save()
