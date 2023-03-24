## Overview of the folder `05-results`

This folder contains the results of the analyses that are used to produce the figures and tables
presented in the manuscript.

0. Folders:

  - `genomes`: contains the FastA files of the metagenome-assembled genomes used for the
    phylogenetic analyses

1. `PREP`: results related to the data preparation

  - `PREP_Nextflow_EAGER_noReads_per_sample.tsv`: the overview of the number of reads per individual
    and library type (all data and non-UDG data only) that were available for *de novo* assembly
  - `PREP_Nextflow_EAGER_noReads_Weyrich2017_Neanderthals.tsv`: the overview of the number of reads
    available for the four Neanderthal calculus samples published by Weyrich *et al.* (2017)
  - `PREP_Nextflow_EAGER_noReads_ElMiron_sediments.tsv`: the overview of the number of reads
    available for the six El Miron sediment samples generated for this study
  - `PREP_Nextflow_EAGER_noReads_ElMiron_toebone.tsv`: the overview of the number of
    reads available for the El Miron toe bone sample, whose metagenomic sequencing data was
    generated for Fu *et al.* (2016), but has not been published
  - `PREP_Nextflow_EAGER_noReads_labcontrols.tsv`: the overview of the number of reads available for
    the extraction and library control samples

2. `ASMB`: results related to the *de novo* assembly and non-reference binning

  - `ASMB_assemblystats_calN50_metaQUAST.tsv`: the overview of the performance of the *de novo*
    assembly summarised using calN50, metaQUAST, Prokka and pyDamage
  - `ASMB_nBins_initialbinning.tsv`: the overview of the number of bins that were returned by the
    individual binners metabat2, maxbin2, and concoct for each sample
  - `ASMB_MAGS_metaWRAP_prefilter.tsv`: the overview of the quality of the MAGs refined by metaWRAP
    prior to filtering for chimeric contigs
  - `ASMB_MAGS_metaWRAP_postfilter.tsv`: the overview of the quality of the MAGs refined by metaWRAP
    after filtering for chimeric contigs
  - `ASMB_rRNA_tRNA_presence.tsv`: the results for the presence of tRNAs and rRNAs in the MAGs that
    fulfill the completeness and contamination estimates for a high-quality MAG following MIMAG
  - `ASMB_representativeMAGs.tsv`: the overview of the representative MAGs obtained by dereplicating
    them using dRep dereplicate
  - `ASMB_oraltaxon_classification.tsv`: the overview of whether or not there was a close genome to
    the representative MAGs in the HOMD database and their respective evolutionary distance
  - `ASMB_taxonomic_profiling_MAGs.tsv`: the overview of the taxonomic assignments of the
    representative MAGs against the GTDB and the SGB database
  - `ASMB_pyDamage_individualContigs.tsv.gz`: the pyDamage results for all individual contigs across
    all ancient dental calculus samples
  - `ASMB_pyDamage_MAGs.tsv`: the summary of the pyDamage results across all contigs of a MAG
    generated from ancient dental calculus samples

3. `QUAL`: results related to the data quality

  - `QUAL_fragmentlength_distribution.tsv`: the overview of the distribution of DNA molecule lengths
    across the sequencing data generated on 2x 75 bp paired-end Illumina sequencing runs
  - `QUAL_freeBayes_nSubsts.tsv`: the overview of the frequency of each allele difference between the
    allele exported by MEGAHIT from the assembly graph and the genotype call inferred by freeBayes
    was observed
  - `QUAL_freeBayes_distributionMAF.tsv`: the distribution of the allele frequencies for the allele
    that freeBayes inferred to be the genotype by 5% frequency bins
  - `QUAL_freeBayes_nContigs.tsv`: the overview of the number of contigs per sample for which a
    difference between the MEGAHIT and freeBayes sequence could be observed
  - `QUAL_pyDamage_summary_qvalue_predaccuracy.tsv`: summary of the pyDamage results regarding the
    fraction of contigs that are considered "ancient" by pyDamage and the distribution of the
    predicted accuracy
  - `QUAL_dentalcalculus_Chlorobiaceae_refalignment.tsv`: the overview of the number of reads
    aligned to the Chlorobiaceae MAG contigs and to all contigs of EMN001 for each dental calculus
    sample
  - `QUAL_damageprofile_Climicola_Abot439.tsv`: the substitution frequency for all five samples with
    HQ Chlorobiaceae MAGs against the genomes of the species *Chlorobium limicola* and Abot439 at
    the 5' end of the reads

4. `REFG`: results related to constructing genomes using reference genomes

  - `REFG_Chlorobiaceae_genomes_snpAD.tsv`: the overview of the performance of the recovery of
    Chlorobiaceae genomes from samples with strong evidence of Chlorobiaceae DNA using the EMN001
    Chlorobiaceae MAG as the reference

5. `PHYL`: results related to phylogenetic analyses

  - `PHYL_Chlorobiaeceae_fastANI_3kb.tsv`: the matrix of the pairwise average nucleotide identity
    values calculated for a subset of the genomes of the family Chlorobiaceae considering only
    contigs with a minimum length of 3 kb
  - `PHYL_Chlorobiaeceae_fastANI_3kb_mapped_fragments.tsv`: The raw fastANI output with percentage of 
    mapped fragments considering only contigs with a minimum length of 3 kb.
  - `PHYL_Chlorobiaeceae_fastANI_1kb.tsv`: the matrix of the pairwise average nucleotide identity
    values calculated for a subset of the genomes of the family Chlorobiaceae considering only
    contigs with a minimum length of 1 kb
  - `PHYL_Chlorobiaeceae_fastANI_1kb_mapped_fragments.tsv`: The raw fastANI output with percentage of 
    mapped fragments considering onlyncontigs with a minimum length of 1 kb.
  - `PHYL_Flexilinea_proteintree_taxa.tsv`: the overview of the comparative genomes for the analysis
    of the family Anaerolineaceae
  - `PHYL_Flexilinea_proteintree_RAxML.tre`: the phylogenetic tree of the family Anaerolineaceae
    generated using RAxML
  - `PHYL_Flexilinea_fastANI.tsv`: the pairwise average nucleotide identity for all members of the
    genus *Flexilinea* and its neighbouring clades after discarding contigs shorter than 2 kb
  - `PHYL_Flexilinea_proteintree_RAxML_bootstrap.tre`: the phylogenetic tree of a subset of the
    family Anaerolineaceae generated using RAxML with the clade support inferred from a bootstrap
    analysis
  - `PHYL_Yersinia_completeness.tsv`: the overview of the fraction of sites per *Yersinia pestis*
    genome that could not be genotyped by Andrades Valtuena *et al.* (2022)
  - `PHYL_Yersinia_fastANI.tsv`: the pairwise average nucleotide identity for all members of the
    species *Yersinia pestis* analysed by Andrades Valtuena *et al.* (2022)
  - `PHYL_Yersinia_proteintree_RAxML_bootstrap.tre`: the phylogenetic tree of the *Yersinia pestis*
    genomes analysed by Andrades Valtuena *et al.* (2022) generated using RAxML with the clade
    support inferred from a bootstrap analysis
  - `PHYL_AFSA_corason.svg`: the result of the corason analysis shoowing the syteny of the 
    Butyrolactone BGC
  - `PHYL_AFSA_family_aa_tree.tre`: the resulting tree from phylophlan analysis using the AFSA 
    family cluster
  - `PHYL_BGC_antismash_presence_absence.txt`: the resulting table from teh antismash anaysis resulting
    in the presence absence of BGCs across modern and ancient genomes
  - `PHYL_BGC_coregenes_alignment.xlsx`: the percent identities of the Butyrolactone BGC core genes 
    using both amino acid and nucleotide sequences
  - `PHYL_bigscape_analysis.xlsx`: The resulting pairwise distances between the ancient and modern
    Butyrolactone BGCs from BIGSCAPE analysis
  - `PHYL_chlorobiales_assembly-accs.txt`: the list of accessions from the Chlorobiales genomes downloaded
    from NCBI
  - `PHYL_chlorobiales_drep_climicola_bootstrap_nucleotide_tree.tre`: the phylogenetic tree constructed using
    *Chlorobium limicola* as the core genes and the dereplicated modern reference Chlorobiales genomes and the
    ancient *Chlorobium* genomes on the nucleotide level 
  - `PHYL_chlorobiales_drep_climicola_bootstrap_protein_tree.tre`: the the phylogenetic tree constructed using
    *Chlorobium limicola* as the core genes and the dereplicated modern reference Chlorobiales genomes and the
    ancient *Chlorobium* genomes on the protein level
  - `PHYL_chlorobiales_drep_cparvum_bootstrap_protein_tree.tre`: the phylogenetic tree constructed using
    *CChlorobaculum parvum* as the core genes and teh dereplicated modern reference Chlorobiales genomes and the
    ancient *Chlorobium* genomes
  - `PHYL_chlorobiales_genomes_annotated_to_spp.txt`: the list of accession numbers of genomes used in the pruned
    Chlorobiales genomes phylogenetic tree
  - `PHYL_chlorobiales_genomes_RaxML_proteintree.tre`: the phylogenetic tree constructed using
    *Chlorobium limicola* as the core genes and the entire modern reference Chlorobiales genomes downloaded from
     the NCBI assembly database and the ancient *Chlorobium* genomes on the nucleotide level
  - `PHYL_eggnog_mapper.xlsx`: the resulting tables from the EGGNOG analysis that functionally classifies the
    contigs of all ancient *Chlorobium* MAGs and 4 modern reference genomes
  - `PHYL_pruned_ancient_clade_boostrap_proteintree.tre`: the pruned tree of Chlorobiales genomes that is
    constructed using all modern genomes classified upto a species level and the ancoent *Chlorobium* genomes
  - `PHYL_AFSA_clustered_family.faa`: a list of modern AFSA that contain a single A-factor domain resembling
    the ancient AFSA

6. `FUNC`: results related to the functional analyses

  - `FUNC_roary_pangenome.tsv`: the results of the pan-genome analysis of the 6 Chlorobium MAGs and
    the four published *Chlorobium limicola* genomes using Roary for identifying clade specific
    genes and determining their function using eggNOG
  - `FUNC_KEGG_specificKO.tsv`: the KEGG pathway annotations of all KEGG orthologs that were
    specific to either the *Chlorobium limicola* or the *Chlorobium* MAGs
