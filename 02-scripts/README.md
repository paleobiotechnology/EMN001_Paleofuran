## Overview of the folder `02-scripts`

This folder contains all the scripts and Jupyter notebooks to prepare and analyse the data.

1. `PREP`: preparing of the data

  - `PREP_download_from_ENA.Snakefile`: download the sequencing data from ENA with the help of the
    Python program [ffq](https://github.com/pachterlab/ffq)
  - `PREP_preprocessing_dentalcalculus_sequencing_data.Snakefile`: prepare input file for
    pre-processing the dental calculus sequecing data with the Nextflow pipeline
    [nf-core/eager](https://nf-co.re/eager) and extract and merge all non-human reads per individual
  - `PREP_preprocessing_publishedNeanderthalcalculus.Snakefile`: download the sequencing data of the
    four Neanderthal dental calculus samples published by Weyrich *et al.* (2017) from ENA,
    pre-process the sequencing data using the Nextflow pipeline
    [nf-core/eager](https://nf-co.re/eager) and extract the non-human reads per individual
  - `PREP_preprocessing_ElMiron_sediments.Snakefile`: download the sequencing data of the six El
    Miron sediment samples generated for this study from ENA, pre-process the sequencing data using
    the Nextflow pipeline [nf-core/eager](https://nf-co.re/eager) and extract the non-human reads
    per individual

2. `ASMB`: analyses related to the *de novo* assembly and the non-reference binning

  - `ASMB_denovo_assembly_binning.Snakefile`: prepare the sample table and config file for
    assembling the ancient and modern DNA samples separately and run the assembly pipeline
    https://github.com/alexhbnr/ancient_metagenome_assembly

3. `QUAL`: analyses related to the data quality

  - `QUAL_fragmentlengths.Snakefile`: summarise the DNA molecule fragment length across the
    sequencing data generated on 2x 75 bp paired-end Illumina sequencing runs for all samples
  - `QUAL_freeBayes_nomismatches_MEGAHIT.Snakefile`: evaluate the number of substitutions (per
    type), the minor allele frequencies of the variant calls, and the number of affected contigs for
    which freeBayes inferred a different major allele than it was exported by MEGAHIT from the
    assembly graph
  - `QUAL_pyDamage_summary.Snakefile`: summarise the pyDamage results regarding the fraction of
    contigs that are considered "ancient" by pyDamage and the distribution of the predicted accuracy
  - `QUAL_refalignment_Chlorobiaceae.Snakefile`: align the short-read sequencing data of all dental
    calculus samples against the contigs of EMN001 and evaluate the number of aligned reads against
    the contigs that were assigned to the Chlorobiaceae MAG

4. `PUBL`: preparing the final figures and tables

  - `PUBL_Figure1.Snakefile`: plot an overview map of all samples, the DNA molecule length
    and contig length distributions for the samples EMN001, PES001, and PLV001, the overview of the
    reconstructed MAGs for the five samples with a HQ Chlorobiaceae MAG, and the aDNA damage for two
    MAGs of EMN001
  - `PUBL_Dataset_S1.Snakefile`: combine the tables containing the information about the samples,
    the sequencing data, and the number of DNA molecules available for *de novo* assembly into a
    single XSLX file
  - `PUBL_Dataset_S2.Snakefile`: combine the tables containing the results about the de novo
    assembly of the ancient and modern dental calculus samples
  - `PUBL_Dataset_S4.Snakefile`: combine the tables containing the results on the alignment of the
    short-read sequencing data against the contigs assembled from EMN001
  - `PUBL_FigureS_freeBayescorrection.Snakefile`: plot the supplementary figure summarising the
    results of the substitution types and minor allele frequencies that were observed when
    correcting the contig sequences using freeBayes
  - `PUBL_FigureS_denovoassembly_metaQUAST_calN50.Snakefile`: plot the supplementary figure
    summarising the evaluation of the *de novo* assembly performance of the ancient and modern
    dental calculus samples
  - `PUBL_FigureS_pyDamage_contigs.Snakefile`: plot the supplementar figure summarising the results
    of pyDamage on contig level

X. `ENVS`: conda environments for the use with Snakemake

  - `ENVS_bioawk.yaml`: [bioawk](https://github.com/lh3/bioawk)
  - `ENVS_bowtie2.yaml`: [bowtie2](https://github.com/BenLangmead/bowtie2)
  - `ENVS_fastp.yaml`: [fastp](https://github.com/OpenGene/fastp)
  - `ENVS_picard.yaml`: [picard](https://broadinstitute.github.io/picard)
  - `ENVS_samtools.yaml`: [samtools](https://github.com/samtools/samtools)
