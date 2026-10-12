<h1 align="center">
   foRNAx
</h1>

<p align="center">
  <img src="https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg" alt="Nextflow">
  <img src="https://img.shields.io/badge/status-active-success.svg" alt="Status">
  <img src="https://img.shields.io/badge/Bioinformatics-RNA--seq-blue" alt="Bioinformatics">
</p>

## Introduction
**FeliceDC/foRNAx** is a comprehensive, modular bioinformatics analysis pipeline used for RNA sequencing data. Developed in Nextflow (DSL2), it automates the entire workflow from raw FASTQ reads to advanced downstream analysis (Differential Expression, Splicing, Fusions, and Deconvolution), ensuring reproducibility and scalability.

The pipeline is built using Docker containers, meaning you don't need to install any bioinformatics tools manually.

<img width="960" height="509" alt="image" src="https://github.com/user-attachments/assets/a3bf5e5d-2e2d-43e2-9ecc-331df4f22128" />

## Pipeline Summary
1. Raw read QC (`FastQC`)
2. Adapter and quality trimming (`Trim Galore!`)
3. Read alignment and indexing (`STAR`)
4. Gene-level quantification (`featureCounts`)
5. Pipeline QC report (`MultiQC`)
6. Differential Expression Analysis (`DESeq2`), followed by pathway enrichment (`EnrichR`)
7. Tumor Deconvolution: Immune and stromal cell infiltration estimation ((`ImmuCellAI`) and (`ImSig`)).
8. Alternative Splicing: Classical statistical splicing analysis (`rMATS`) complete with automated Volcano, Bar, and Sashimi plots.
9. Gene Fusions: Structural variant detection (`Arriba`)


## Usage

To run the pipeline on your own samples, you need to provide:
1. Your raw fastq.gz files
2. A reference genome
3. An annotation file
4. A design matrix (named "samplesheet").

The samplesheet must be a comma-separated values file (.csv). The first column (called "sample") must match the FASTQ file names.The other column should contain the variable you want to use for differential analysis. Optionally, you can perform differential analysis using mutiple variables if needed.
Example:

**samplesheet.csv**
```bash
sample,condition,age,library_selection
SRR8518319,normal_adiacent,52,cDNA
SRR8518327,tumor,37,cDNA
SRR8518335,normal_adiacent,62,cDNA
SRR8518360,tumor,54,cDNA
```

Once the samplesheet has been created, make sure you have ready the samplesheet, FASTQ files, genome, and GTF files paths.

Now you should be ready to run the pipeline.
>[!NOTE]
>An example running code is
>```bash
>Nextflow run FeliceDC/foRNAx --input_reads "/Your/Files/Path/*fastq.gz" --fasta "/Your/Genome/Path" --gtf "/Your/Annotations/Path" --design "condition" --samplesheet "/Your/File/Path"
>```
>
>If you want, you can run Deseq2 with two variables. Then you have to write --design "variable1 + variable2" <p>


There are two list of the flags that you have to specify in the command string. <p>
**Necessary flags** <p>
 | Flag | Description |
| :--- | :--- |
| `--input_reads` | Specifies the path to the directory containing the `fastq.gz` files |
| `--fasta` | Specifies the path to the reference genome file |
| `--gtf` | Specifies the path to the genomic annotation file |
| `--samplesheet` | Specifies the path to the samplesheet file |

**Optional flags**  <p>

 | Flag | Description |
| :--- | :--- |
| `--skip_fusions true` | Skip the gene fusion analysis (Arriba) |
| `--skip_deconvolution true` | Skip the estimation of cellular infiltration (ImSig / ImmucellAI |
| `--skip_differential true` | Skip the differential and pathways analysis (DESeq2 / EnrichR) |
| `--skip_splicing true` | Skip the alternative splicing analysis (rMATS-turbo) |
| `--g gene_name` | This parameter controls the execution of FeatureCounts. By default, the pipeline is configured to use `gene_id`. If this type of output is desired instead, add this parameter |
| `--strandedness 0/1/2` | The pipeline is configured by default to use strandedness 0, but you can select any value among 0, 1, and 2 |
| `--immucell_ref blood` | Allows you to specify the reference database used by the pipeline for ImSig immune deconvolution. By default, the tool uses `load_tumor_reference_data()`. To use `load_blood_reference_data()` instead, set this parameter |
| `--deseq2_pvalue` | Allows you to set the p-value threshold used to filter DESeq2 results. By default, the pipeline uses `pvalue = 0.05` |
| `--deseq2_logfc` | Allows you to set the log fold change threshold used to filter DESeq2 results. By default, the pipeline uses `logfc = 1.5` |
| `--single_end true` |  The pipeline is configured to automatically handle PAIRED-END data. If you are using SINGLE-END data instead, add this parameter |
| `--enrichr_database` | Specifies which databases to use among those natively supported by Enrichr |
| `--outdir` | Allows you to specify a different directory where the results will be saved |
| `--max_cpus` | Allows you to specify the maximum number of CPUs to use. The default value is set to 16 |
| `--max_memory` | Allows you to specify the amount of RAM to allocate. The default value is set to 64 GB |
| `--max_time` | Allows you to specify the maximum amount of time allocated for the analysis. The default value is set to 24 hours |




</select>

>[!WARNING]
>Running the pipeline on full human datasets requires significant computational resources. It is highly recommended to check your machine and specify an appropriate --max_cpus limit.

## Output Structure
By default, the pipeline creates a results/ directory containing the following sub-directories:

- fastqc/ and multiqc/: Interactive HTML quality reports.

- star/: Sorted .bam files ready for IGV visualization.

- featurecounts/: Raw count matrices.

- deseq2/: CSV tables with statistically significant Differentially Expressed Genes (DEGs) and related plots (MA plot, PCA, Volcano plot, Heatmap).

- enrichr/: Pathway enrichment tables and bar plots.

- splicing/: Separated results for rmats/ and darts_ai/, including raw tables, .pdf summary plots, and genomic Sashimi plots.

- fusions/: Arriba fusion tables and circular/linear .pdf visualizations.

- deconvolution/: Infiltration abundance matrices and comparative plots from ImmuCellAI and ImSig.


## Author

**Felice Di Casola** <p>
Laboratory of Molecular Medicine and Genomics, Department of Medicine, Surgery and Dentistry "Scuola Medica Salernitana", University of Salerno, 84081, Baronissi, SA, Italy.
