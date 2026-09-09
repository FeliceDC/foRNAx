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

<img width="2624" height="1682" alt="mermaid-diagram" src="https://github.com/user-attachments/assets/d12323ad-bb0e-469f-bf15-334c0a9c1245" />


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
>If you want, you can run Deseq2 with two variables. Then you have to write --design "variable1 + variable2"


 | Flag | Description |
| :--- | :--- |
| `--input_reads` | Serve per specificare il percorso in cui si trovano i file fastq.gz |
| `--fasta` | Serve per specificare il percorso in cui si trova il genoma di riferimento |
| `--gtf` | Serve per specificare il percorso in cui si trova il file per le annotazioni genomiche |
| `--samplesheet` | Serve per specificare il percorso in cui si trova il samplesheet di riferimento |

**optional flags**

 | Flag | Description |
| :--- | :--- |
| `--skip_fusions` | Salta l'analisi delle fusioni geniche (Arriba) |
| `--skip_deconvolution` | Salta la stima dell'infiltrazione cellulare |
| `--skip_differential` | Salta l'analisi differenziale (DESeq2/EnrichR) |
| `--skip_splicing` | Salta l'analisi dello splicing alternativo (RMats/LeafCutter) |
| `--g gene_name` | Questo parametro riguarda l'esecuzione di Feature Count. La pipeline è impostata di default ad usare `gene_id``. Se invece si desidera ottenere questo tipo di output, aggiungere questo parametro |
| `--strandedness` | La pipeline è impostata di default per utilizzare strandedness 0, ma puoi selezionare a piacere tra 0, 1 e 2 |
| `--immucell_ref blood` | Permette di impostare il database di riferimento utilizzato dalla pipeline per realizzare la deconvoluzione immunitaria di Imsig. Di base, il tool utilizza “load_tumor_reference_data()”, se invece si desidera  “load_blood_reference_data()” bisogna utilizzare questo parametro |
| `--deseq2_pvalue` | Permette di impostare il valore di p-value per il filtraggio dei risultati di Deseq2. Di base, la pipeline utilizza pvalue=0.05 |
| `--deseq2_logfc` | Permette di impostare il valore di log fold change per il filtraggio dei risultati di Deseq2. Di base, la pipeline utilizza logfc=1.5 |
| `--single_end true` |  La pipeline è impostata per gestire autonomamente dati del tipo PAIRED END. Se invece si stanno usando dati del tipo SINGLE END aggiungere questo parametro |
| `--enrichr_database` | Spericificare i database da voler utilizzare tra quelli nativamente supportati da EnrichR |
| `--outdir` | Puoi specificare una cartella diversa in cui verranno salvati i risultati |
| `--max_cpus` | Puoi selezionare il numero massimo di cpu da impiegare. Il valore di default è impostato a 16 |
| `--max_memory` | Puoi selezionare la quantità di ram da impiegare. Il valore di default è impostato a 64 GB |
| `--max_time` | Puoi selezionare la quantità massima di tempo da dedicare all'analisi. Il valore di default è impostato a 24 h |




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
