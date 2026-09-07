nextflow.enable.dsl=2

include { FASTQC } from '../modules/fastqc'
include { TRIMGALORE } from '../modules/trimgalore'
include { STAR_INDEX; STAR_ALIGN } from '../modules/star'
include { FEATURECOUNTS } from '../modules/subread'
include { MULTIQC } from '../modules/multiqc'
include { DESEQ2 } from '../modules/deseq2'
include { ENRICHR } from '../modules/enrichr'
include { IMMUCELLAI } from '../modules/immucellai'
include { PLOT_DECONVOLUTION } from '../modules/plot_deconvolution'
include { IMSIG } from '../modules/imsig'
include { ARRIBA } from '../modules/arriba'
include { RMATS } from '../modules/rmats'
include { SPLICING_PLOTS as RMATS_PLOT } from '../modules/splicing_plots'
include { SASHIMI_PLOT as RMATS_SASHIMI } from '../modules/splicing_sashimi_plot'

workflow RNA_SEQ_ANALYSIS {
    
    main:
        log.info "RNA-seq analysis started..."

        if (params.single_end) {
            ch_reads = Channel.fromPath(params.input_reads, checkIfExists: true)
                              .map {file -> tuple(file.simpleName, [file]) }
        } else {
            ch_reads = Channel.fromFilePairs(params.input_reads, checkIfExists: true)
        }

        FASTQC(ch_reads)
        TRIMGALORE(ch_reads)
        ch_fasta = file(params.fasta)
        ch_gtf   = file(params.gtf)
        
        STAR_INDEX(ch_fasta, ch_gtf)
        STAR_ALIGN(STAR_INDEX.out.index, TRIMGALORE.out.reads)
        
        ch_bams_raccolti = STAR_ALIGN.out.bam.map { it[1] }.collect()
        FEATURECOUNTS(ch_gtf, ch_bams_raccolti)

        // Inizializzazione Canali di Sicurezza 
        ch_arriba_fusions      = Channel.empty()
        ch_arriba_discarded    = Channel.empty()
        ch_arriba_plots        = Channel.empty()
        ch_arriba_multiqc      = Channel.empty()
        
        ch_deseq2_results      = Channel.empty()
        ch_deseq2_multiqc      = Channel.empty()
        ch_enrichr_results     = Channel.empty()
        ch_enrichr_multiqc     = Channel.empty()
        
        ch_immucellai_results  = Channel.empty()
        ch_deconvolution_plots = Channel.empty()
        ch_deconv_multiqc      = Channel.empty()
        ch_imsig_results       = Channel.empty()
        ch_imsig_plot          = Channel.empty()
        ch_imsig_multiqc       = Channel.empty()
        
        ch_rmats_results       = Channel.empty() 
        ch_rmats_plots         = Channel.empty()
        ch_rmats_sashimi       = Channel.empty()
        ch_rmats_multiqc       = Channel.empty()

        // ESECUZIONI CONDIZIONALI
        if (!params.skip_fusions) {
            ARRIBA(STAR_ALIGN.out.bam, ch_fasta, ch_gtf)
            ch_arriba_fusions   = ARRIBA.out.fusions
            ch_arriba_discarded = ARRIBA.out.discarded
            ch_arriba_plots     = ARRIBA.out.plots
            ch_arriba_multiqc   = ARRIBA.out.multiqc_counts
        }

if (!params.skip_differential) {
            DESEQ2(FEATURECOUNTS.out.counts, file(params.samplesheet))
            
            // Passiamo SOLTANTO i file filtrati a Enrichr, appiattiti
            ENRICHR(DESEQ2.out.filtered_results.flatten())
            
            // Aggiorna anche le variabili di output per la fine del workflow
            ch_deseq2_results  = DESEQ2.out.raw_results.mix(DESEQ2.out.complete_tables, DESEQ2.out.filtered_results, DESEQ2.out.results_pdf)
            ch_deseq2_multiqc  = DESEQ2.out.multiqc_png.mix(DESEQ2.out.pca_png)
            ch_enrichr_results = ENRICHR.out.enrichr_results
            ch_enrichr_multiqc = ENRICHR.out.multiqc_png
        }

        if (!params.skip_deconvolution) {
            IMMUCELLAI(FEATURECOUNTS.out.counts, params.immucell_ref)
            PLOT_DECONVOLUTION(IMMUCELLAI.out.fractions)
            IMSIG(FEATURECOUNTS.out.counts)
            
            ch_immucellai_results  = IMMUCELLAI.out.tpm_matrix.mix(IMMUCELLAI.out.fractions)
            ch_deconvolution_plots = PLOT_DECONVOLUTION.out.plots
            ch_deconv_multiqc      = PLOT_DECONVOLUTION.out.multiqc_png
            ch_imsig_results       = IMSIG.out.results
            ch_imsig_plot          = IMSIG.out.plot
            ch_imsig_multiqc       = IMSIG.out.multiqc_png
        }

        if (!params.skip_splicing) {
            RMATS(ch_bams_raccolti, file(params.samplesheet), ch_gtf)
            RMATS_PLOT(RMATS.out.splicing_results, 'rMATS')
            RMATS_SASHIMI(ch_bams_raccolti, file(params.samplesheet), RMATS.out.splicing_results)
            
            ch_rmats_results = RMATS.out.splicing_results
            ch_rmats_plots   = RMATS_PLOT.out.plots
            ch_rmats_sashimi = RMATS_SASHIMI.out.pdf_plots
            ch_rmats_multiqc = RMATS_PLOT.out.multiqc_png
        }


        // MULTIQC
        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yaml", checkIfExists: true)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQC.out.zip,
            TRIMGALORE.out.log,
            STAR_ALIGN.out.log,
            FEATURECOUNTS.out.summary,
            ch_deseq2_multiqc,
            ch_enrichr_multiqc,
            ch_arriba_multiqc,
            ch_imsig_multiqc,
            ch_deconv_multiqc,
            ch_rmats_multiqc
        )

        MULTIQC( ch_multiqc_files.collect(), ch_multiqc_config )

    emit:
        fastqc_results        = FASTQC.out.html.mix(FASTQC.out.zip)
        trimgalore_results    = TRIMGALORE.out.reads.mix(TRIMGALORE.out.log)
        star_index_results    = STAR_INDEX.out.index
        star_align_results    = STAR_ALIGN.out.bam.mix(STAR_ALIGN.out.log)
        featurecounts_results = FEATURECOUNTS.out.counts.mix(FEATURECOUNTS.out.summary)
        multiqc_results       = MULTIQC.out.report.mix(MULTIQC.out.data)
        
        arriba_fusions        = ch_arriba_fusions.flatten()
        arriba_discarded      = ch_arriba_discarded.flatten()
        arriba_plots          = ch_arriba_plots.flatten()
        
        deseq2_results        = ch_deseq2_results.flatten()
        enrichr_results       = ch_enrichr_results.flatten()
        
        immucellai_results    = ch_immucellai_results.flatten()
        deconvolution_plots   = ch_deconvolution_plots.flatten()
        imsig_results         = ch_imsig_results.flatten()
        imsig_plot            = ch_imsig_plot.flatten()
        
        rmats_results         = ch_rmats_results.flatten()
        rmats_plots           = ch_rmats_plots.flatten()
        rmats_sashimi         = ch_rmats_sashimi.flatten()
}
