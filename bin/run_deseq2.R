#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
metadata_file <- args[2]
user_design <- args[3]
user_pvalue <- as.numeric(args[4])
user_logfc  <- as.numeric(args[5])

library(DESeq2)
library(ggplot2)

counts <- read.table(counts_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, row.names = 1, stringsAsFactors = TRUE)

counts <- counts[, 6:ncol(counts)]

col_names <- colnames(counts)
col_names <- sub("^X", "", col_names)
col_names <- sub("\\..*", "", col_names)

true_samples <- rownames(meta)

for (ts in true_samples) {
    match_idx <- grep(paste0("^", ts), col_names)
    if (length(match_idx) == 0) {
        match_idx <- grep(ts, col_names, fixed = TRUE)
    }
    if (length(match_idx) == 1) {
        col_names[match_idx] <- ts
    }
}

colnames(counts) <- col_names
common_samples <- intersect(colnames(counts), rownames(meta))

if (length(common_samples) == 0) {
    message <- paste("\n\n#####################################################\n",
                       "ERROR: The sample names do not correspond!\n",
                       "Matrix names: ", paste(colnames(counts), collapse=" , "), "\n",
                       "Samplesheet names: ", paste(rownames(meta), collapse=" , "), "\n",
                       "#####################################################\n\n")
    stop(message)
}

counts <- counts[, common_samples]
meta <- meta[common_samples, , drop=FALSE]

design_formula <- as.formula(paste("~", user_design))
pca_groups <- trimws(unlist(strsplit(user_design, "\\+")))
main_condition <- pca_groups[length(pca_groups)]


meta[[main_condition]] <- as.factor(meta[[main_condition]])

if (length(unique(as.list(counts))) == 1) {
    message("\n========================================================")
    message("DETECTED: Samples Are Identical (Profile Test).")
    message("Adding artificial batch effect in order to test the pipeline.")
    message("========================================================\n")
    set.seed(42)
    rumore <- matrix(rpois(nrow(counts) * ncol(counts), lambda = 5), nrow = nrow(counts))
    counts <- counts + rumore
}

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- tryCatch({
    DESeq(dds, sfType="poscounts")
}, error = function(e) {
    message("\nStandard curve fitting failed (too few genes).")
    message("Unpacking the algorithm and forcing manual calculations...\n")
    dds <- estimateSizeFactors(dds, type="poscounts")
    dds <- estimateDispersionsGeneEst(dds)
    dispersions(dds) <- mcols(dds)$dispGeneEst
    dds <- nbinomWaldTest(dds)
    return(dds)
})

# Pre-calcolo VST/NormTransform per PCA e Heatmap
if (nrow(dds) < 1000) {
    message("Less than 1000 genes: using normTransform instead of vst for plots")
    vsd <- normTransform(dds)
} else {
    vsd <- vst(dds, blind=FALSE)
}

# --- PLOTS ---
pdf("deseq2_plots.pdf", width = 16, height = 7)

pca_data <- plotPCA(vsd, intgroup = colnames(meta), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

for (var_col in colnames(meta)) {
    pca_plot <- ggplot(pca_data, aes_string(x = "PC1", y = "PC2", color = var_col)) +
        geom_point(size = 4, alpha = 0.85) +
        geom_text(aes(label = name), size = 3.5, vjust = -1, show.legend = FALSE) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(paste("PCA - Colored by", var_col)) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            legend.title = element_text(face = "bold", size = 11),
            legend.text = element_text(size = 10),
            axis.title = element_text(size = 11, face = "bold")
        )
    print(pca_plot)
    ggsave(paste0("deseq2_pca_", var_col, "_mqc.png"), plot = pca_plot, width = 10, height = 7, dpi = 300)
}

# --- ANALISI PER TUTTI I CONTRASTI ---
cond_levels <- levels(dds[[main_condition]])
contrast_pairs <- combn(cond_levels, 2, simplify = FALSE)

plot_colors <- c("#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a", "#b15928")
dynamic_colors <- rep(plot_colors, length.out = length(cond_levels))

for (pair in contrast_pairs) {
    cond1 <- pair[1]
    cond2 <- pair[2]
    c_name <- paste0(cond1, "_vs_", cond2)
    message(paste("\nProcessing contrast:", c_name))
    
    # 1. Estrazione risultati
    res <- results(dds, contrast=c(main_condition, cond1, cond2))
    
    # Salvataggio Tabelle
    write.table(as.data.frame(res), file=paste0("deseq2_results_", c_name, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
    
    res_clean <- res[!is.na(res$padj), ]
    res_filt <- res_clean[res_clean$padj < user_pvalue & abs(res_clean$log2FoldChange) > user_logfc, ]
    res_filt_df <- as.data.frame(res_filt)
    res_filt_df$Gene_Name <- rownames(res_filt_df)
    write.table(res_filt_df, file=paste0("filtered_results_", c_name, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
    
    normalized_counts <- counts(dds, normalized=TRUE)
    complete_table <- merge(as.data.frame(res), normalized_counts, by="row.names", all=TRUE)
    colnames(complete_table)[1] <- "Gene_Name"
    complete_table <- complete_table[order(complete_table$padj), ]
    write.table(as.data.frame(complete_table), file=paste0("complete_table_", c_name, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
    
    # 2. MA Plot
    plotMA(res, main=paste("MA Plot:", c_name))
    
    # 3. Volcano Plot (dentro il PDF principale)
    max_fc <- max(abs(res$log2FoldChange), na.rm=TRUE)
    limite_x <- ifelse(is.finite(max_fc), max_fc * 1.1, 5) 
    
    with(res, plot(log2FoldChange, -log10(padj), pch=20, main=paste("Volcano Plot:", c_name), col="darkgrey", xlim=c(-limite_x, limite_x)))
    with(subset(res, padj < user_pvalue & abs(log2FoldChange) > user_logfc), points(log2FoldChange, -log10(padj), pch=20, col="red"))
    abline(v=c(-user_logfc, user_logfc), col="blue", lty=2)
    abline(h=-log10(user_pvalue), col="blue", lty=2)
    
    # Export Volcano PNG (separato per MultiQC)
    png(paste0("deseq2_volcano_", c_name, "_mqc.png"), width = 1400, height = 1100, res = 150)
    with(res, plot(log2FoldChange, -log10(padj), pch=20, main=paste("Volcano Plot:", c_name), col="darkgrey", xlim=c(-limite_x, limite_x)))
    with(subset(res, padj < user_pvalue & abs(log2FoldChange) > user_logfc), points(log2FoldChange, -log10(padj), pch=20, col="red"))
    abline(v=c(-user_logfc, user_logfc), col="blue", lty=2)
    abline(h=-log10(user_pvalue), col="blue", lty=2)
    dev.off()
    
 # 4. Heatmap Top 50 
    top_genes <- head(order(res$padj), 50)
    if (length(top_genes) > 1) {
        mat <- assay(vsd)[top_genes, ]
        mat <- mat - rowMeans(mat)
        colori_heatmap <- colorRampPalette(c("blue", "white", "red"))(256)
    
        heatmap(mat, scale="none", col=colori_heatmap, margins=c(12, 12), cexCol=0.3, cexRow=0.5, main=paste("Heatmap Top 50:", c_name))
    }

    
    # 5. Top 6 geni (specifici per il contrasto)
    top6_genes <- head(order(res$padj), 6)
    if (length(top6_genes) > 0) {
        par(mfrow=c(3,2), las=1)
        top6_counts <- counts(dds, normalized=TRUE)[top6_genes, , drop = FALSE]
        absolute_max <- max(top6_counts)
        y_limit <- c(0.5, absolute_max + (absolute_max * 0.1))
        
        for (i in top6_genes) {
            gene_name <- rownames(res)[i]
            plotCounts(dds, 
                   gene = gene_name, 
                   intgroup = main_condition, 
                   main = paste("Exp:", gene_name, "\n(", c_name, ")"),
                   col = dynamic_colors[as.numeric(dds[[main_condition]])], 
                   pch = 16,                             
                   cex = 1.5,
                   ylim = y_limit 
            )
        }
        par(mfrow=c(1,1))
    }
}

dev.off()
