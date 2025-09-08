#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

option_list <- list(
  make_option("--counts", type="character"),
  make_option("--meta",   type="character"),
  make_option("--outdir", type="character", default=".")
)
opt <- parse_args(OptionParser(option_list=option_list))

counts_path <- opt$counts
meta_path   <- opt$meta
outdir      <- opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Read featureCounts table
fc <- read.delim(counts_path, comment.char="#", check.names=FALSE)
# featureCounts format: first columns are annotation; samples start at column 7 by default
ann_cols <- 1:6
cts <- fc[ , -(ann_cols)]
rownames(cts) <- fc$Geneid

# Metadata (CSV with columns: sample,condition)
meta <- read.csv(meta_path, header=FALSE, col.names=c("sample","condition"), stringsAsFactors=TRUE)
# Ensure sample order matches counts columns
cts <- cts[ , meta$sample, drop=FALSE]

# edgeR
group <- meta$condition
dge <- DGEList(counts=cts, group=group)
keep <- filterByExpr(dge, group=group)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# If two groups: first vs second
if (ncol(design) == 2) {
  cn <- paste0(colnames(design)[1], "-", colnames(design)[2])
  ql <- glmQLFTest(fit, contrast=makeContrasts(!!as.name(cn), levels=design))
} else {
  stop("edgeR script expects exactly 2 conditions for this template.")
}

tab <- topTags(ql, n=Inf)$table
tab$gene_id <- rownames(tab)
tab <- tab[ , c("gene_id","logFC","logCPM","F","PValue","FDR")]
write.csv(tab, file=file.path(outdir,"edgeR_DE_results.csv"), row.names=FALSE)

# Quick MDS plot
png(file.path(outdir,"edgeR_MDS.png"), width=1200, height=900, res=150)
plotMDS(dge, labels=meta$sample, col=as.integer(group))
legend("topright", legend=levels(group), col=seq_along(levels(group)), pch=16)
dev.off()

# Quick MA plot
png(file.path(outdir,"edgeR_MAplot.png"), width=1200, height=900, res=150)
with(tab, {
  plot(logCPM, logFC, pch=16, col=ifelse(FDR < 0.05, "red", "grey60"),
       xlab="logCPM", ylab="logFC", main="edgeR MA-plot")
  abline(h=0, col="blue", lty=2)
})
dev.off()

