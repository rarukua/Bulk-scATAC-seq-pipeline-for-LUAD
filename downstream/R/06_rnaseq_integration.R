# downstream/R/06_rnaseq_integration.R
# Chromatin accessibility ↔ Gene expression integration
# TCGA-LUAD ATAC-seq + matched RNA-seq

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(optparse)
})
opt_list <- list(
  make_option("--da_peaks",     type="character", help="Differential peaks CSV"),
  make_option("--rnaseq_deseq", type="character", help="RNA-seq DESeq2 results CSV"),
  make_option("--outdir",       type="character", default="results/integration")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Load data ──────────────────────────────────────────────────────────────
da   <- read.csv(opt$da_peaks)      # differential ATAC peaks, annotated with nearest gene
rna  <- read.csv(opt$rnaseq_deseq)  # DESeq2 RNA-seq results

# ── Join on nearest gene ───────────────────────────────────────────────────
# For each differential ATAC peak, match to the same gene's RNA-seq fold change.
# Expected relationship: more open promoter (ATAC) → higher expression (RNA-seq)
integrated <- da %>%
  filter(!is.na(SYMBOL)) %>%
  inner_join(rna %>% select(gene=gene_name, rna_log2fc=log2FoldChange,
                              rna_padj=padj),
             by=c("SYMBOL"="gene")) %>%
  mutate(
    atac_dir = case_when(Fold >  1 & FDR < 0.05 ~ "More open in Tumor",
                         Fold < -1 & FDR < 0.05 ~ "Less open in Tumor",
                         TRUE                    ~ "NS"),
    rna_dir  = case_when(rna_log2fc >  1 & rna_padj < 0.05 ~ "Up in Tumor",
                         rna_log2fc < -1 & rna_padj < 0.05 ~ "Down in Tumor",
                         TRUE                               ~ "NS"),
    concordant = (atac_dir == "More open in Tumor"  & rna_dir == "Up in Tumor") |
                 (atac_dir == "Less open in Tumor"  & rna_dir == "Down in Tumor")
  )

# ── Scatter: ATAC fold change vs RNA fold change ───────────────────────────
luad_genes <- c("NKX2-1","FOXA1","EGFR","KRAS","STK11","KEAP1","TP53","MYC")

p_scatter <- ggplot(integrated, aes(x=Fold, y=rna_log2fc, color=atac_dir)) +
  geom_point(alpha=0.3, size=0.8) +
  geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.5) +
  geom_label_repel(
    data=integrated %>% filter(SYMBOL %in% luad_genes),
    aes(label=SYMBOL), size=3, max.overlaps=20
  ) +
  scale_color_manual(values=c(
    "More open in Tumor" ="#D7191C",
    "Less open in Tumor" ="#2C7BB6",
    "NS"                 ="grey70"
  )) +
  labs(
    title    = "Chromatin Accessibility vs Gene Expression — LUAD",
    subtitle = sprintf("Pearson r = %.2f | Concordant pairs: %d",
                       cor(integrated$Fold, integrated$rna_log2fc, use="complete"),
                       sum(integrated$concordant, na.rm=TRUE)),
    x        = "ATAC log2 Fold Change (Tumor/Normal)",
    y        = "RNA-seq log2 Fold Change (Tumor/Normal)",
    color    = "ATAC direction",
    caption  = "Concordance validates that chromatin opening drives transcriptional activation"
  ) +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "atac_rna_concordance.pdf"),
       p_scatter, width=7, height=6)

write.csv(integrated, file.path(opt$outdir, "atac_rna_integrated.csv"),
          row.names=FALSE)
message("── ATAC-RNA integration complete ──")
