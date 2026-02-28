# downstream/R/04_tobias_footprinting.R
# TF Footprinting analysis using TOBIAS output
# ─────────────────────────────────────────────────────────────
# TOBIAS BINDetect classifies each TF motif occurrence as:
#   BOUND    — protected from Tn5 within open region (footprint present)
#   UNBOUND  — accessible but no TF protection
#   CHANGED  — differential binding between Tumor and Normal
#
# This is fundamentally different from motif ENRICHMENT (HOMER):
#   HOMER:   "is this motif sequence over-represented in open peaks?"
#   TOBIAS:  "is this TF actually bound at its motifs right now?"
#
# An enriched motif (HOMER) with no footprint (TOBIAS) means the
# chromatin is open but the TF is NOT bound there — a key distinction
# for understanding active vs permissive regulatory states.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(tidyr)
  library(optparse)
})

opt_list <- list(
  make_option("--bindetect_dir", type="character", help="TOBIAS BINDetect output dir"),
  make_option("--samplesheet",   type="character"),
  make_option("--outdir",        type="character", default="results/footprinting")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Load BINDetect results ─────────────────────────────────────────────────
# BINDetect_results.txt: one row per TF motif
# Columns: motif_id, name, output_prefix, motif_file,
#          <condition1>_mean_score, <condition2>_mean_score,
#          <condition1>_bound, <condition2>_bound,
#          change, pvalue (binomial test)
message("── Loading TOBIAS BINDetect results ──")

bindetect_file <- file.path(opt$bindetect_dir, "bindetect_results.txt")
bd <- read_tsv(bindetect_file, comment="#")

# Standardize column names (TOBIAS uses condition names in col headers)
# Assumes Tumor vs Normal comparison
tumor_score_col  <- grep("Tumor.*mean_score",  colnames(bd), value=TRUE)[1]
normal_score_col <- grep("Normal.*mean_score", colnames(bd), value=TRUE)[1]
tumor_bound_col  <- grep("Tumor.*bound",       colnames(bd), value=TRUE)[1]
normal_bound_col <- grep("Normal.*bound",      colnames(bd), value=TRUE)[1]

bd <- bd %>%
  rename(
    tumor_score  = all_of(tumor_score_col),
    normal_score = all_of(normal_score_col),
    tumor_bound  = all_of(tumor_bound_col),
    normal_bound = all_of(normal_bound_col)
  ) %>%
  mutate(
    score_change  = tumor_score - normal_score,
    pvalue_adj    = p.adjust(pvalue, method="BH"),
    tf_status = case_when(
      pvalue_adj < 0.05 & score_change >  0.1 ~ "Gained in Tumor",
      pvalue_adj < 0.05 & score_change < -0.1 ~ "Lost in Tumor",
      TRUE                                     ~ "Unchanged"
    )
  )

message(sprintf("  TFs gained in tumor: %d", sum(bd$tf_status=="Gained in Tumor")))
message(sprintf("  TFs lost in tumor:   %d", sum(bd$tf_status=="Lost in Tumor")))

# LUAD-relevant TFs to highlight
luad_tfs <- c("NKX2-1","FOXA1","FOXA2","TP63","SOX2","CEBPA",
              "AP-1","JUNB","FOSL2","KLF4","NFI","TEAD")

# ── Plot 1: BINDetect volcano plot ─────────────────────────────────────────
# X-axis: change in footprint score (tumor - normal)
# Y-axis: -log10(adjusted p-value)
# High footprint score change + significant = TF with differential binding
p_binddetect <- ggplot(bd,
  aes(x=score_change, y=-log10(pvalue_adj), color=tf_status)) +
  geom_point(alpha=0.5, size=1.2) +
  geom_label_repel(
    data=bd %>% filter(tf_status != "Unchanged",
                       grepl(paste(luad_tfs, collapse="|"), name, ignore.case=TRUE)),
    aes(label=name), size=3, max.overlaps=20
  ) +
  scale_color_manual(values=c(
    "Gained in Tumor" = "#D7191C",
    "Lost in Tumor"   = "#2C7BB6",
    "Unchanged"       = "grey70"
  )) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(0.05),  linetype="dashed", color="grey40") +
  labs(
    title    = "TOBIAS BINDetect: TF Footprint Changes in LUAD",
    subtitle = "X-axis = change in footprint score (positive = more binding in tumor)\nEach point = one TF motif",
    x        = "Footprint score change (Tumor − Normal)",
    y        = "-log10(adjusted p-value)",
    color    = NULL,
    caption  = "Footprinting distinguishes BOUND from UNBOUND motifs — unlike motif enrichment alone"
  ) +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "bindetect_volcano.pdf"),
       p_binddetect, width=9, height=6)

# ── Plot 2: Footprint profile at key LUAD TF ──────────────────────────────
# Each TF directory in BINDetect output contains:
#   <TF>_overview.png: aggregate footprint profile across all bound sites
#   This is the "footprint" — a dip in Tn5 signal at the motif center
#   flanked by high accessibility = hallmark of a bound TF
plot_footprint <- function(tf_name, bindetect_dir) {
  # Find TF directory
  tf_dirs <- list.dirs(bindetect_dir, recursive=FALSE)
  tf_dir  <- tf_dirs[grepl(tf_name, basename(tf_dirs), ignore.case=TRUE)][1]
  if (is.na(tf_dir)) return(NULL)

  # Read per-position footprint scores (TOBIAS outputs these as bigWig)
  # Here we read the aggregate score file if available
  agg_file <- file.path(tf_dir, paste0(tf_name, "_aggregate.txt"))
  if (!file.exists(agg_file)) return(NULL)

  agg <- read_tsv(agg_file, col_names=c("position","tumor","normal"),
                  comment="#")

  agg_long <- agg %>%
    pivot_longer(cols=c(tumor, normal),
                 names_to="condition", values_to="score") %>%
    mutate(condition=factor(condition, levels=c("normal","tumor")))

  ggplot(agg_long, aes(x=position, y=score, color=condition)) +
    geom_line(linewidth=0.8) +
    geom_vline(xintercept=0, linetype="dashed", color="black", linewidth=0.3) +
    scale_color_manual(values=c(tumor="#D7191C", normal="#2C7BB6")) +
    labs(title=sprintf("TF Footprint: %s", tf_name),
         subtitle="Dip at center = nucleosome-free + TF protection from Tn5",
         x="Position relative to motif center (bp)",
         y="TOBIAS footprint score",
         color="Condition") +
    theme_bw(base_size=12)
}

# Generate footprint profiles for top LUAD TFs
for (tf in c("NKX2-1","FOXA1","AP-1")) {
  p <- plot_footprint(tf, opt$bindetect_dir)
  if (!is.null(p)) {
    ggsave(file.path(opt$outdir, sprintf("footprint_%s.pdf", tf)),
           p, width=6, height=4)
  }
}

# ── Plot 3: Heatmap — TF binding across samples ──────────────────────────
# For each sample, TOBIAS reports per-TF bound fraction
# This heatmap shows which TFs are differentially active across samples
sample_scores <- bd %>%
  filter(tf_status != "Unchanged") %>%
  arrange(score_change) %>%
  select(name, tumor_score, normal_score) %>%
  pivot_longer(cols=c(tumor_score, normal_score),
               names_to="condition", values_to="score") %>%
  mutate(condition=sub("_score","",condition))

p_heatmap_data <- bd %>%
  filter(tf_status != "Unchanged") %>%
  slice_max(abs(score_change), n=30) %>%
  select(name, tumor_score, normal_score) %>%
  tidyr::pivot_longer(cols=-name, names_to="condition", values_to="score") %>%
  mutate(condition = recode(condition,
    "tumor_score"="Tumor", "normal_score"="Normal"))

p_tf_heatmap <- ggplot(p_heatmap_data,
                        aes(x=condition, y=reorder(name, score), fill=score)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C",
                       midpoint=median(p_heatmap_data$score)) +
  labs(title="Top 30 Differentially Active TFs (TOBIAS Footprint Score)",
       x=NULL, y="Transcription Factor", fill="Footprint\nscore") +
  theme_bw(base_size=11) +
  theme(axis.text.y=element_text(size=8))

ggsave(file.path(opt$outdir, "tf_footprint_heatmap.pdf"),
       p_tf_heatmap, width=5, height=10)

write.csv(bd, file.path(opt$outdir, "bindetect_full_results.csv"),
          row.names=FALSE)

message("── TOBIAS footprinting analysis complete ──")
