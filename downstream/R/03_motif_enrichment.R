# downstream/R/03_motif_enrichment.R
# Motif enrichment analysis — HOMER output parsing and visualisation
# ─────────────────────────────────────────────────────────────
# HOMER findMotifsGenome is run as a Nextflow process (see modules/).
# This script parses the output and produces publication-quality figures.
#
# HOMER vs TOBIAS — important distinction to articulate in interviews:
#   HOMER:  sequence-based — "is this motif over-represented in open peaks?"
#           Does NOT tell you whether the TF is actually bound.
#           Fast, widely used, good for discovery.
#
#   TOBIAS: signal-based — "is there a footprint (Tn5 protection dip)
#           at this motif, suggesting the TF is occupying its binding site?"
#           Requires high-quality ATAC data. Distinguishes bound from unbound.
#
# Use HOMER first for discovery, then TOBIAS to validate binding.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(stringr)
  library(optparse)
})

opt_list <- list(
  make_option("--homer_tumor",  type="character",
              help="HOMER output dir for tumor-specific peaks"),
  make_option("--homer_normal", type="character",
              help="HOMER output dir for normal-specific peaks"),
  make_option("--outdir",       type="character", default="results/motif_enrichment")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Parse HOMER knownResults.txt ───────────────────────────────────────────
parse_homer <- function(homer_dir, condition) {
  known_file <- file.path(homer_dir, "knownResults.txt")
  if (!file.exists(known_file)) stop(paste("Not found:", known_file))

  df <- read.delim(known_file, stringsAsFactors=FALSE) %>%
    # HOMER column names have spaces — clean them up
    setNames(make.names(colnames(.))) %>%
    transmute(
      motif_name  = Motif.Name,
      tf_family   = str_extract(Motif.Name, "(?<=\\().*(?=\\))"),
      log_pvalue  = as.numeric(gsub("1e-","1e-",Log.P.value)),
      pvalue      = exp(log_pvalue),
      fold_enrich = as.numeric(gsub("%","",X..of.Target.Sequences.with.Motif)) /
                   pmax(as.numeric(gsub("%","",X..of.Background.Sequences.with.Motif)), 0.1),
      pct_target  = as.numeric(gsub("%","",X..of.Target.Sequences.with.Motif)),
      pct_bg      = as.numeric(gsub("%","",X..of.Background.Sequences.with.Motif)),
      condition   = condition
    ) %>%
    mutate(padj = p.adjust(pvalue, method="BH")) %>%
    arrange(pvalue)
  return(df)
}

message("── Parsing HOMER results ──")
tumor_motifs  <- parse_homer(opt$homer_tumor,  "Tumor-specific peaks")
normal_motifs <- parse_homer(opt$homer_normal, "Normal-specific peaks")
all_motifs    <- bind_rows(tumor_motifs, normal_motifs)

# ── LUAD-relevant TF families to highlight ────────────────────────────────
# These are the key transcription factors in LUAD biology:
#   NKX2-1:  master lineage TF, amplified/overexpressed in LUAD
#   FOXA1/2: pioneer factors that open chromatin at lung-specific enhancers
#   AP-1:    JUN/FOS family, oncogenic in LUAD (RAS pathway downstream)
#   TP63:    squamous vs adenocarcinoma lineage marker
#   CEBP:    myeloid/lung lineage, LUAD survival marker
#   NFI:     nuclear factor I, lung development
highlight_tfs <- c("NKX2-1","FOXA1","FOXA2","AP-1","JUNB","FOSL2",
                    "TP63","CEBPA","CEBPB","NFI","KLF","TEAD","ETS")

# ── Plot 1: Dot plot — top enriched motifs per condition ─────────────────
top_per_condition <- all_motifs %>%
  filter(padj < 0.05, fold_enrich > 1.5) %>%
  group_by(condition) %>%
  slice_min(pvalue, n=20) %>%
  ungroup() %>%
  mutate(
    tf_short  = str_extract(motif_name, "^[^/\\(]+") %>% trimws(),
    highlight = str_detect(toupper(tf_short),
                            paste(toupper(highlight_tfs), collapse="|"))
  )

p_dotplot <- ggplot(top_per_condition,
  aes(x=condition,
      y=reorder(tf_short, -log10(pvalue)),
      size=-log10(pvalue),
      color=log2(fold_enrich))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(range=c(2, 9), name="-log10(p-value)") +
  scale_color_gradient2(
    low="grey80", mid="#FDAE61", high="#D7191C",
    midpoint=1, name="log2(Fold\nEnrichment)"
  ) +
  facet_wrap(~condition, scales="free_x") +
  labs(
    title    = "HOMER Motif Enrichment — LUAD Tumor vs Normal Peaks",
    subtitle = "Top 20 enriched motifs per peak set (FDR < 0.05)",
    x        = NULL, y = "TF Motif",
    caption  = "Size = significance; Color = fold enrichment over background"
  ) +
  theme_bw(base_size=11) +
  theme(axis.text.y=element_text(size=8),
        strip.background=element_rect(fill="grey90"))

ggsave(file.path(opt$outdir, "homer_motif_dotplot.pdf"),
       p_dotplot, width=10, height=10)

# ── Plot 2: Comparison scatter — tumor vs normal enrichment ───────────────
# For each motif, plot enrichment in tumor vs normal peaks.
# Motifs in upper-left = specific to tumor; lower-right = specific to normal
motif_wide <- all_motifs %>%
  select(motif_name, tf_short=motif_name, condition, log_pvalue) %>%
  mutate(tf_short=str_extract(motif_name, "^[^/\\(]+") %>% trimws()) %>%
  pivot_wider(id_cols=tf_short,
              names_from=condition,
              values_from=log_pvalue,
              values_fill=0) %>%
  setNames(c("tf","tumor_logp","normal_logp")) %>%
  mutate(
    tumor_specific  = tumor_logp < -5 & tumor_logp < normal_logp * 0.5,
    normal_specific = normal_logp < -5 & normal_logp < tumor_logp * 0.5,
    category = case_when(
      tumor_specific  ~ "Tumor-specific",
      normal_specific ~ "Normal-specific",
      TRUE            ~ "Shared"
    )
  )

p_comparison <- ggplot(motif_wide,
  aes(x=-normal_logp, y=-tumor_logp, color=category)) +
  geom_point(alpha=0.5, size=1.2) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_label_repel(
    data=motif_wide %>%
      filter(category != "Shared") %>%
      arrange(tumor_logp + normal_logp) %>%
      head(20),
    aes(label=tf), size=2.5, max.overlaps=15
  ) +
  scale_color_manual(values=c(
    "Tumor-specific"  ="#D7191C",
    "Normal-specific" ="#2C7BB6",
    "Shared"          ="grey70"
  )) +
  labs(
    title    = "Motif Enrichment: Tumor vs Normal Peaks",
    subtitle = "Points above diagonal = tumor-enriched motifs",
    x        = "-log10(p) in Normal-specific peaks",
    y        = "-log10(p) in Tumor-specific peaks",
    color    = NULL
  ) +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "homer_tumor_vs_normal_scatter.pdf"),
       p_comparison, width=7, height=6)

write.csv(all_motifs, file.path(opt$outdir, "homer_all_motifs.csv"),
          row.names=FALSE)

message("── Motif enrichment analysis complete ──")
message("  Next: run 04_tobias_footprinting.R to validate which motifs are BOUND")
