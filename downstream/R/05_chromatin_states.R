# downstream/R/05_chromatin_states.R
# Chromatin state annotation of ATAC-seq peaks
# ─────────────────────────────────────────────────────────────
# chromHMM learns a Hidden Markov Model from multiple histone
# mark ChIP-seq tracks to define chromatin states:
#   State 1 (TssA):        Active TSS — H3K4me3 + H3K27ac
#   State 2 (TssAFlnk):    Flanking active TSS
#   State 3 (TxFlnk):      Transcr. at gene 5' and 3'
#   State 4 (Tx):          Strong transcription — H3K36me3
#   State 5 (TxWk):        Weak transcription
#   State 6 (EnhG):        Genic enhancers — H3K4me1 + H3K36me3
#   State 7 (Enh):         Enhancers — H3K4me1 + H3K27ac
#   State 8 (ZNF/Rpts):    ZNF genes & repeats
#   State 9 (Het):         Heterochromatin — H3K9me3
#   State 10 (TssBiv):     Bivalent/Poised TSS — H3K4me3 + H3K27me3
#   State 11 (BivFlnk):    Flanking bivalent TSS/Enh
#   State 12 (EnhBiv):     Bivalent Enhancer
#   State 13 (ReprPC):     Repressed PolyComb — H3K27me3
#   State 14 (ReprPCWk):   Weak Repressed PolyComb
#   State 15 (Quies):      Quiescent/Low signal
#
# Here we use the ENCODE/Roadmap pre-computed 15-state model for
# lung tissue (E096: Lung) and overlap our DA ATAC peaks with it.
#
# Key question: Are tumor-specific open chromatin sites enriched at
# ENHANCERS vs PROMOTERS? This distinguishes transcriptional
# reprogramming (enhancer rewiring) from direct TSS activation.
# ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
  library(optparse)
})

opt_list <- list(
  make_option("--da_peaks",      type="character", help="DA peaks CSV from DiffBind"),
  make_option("--chromhmm_bed",  type="character",
              help="chromHMM BED file (e.g. E096_15_coreMarks_hg38lift_mnemonics.bed.gz)"),
  make_option("--consensus_peaks",type="character", help="Consensus peaks BED"),
  make_option("--outdir",        type="character", default="results/chromatin_states")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── chromHMM state colour palette (standard Roadmap colours) ──────────────
state_colours <- c(
  "TssA"     = "#FF0000", "TssAFlnk" = "#FF6E00", "TxFlnk"   = "#32CD32",
  "Tx"       = "#008000", "TxWk"     = "#006400", "EnhG"     = "#C2E105",
  "Enh"      = "#FFFC04", "ZNF/Rpts" = "#66CDAA", "Het"      = "#8A91D0",
  "TssBiv"   = "#CD5C5C", "BivFlnk"  = "#E9967A", "EnhBiv"   = "#BDB76B",
  "ReprPC"   = "#808080", "ReprPCWk" = "#C0C0C0", "Quies"    = "#FFFFFF"
)

# ── State groupings for summary analysis ─────────────────────────────────
state_groups <- list(
  "Active Promoter"  = c("TssA", "TssAFlnk"),
  "Enhancer"         = c("Enh", "EnhG", "EnhBiv"),
  "Transcribed"      = c("Tx", "TxFlnk", "TxWk"),
  "Bivalent/Poised"  = c("TssBiv", "BivFlnk"),
  "Repressed"        = c("ReprPC", "ReprPCWk"),
  "Heterochromatin"  = c("Het", "ZNF/Rpts"),
  "Quiescent"        = "Quies"
)

# ── Load data ──────────────────────────────────────────────────────────────
message("── Loading chromHMM annotations ──")

# Load Roadmap chromHMM BED (lung E096, hg38 liftover)
# Download: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/
chromhmm <- import(opt$chromhmm_bed, format="BED")
chromhmm$state <- gsub("^[0-9]+_", "", chromhmm$name)   # strip state number prefix

# Load DA peaks
da <- read.csv(opt$da_peaks) %>%
  filter(!is.na(seqnames)) %>%
  mutate(direction = case_when(
    FDR < 0.05 & Fold >  1 ~ "Tumor_open",
    FDR < 0.05 & Fold < -1 ~ "Normal_open",
    TRUE                   ~ "Background"
  ))

da_gr <- makeGRangesFromDataFrame(da, keep.extra.columns=TRUE)

# Also load all consensus peaks as background
consensus <- import(opt$consensus_peaks, format="BED")
consensus$direction <- "Background"

# ── Overlap ATAC peaks with chromHMM states ────────────────────────────────
message("── Overlapping peaks with chromHMM states ──")

annotate_with_chromhmm <- function(peaks_gr, chromhmm_gr) {
  hits   <- findOverlaps(peaks_gr, chromhmm_gr, select="first")
  states <- chromhmm_gr$state[hits]
  states[is.na(states)] <- "Quies"
  return(states)
}

da_gr$chromhmm_state <- annotate_with_chromhmm(da_gr, chromhmm)
consensus$chromhmm_state <- annotate_with_chromhmm(consensus, chromhmm)

# ── Plot 1: State distribution by DA peak direction ───────────────────────
state_df <- as.data.frame(da_gr) %>%
  bind_rows(as.data.frame(consensus) %>%
              mutate(direction="Background")) %>%
  count(direction, chromhmm_state) %>%
  group_by(direction) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    direction = factor(direction,
                       levels=c("Tumor_open","Normal_open","Background")),
    state_group = case_when(
      chromhmm_state %in% state_groups[["Active Promoter"]]  ~ "Active Promoter",
      chromhmm_state %in% state_groups[["Enhancer"]]         ~ "Enhancer",
      chromhmm_state %in% state_groups[["Transcribed"]]      ~ "Transcribed",
      chromhmm_state %in% state_groups[["Bivalent/Poised"]]  ~ "Bivalent/Poised",
      chromhmm_state %in% state_groups[["Repressed"]]        ~ "Repressed",
      chromhmm_state %in% state_groups[["Heterochromatin"]]  ~ "Heterochromatin",
      TRUE                                                    ~ "Quiescent"
    )
  )

# Grouped bar chart by state
p_states <- ggplot(state_df,
  aes(x=direction, y=pct, fill=chromhmm_state)) +
  geom_col(position="stack", width=0.7) +
  scale_fill_manual(values=state_colours, na.value="grey80") +
  labs(
    title    = "Chromatin State Composition of DA ATAC Peaks (LUAD)",
    subtitle = "Roadmap Epigenomics 15-state model — Lung tissue (E096)",
    x        = NULL,
    y        = "Percentage of peaks (%)",
    fill     = "chromHMM state",
    caption  = "Tumor-open peaks enriched at enhancers = transcriptional rewiring\nTumor-open peaks at TssA = direct promoter activation"
  ) +
  theme_bw(base_size=13) +
  theme(legend.text=element_text(size=8))

ggsave(file.path(opt$outdir, "chromhmm_state_barplot.pdf"),
       p_states, width=9, height=6)

# ── Plot 2: Enrichment of each state in Tumor-open vs Background ─────────
# Fisher's test for enrichment of each state in tumor-specific peaks
message("── Computing state enrichments ──")

state_enrichment <- lapply(unique(state_df$chromhmm_state), function(s) {
  tumor_in  <- sum(as.data.frame(da_gr)$chromhmm_state == s &
                     as.data.frame(da_gr)$direction == "Tumor_open")
  tumor_tot <- sum(as.data.frame(da_gr)$direction == "Tumor_open")
  bg_in     <- sum(as.data.frame(consensus)$chromhmm_state == s)
  bg_tot    <- length(consensus)

  ft <- fisher.test(matrix(c(tumor_in, tumor_tot - tumor_in,
                              bg_in,    bg_tot    - bg_in), nrow=2))
  data.frame(
    state   = s,
    odds_ratio = ft$estimate,
    pvalue     = ft$p.value,
    pct_tumor  = tumor_in / max(tumor_tot, 1) * 100,
    pct_bg     = bg_in    / max(bg_tot, 1) * 100
  )
}) %>% bind_rows() %>%
  mutate(padj = p.adjust(pvalue, method="BH"),
         log2_or = log2(pmax(odds_ratio, 0.01)),
         significant = padj < 0.05)

p_enrichment <- ggplot(state_enrichment,
  aes(x=reorder(state, log2_or), y=log2_or,
      fill=significant)) +
  geom_col(width=0.7) +
  geom_hline(yintercept=0, linewidth=0.4, color="black") +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="#D7191C", "FALSE"="grey70")) +
  labs(
    title    = "chromHMM State Enrichment in Tumor-Specific Open Peaks",
    subtitle = "Positive = enriched in tumor-open peaks vs all consensus peaks",
    x        = "chromHMM state",
    y        = "log2(Odds Ratio)",
    fill     = "FDR < 0.05"
  ) +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "chromhmm_enrichment_lollipop.pdf"),
       p_enrichment, width=7, height=6)

# ── Plot 3: Bivalent → Active transition (key LUAD biology) ──────────────
# Bivalent chromatin (H3K4me3 + H3K27me3) in normal tissue marks
# developmental genes held in a "poised" state.
# In cancer, these regions often lose H3K27me3 (polycomb repression)
# and become fully active — a key epigenetic driver of oncogenesis.
# Peaks that are BIVALENT in normal but OPEN in tumor capture this.

bivalent_in_normal <- as.data.frame(da_gr) %>%
  filter(direction == "Tumor_open",
         chromhmm_state %in% c("TssBiv","BivFlnk","EnhBiv"))

message(sprintf("  Peaks bivalent in lung but tumor-open: %d", nrow(bivalent_in_normal)))
message(sprintf("  These may represent polycomb derepression events"))

write.csv(bivalent_in_normal,
          file.path(opt$outdir, "bivalent_to_active_peaks.csv"),
          row.names=FALSE)

# ── Summary table ─────────────────────────────────────────────────────────
write.csv(state_enrichment,
          file.path(opt$outdir, "chromhmm_state_enrichment.csv"),
          row.names=FALSE)

message("── Chromatin state analysis complete ──")
