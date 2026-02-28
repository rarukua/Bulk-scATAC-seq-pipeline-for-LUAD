# ATAC-seq QC Metrics — What They Mean and Why They Matter

## TSS Enrichment Score
**What it is:** Fold enrichment of ATAC signal at transcription start sites
compared to flanking regions (2000bp window around TSS vs 100bp distal flanks).

**Why it matters:** TSSs are universally accessible in active cells because RNA
Pol II needs to bind there. High TSS enrichment confirms Tn5 successfully
accessed open chromatin rather than random genomic DNA.

**Thresholds:**
- < 4: Failed library — do not proceed with analysis
- 4–6: Poor quality — use with caution, flag in results
- 6–10: Good quality
- > 10: Excellent quality

**Computed by:** ataqv (this pipeline) or deepTools plotProfile

---

## FRiP Score (Fraction of Reads in Peaks)
**What it is:** Proportion of total aligned reads that overlap called peaks.

**Why it matters:** Low FRiP means most reads are in "background" regions —
either the library is noisy or peak calling was too stringent.

**Thresholds (ENCODE ATAC guidelines):**
- < 0.10: Concerning — check library prep and alignment rate
- 0.10–0.20: Acceptable
- > 0.20: Good
- > 0.40: Excellent

**Note:** FRiP is sensitive to peak caller parameters. Relaxing MACS2 `-q`
threshold will increase peak count and FRiP but reduce specificity.

---

## Fragment Size Distribution
**What it is:** Histogram of insert sizes from paired-end sequencing.

**Why it matters:** ATAC-seq should show a characteristic ladder pattern
reflecting nucleosome organisation:
```
< 200bp:   Nucleosome-free regions (NFR) — open chromatin ← most important
200-400bp: Mononucleosomal fragments (1 nucleosome)
400-600bp: Dinucleosomal fragments (2 nucleosomes)
600-800bp: Trinucleosomal (3 nucleosomes)
```
A library without the NFR peak indicates failed Tn5 insertion or poor
open chromatin accessibility. The relative height of NFR vs mono-nucleosomal
peak indicates library quality.

**Enrichment score:** NFR/(mono-nucleosomal) ratio; should be > 1.

---

## Mitochondrial Read Fraction
**What it is:** Percentage of aligned reads mapping to chrM.

**Why it matters:** Mitochondrial DNA is not chromatinised (no nucleosomes),
so Tn5 inserts into it very efficiently. High mito fraction = less nuclear
open chromatin captured.

**Thresholds:**
- < 20%: Good
- 20–50%: Typical for most ATAC protocols
- > 80%: Library failure — likely insufficient nuclear input

**Note:** Mito fraction is removed BEFORE peak calling in this pipeline.
Report it for QC only — it does not affect downstream analysis quality
as long as nuclear reads are sufficient.

---

## Duplication Rate
**What it is:** Percentage of reads flagged as PCR duplicates by Picard.

**Why it matters:** High duplication = low library complexity = fewer unique
Tn5 insertion events captured. Heavily duplicated libraries underrepresent
rare accessible sites.

**Thresholds:**
- < 20%: Excellent library complexity
- 20–50%: Typical
- > 70%: Poor — consider re-sequencing at higher depth or re-library

**Note:** Unlike RNA-seq, some ATAC duplicates are genuine (same Tn5
insertion site in different cells). The standard approach (Picard MarkDuplicates)
may slightly over-remove true signal. Single-cell ATAC tools (Cell Ranger, ArchR)
handle this differently — deduplication is per-barcode.

---

## Alignment Rate
**What it is:** Percentage of trimmed reads that aligned to the reference genome.

**Expected:** > 95% for good-quality human ATAC data with Bowtie2 `--very-sensitive`.

**Low alignment causes:**
- Wrong genome (mouse ATAC aligned to human)
- Heavy adapter contamination (check Trim Galore log)
- Library degradation (check fragment size distribution)
- Contamination (FASTQ screen against common contaminants)
