# =============================================================================
# STEP 2 OF 2 — Normalize k-mer counts and export genotype table
# =============================================================================
#
# PURPOSE:
#   Read the two wide TSV files produced by step1_collect_kmer_counts.R,
#   normalize mapped_kmers counts to reads-per-billion (RPB), filter variants
#   by coefficient of variation (CV), scale values to 0–2, and write the
#   final genotype matrix as a TSV.
#
# INPUTS (produced by step1_collect_kmer_counts.R):
#   wide_mapped_kmers.tsv  — wide matrix: rows = variants, cols = accessions
#                            values = raw mapped_kmers counts
#   wide_count.tsv         — same shape but values = raw count
#                            (available for optional downstream use)
#
# OUTPUTS:
#   Pan_WGS_Genotype_Normalized.tsv
#       Wide matrix: rows = accessions, columns = variant IDs
#       Values are RPB-normalized, CV-filtered, and 0–2 scaled
#
# NORMALIZATION LOGIC (applied to mapped_kmers):
#   1. RPB normalization: per accession, divide each variant's mapped_kmers
#      by the column total, then multiply by 1e9. This corrects for differences
#      in sequencing depth between accessions.
#   2. CV filter: keep only variants where (sd / mean) >= CV_CUTOFF across
#      all accessions. Low-CV variants are nearly constant and carry little
#      genotypic signal.
#   3. 0–2 scaling: for each retained variant, linearly rescale values so the
#      minimum becomes 0 and the maximum becomes 2. This matches diploid dosage
#      encoding used by GWAS tools (0 = ref, 1 = het, 2 = alt).
#
# HOW TO RUN:
#   Rscript step2_normalize_genotype.R
#   Adjust the CONFIG section below as needed.
# =============================================================================

options(datatable.print.class = FALSE)
options(scipen = 999)
options(digits  = 4)

library(data.table)

# =============================================================================
# CONFIG — edit these paths/parameters before running
# =============================================================================

MAPPED_KMERS_FILE <- "wide_mapped_kmers.tsv"   # output from step 1
COUNT_FILE        <- "wide_count.tsv"           # output from step 1 (loaded but not used in normalization; available if needed)
OUTPUT_FILE       <- "Pan_WGS_Genotype_Normalized.tsv"
CV_CUTOFF         <- 0.3    # minimum coefficient of variation to retain a variant

# =============================================================================
# STEP 2A: Load wide tables
# =============================================================================

cat("Loading mapped_kmers table...\n")
w <- fread(MAPPED_KMERS_FILE, sep = "\t")
cat("  Shape:", nrow(w), "variants x", ncol(w) - 1, "accessions\n")

# Identify accession columns (everything except v.id)
acc_cols <- setdiff(names(w), "v.id")
cat("  Accessions:", length(acc_cols), "\n")

# Load count table — available for reference / downstream use
cat("Loading count table (for reference)...\n")
w_count <- fread(COUNT_FILE, sep = "\t")
cat("  Shape:", nrow(w_count), "variants x", ncol(w_count) - 1, "accessions\n")

# =============================================================================
# STEP 2B: RPB normalization of mapped_kmers
# =============================================================================
# For each accession column:
#   normalized_value = (mapped_kmers / sum_of_column) * 1,000,000,000
#
# This converts raw counts to "reads per billion", making values comparable
# across accessions regardless of sequencing depth.

cat("Applying RPB normalization (reads per billion)...\n")
for (col in acc_cols) {
  col_sum <- sum(w[[col]], na.rm = TRUE)
  if (col_sum == 0) {
    warning("  Accession '", col, "' has zero total mapped kmers — column left as 0.")
  } else {
    w[, (col) := w[[col]] / col_sum * 1e9]
  }
}
cat("  Normalization complete.\n")

# =============================================================================
# STEP 2C: Compute per-variant statistics
# =============================================================================
# For each variant (row), compute:
#   mean   — mean RPB across all accessions
#   sd     — standard deviation across accessions
#   cv     — coefficient of variation = sd / mean
#            measures relative variability; higher CV = more discriminating

cat("Computing per-variant mean, SD, and CV...\n")
w[, mean_rpb := rowMeans(.SD, na.rm = TRUE), .SDcols = acc_cols]
w[, sd_rpb   := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = acc_cols]
w[, cv       := sd_rpb / mean_rpb]

cat("  CV distribution (quantiles):\n")
print(quantile(w$cv, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1), na.rm = TRUE))

# =============================================================================
# STEP 2D: Parse chromosome and position from variant IDs
# =============================================================================
# Variant IDs have the format: <source>|<chr>:<start>-<end>
# e.g. "Elliot_HAP1|Chr01:10004352-10004474"
# Extract chr, start, end for position-aware output.

cat("Parsing chromosomal coordinates from variant IDs...\n")
w[, c("source", "coords") := tstrsplit(v.id, "|", fixed = TRUE)]
w[, c("chr", "pos_range") := tstrsplit(coords, ":", fixed = TRUE)]
w[, c("start", "end")     := tstrsplit(pos_range, "-", fixed = TRUE)]
w[, start := as.integer(start)]
w[, c("source", "coords", "pos_range") := NULL]

# =============================================================================
# STEP 2E: Resolve duplicate chr:start positions
# =============================================================================
# Multiple variants may share the same chr:start (overlapping k-mer windows).
# Assign unique positions by incrementing duplicates by 1 until position is free.
# This ensures each variant gets a unique genomic coordinate for GWAS tools.

cat("Resolving duplicate chromosomal positions...\n")
resolve_duplicates <- function(chr_vec, pos_vec) {
  seen         <- new.env(hash = TRUE)
  resolved_pos <- integer(length(pos_vec))
  for (i in seq_along(pos_vec)) {
    chr <- chr_vec[i]
    pos <- pos_vec[i]
    key <- paste(chr, pos, sep = ":")
    while (!is.null(seen[[key]])) {
      pos <- pos + 1L
      key <- paste(chr, pos, sep = ":")
    }
    seen[[key]]    <- TRUE
    resolved_pos[i] <- pos
  }
  resolved_pos
}

w[, Position   := resolve_duplicates(chr, start)]
w[, Chromosome := as.integer(factor(chr, levels = sort(unique(chr))))]

cat("  Marker count per chromosome (before CV filter):\n")
print(w[, .N, by = .(Chromosome, chr)][order(Chromosome)])

# =============================================================================
# STEP 2F: Filter variants by CV
# =============================================================================
# Retain only variants with CV >= CV_CUTOFF.
# Variants with very low CV are nearly invariant across accessions and add
# noise without genotypic signal. CV >= 0.3 is a commonly used threshold.

cat("Filtering variants with CV >=", CV_CUTOFF, "...\n")
n_before <- nrow(w)
gw <- w[cv >= CV_CUTOFF]
n_after  <- nrow(gw)
cat("  Retained", n_after, "of", n_before, "variants (", round(n_after/n_before*100, 1), "%).\n")

cat("  Marker count per chromosome (after CV filter):\n")
print(gw[, .N, by = .(Chromosome, chr)][order(Chromosome)])

# =============================================================================
# STEP 2G: Scale RPB values to 0–2 range per variant
# =============================================================================
# For each variant independently:
#   scaled = 2 * (value - min) / (max - min)
#
# Result is in [0, 2], compatible with diploid dosage genotype encoding:
#   0 = homozygous reference-like (lowest expression/presence)
#   2 = homozygous alternate-like (highest expression/presence)
#
# Work in long format for clean per-variant scaling, then cast back to wide.

cat("Scaling values to 0–2 per variant...\n")

# Drop stat/coordinate columns before melting — keep only v.id + accession cols
gw[, c("mean_rpb", "sd_rpb", "cv", "chr", "start", "end",
       "Chromosome", "Position") := NULL]

long <- melt(
  gw,
  id.vars       = "v.id",
  measure.vars  = acc_cols,
  variable.name = "accession",
  value.name    = "rpb",
  variable.factor = FALSE
)

long[, scaled := {
  mn <- min(rpb, na.rm = TRUE)
  mx <- max(rpb, na.rm = TRUE)
  if (mx == mn) rep(0, .N) else 2 * (rpb - mn) / (mx - mn)
}, by = v.id]

cat("  Scaling complete.\n")

# =============================================================================
# STEP 2H: Cast to final wide genotype matrix and write TSV
# =============================================================================
# Final shape: rows = accessions, columns = variant IDs (+ accession label col)
# Values are 0–2 scaled, RPB-normalized pseudo-dosage genotypes.

cat("Casting to final wide genotype matrix...\n")
genotype_matrix <- dcast(long, accession ~ v.id, value.var = "scaled")

cat("  Final matrix:", nrow(genotype_matrix), "accessions x",
    ncol(genotype_matrix) - 1, "variants\n")

fwrite(genotype_matrix, file = OUTPUT_FILE, sep = "\t", quote = FALSE)

cat("\n=============================================================\n")
cat("Genotype table written to:", OUTPUT_FILE, "\n")
cat("Dimensions            :", nrow(genotype_matrix), "accessions x",
    ncol(genotype_matrix) - 1, "variants\n")
cat("CV cutoff applied     :", CV_CUTOFF, "\n")
cat("Value range           : 0 to 2 (RPB-normalized, per-variant scaled)\n")
cat("=============================================================\n")