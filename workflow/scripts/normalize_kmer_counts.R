# =============================================================================
# STEP 2 — Normalize k-mers and export genotype + SNP map (CLEAN VERSION)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

options(datatable.print.class = FALSE)
options(scipen = 999)

# -----------------------------------------------------------------------------
# Args
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-m", "--cm"), type="character"),
  make_option(c("-c", "--co"), type="character"),
  make_option(c("-o", "--outdir"), type="character"),
  make_option(c("-v", "--cv"), type="numeric", default=0.3)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cm) || is.null(opt$co) || is.null(opt$outdir)) {
  print_help(OptionParser(option_list = option_list))
  quit()
}

MAPPED_KMERS_FILE <- opt$cm
COUNT_FILE        <- opt$co
OUTMATRIX       <- paste0(opt$outdir , "/genotype_matrix.tsv.gz")
OUTMAP          <- paste0(opt$outdir , "/genotype_map.tsv.gz")
CV_CUTOFF         <- opt$cv

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

cat("Loading mapped k-mers...\n")
w <- fread(cmd = paste("gzip -kdc", MAPPED_KMERS_FILE), sep="\t")

acc_cols <- setdiff(names(w), "v.id")

cat("Loading counts (reference only)...\n")
w_count <- fread(cmd = paste("gzip -kdc", COUNT_FILE), sep="\t")

# -----------------------------------------------------------------------------
# RPB normalization
# -----------------------------------------------------------------------------

cat("RPB normalization...\n")

for (col in acc_cols) {
  col_sum <- sum(w[[col]], na.rm = TRUE)
  
  if (col_sum > 0) {
    w[, (col) := (get(col) / col_sum) * 1e9]
  } else {
    w[, (col) := 0]
  }
}

# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

cat("Computing CV...\n")

w[, mean_rpb := rowMeans(.SD, na.rm=TRUE), .SDcols=acc_cols]
w[, sd_rpb   := apply(.SD, 1, sd, na.rm=TRUE), .SDcols=acc_cols]
w[, cv       := sd_rpb / mean_rpb]

# -----------------------------------------------------------------------------
# Parse coordinates
# -----------------------------------------------------------------------------

w[, c("source","coords") := tstrsplit(v.id, "|", fixed=TRUE)]
w[, c("chr","pos") := tstrsplit(coords, ":", fixed=TRUE)]
w[, c("start","end") := tstrsplit(pos, "-", fixed=TRUE)]
w[, start := as.integer(start)]

# numeric chromosome index (stable ordering)
chr_levels <- sort(unique(w$chr))
w[, Chromosome := match(chr, chr_levels)]

# -----------------------------------------------------------------------------
# Resolve duplicate positions
# -----------------------------------------------------------------------------

resolve_duplicates <- function(chr_vec, pos_vec) {
  seen <- new.env(hash=TRUE)
  out <- integer(length(pos_vec))
  
  for (i in seq_along(pos_vec)) {
    k <- paste(chr_vec[i], pos_vec[i], sep=":")
    while (!is.null(seen[[k]])) {
      pos_vec[i] <- pos_vec[i] + 1L
      k <- paste(chr_vec[i], pos_vec[i], sep=":")
    }
    seen[[k]] <- TRUE
    out[i] <- pos_vec[i]
  }
  out
}

w[, Position := resolve_duplicates(chr, start)]

# -----------------------------------------------------------------------------
# FILTER (THIS IS THE SINGLE SOURCE OF TRUTH)
# -----------------------------------------------------------------------------

cat("Filtering CV >=", CV_CUTOFF, "...\n")

gw <- w[cv >= CV_CUTOFF]

cat("Variants kept:", nrow(gw), "\n")

# -----------------------------------------------------------------------------
# DEFINE FINAL ORDER (CRITICAL FIX)
# -----------------------------------------------------------------------------

setorder(gw, Chromosome, Position)

variant_order <- gw$v.id

# -----------------------------------------------------------------------------
# BUILD SNP MAP (FROM SAME OBJECT)
# -----------------------------------------------------------------------------

cat("Writing SNP map...\n")

snp_map <- gw[, .(
  snp.id = v.id,
  Chromosome,
  Position
)]

fwrite(snp_map, OUTMAP, sep="\t", quote=FALSE)

# -----------------------------------------------------------------------------
# BUILD LONG FORMAT FOR MATRIX
# -----------------------------------------------------------------------------

gw[, c("mean_rpb","sd_rpb","cv","chr","start","end") := NULL]

long <- melt(
  gw,
  id.vars = "v.id",
  measure.vars = acc_cols,
  variable.name = "accession",
  value.name = "rpb"
)

# -----------------------------------------------------------------------------
# SCALE PER VARIANT
# -----------------------------------------------------------------------------

long[, scaled := {
  mn <- min(rpb)
  mx <- max(rpb)
  if (mx == mn) rep(0, .N)
  else 2 * (rpb - mn) / (mx - mn)
}, by = v.id]

# -----------------------------------------------------------------------------
# CAST MATRIX
# -----------------------------------------------------------------------------

cat("Building genotype matrix...\n")

genotype_matrix <- dcast(long, accession ~ v.id, value.var="scaled")

# -----------------------------------------------------------------------------
# FORCE CORRECT COLUMN ORDER (CRITICAL)
# -----------------------------------------------------------------------------

genotype_matrix <- genotype_matrix[, c("accession", variant_order), with=FALSE]

# -----------------------------------------------------------------------------
# ROUND (FINAL STEP ONLY)
# -----------------------------------------------------------------------------

acc_cols <- setdiff(names(genotype_matrix), "accession")

genotype_matrix[, (acc_cols) := lapply(.SD, round, 2), .SDcols=acc_cols]

# -----------------------------------------------------------------------------
# WRITE OUTPUT
# -----------------------------------------------------------------------------

fwrite(genotype_matrix, OUTMATRIX, sep="\t", quote=FALSE)

cat("\nDONE\n")
cat("Genotype:", OUTMATRIX, "\n")
cat("SNP map :", OUTMAP, "\n")
cat("Markers :", length(variant_order), "\n")