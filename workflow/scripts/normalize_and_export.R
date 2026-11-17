# =============================================================================
# STEP 2 — Normalize k-mers, filter zero-variance, and export both outputs
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
  make_option(c("-m", "--cm"), type="character", help="Mapped k-mers gzipped file"),
  make_option(c("-c", "--co"), type="character", help="Count gzipped file"),
  make_option(c("-o", "--outprefix"), type="character", help="Output file prefix (path + prefix)"),
  make_option(c("-v", "--cv"), type="numeric", default=0.3, help="CV cutoff threshold")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cm) || is.null(opt$co) || is.null(opt$outprefix)) {
  print_help(OptionParser(option_list = option_list))
  quit()
}

# Ensure the output directory exists if outprefix includes a path
out_dir <- dirname(opt$outprefix)
if (out_dir != "." && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", out_dir, "\n")
}

MAPPED_KMERS_FILE <- opt$cm
COUNT_FILE        <- opt$co
OUT_RDATA         <- paste0(opt$outprefix, ".GAPIT.RData")
OUT_TSV           <- paste0(opt$outprefix, ".genotype.tsv")
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
# Statistics & Coordinate Parsing
# -----------------------------------------------------------------------------

cat("Computing CV and parsing coordinates...\n")

w[, mean_rpb := rowMeans(.SD, na.rm=TRUE), .SDcols=acc_cols]
w[, sd_rpb   := apply(.SD, 1, sd, na.rm=TRUE), .SDcols=acc_cols]
w[, cv       := sd_rpb / mean_rpb]

# Parse coordinates to extract CHR, START, END
w[, c("source","coords") := tstrsplit(v.id, "|", fixed=TRUE)]
w[, c("chr","pos") := tstrsplit(coords, ":", fixed=TRUE)]

# Clean chromosome name to strictly numeric (e.g., "Chr16" -> 16)
w[, chr := as.integer(gsub("[^0-9]", "", chr))]

# Clean 'pos' to remove any trailing suffixes (e.g., _L, _R) before splitting start-end
w[, pos := sub("_.*$", "", pos)]

w[, c("start","end") := tstrsplit(pos, "-", fixed=TRUE)]
w[, start := as.integer(start)]
w[, end   := as.integer(end)]

# Drop rows with invalid/missing coordinates
w <- w[!is.na(chr) & !is.na(start) & !is.na(end)]

# -----------------------------------------------------------------------------
# FILTER
# -----------------------------------------------------------------------------

cat("Filtering CV >=", CV_CUTOFF, "...\n")
gw <- w[cv >= CV_CUTOFF]
cat("Variants kept after CV filter:", nrow(gw), "\n")

# Sort variants by chromosome and start position
setorder(gw, chr, start)
variant_order <- gw$v.id

# -----------------------------------------------------------------------------
# SCALE PER VARIANT
# -----------------------------------------------------------------------------

gw[, c("mean_rpb","sd_rpb","cv","source","coords","pos") := NULL]

long <- melt(
  gw,
  id.vars = c("v.id", "chr", "start", "end"),
  measure.vars = acc_cols,
  variable.name = "accession",
  value.name = "rpb"
)

long[, scaled := {
  mn <- min(rpb)
  mx <- max(rpb)
  if (mx == mn) rep(0, .N)
  else 2 * (rpb - mn) / (mx - mn)
}, by = v.id]

long[, scaled := round(scaled, 2)]

# -----------------------------------------------------------------------------
# BUILD WIDE FORMAT & FILTER ZERO-VARIANCE MARKERS
# -----------------------------------------------------------------------------

cat("Formatting genotype matrix and removing zero-variance markers...\n")

wide <- dcast(long, v.id + chr + start + end ~ accession, value.var = "scaled")
setnames(wide, c("v.id", "chr", "start", "end"), c("ID", "CHR", "START", "END"))

# Ensure sample columns are strictly numeric
sample_cols <- setdiff(names(wide), c("ID", "CHR", "START", "END"))
wide[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]

# Remove zero-variance markers (where all samples have the identical score)
wide[, marker_sd := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = sample_cols]
initial_count <- nrow(wide)
wide <- wide[marker_sd > 0]
wide[, marker_sd := NULL]

cat("Removed", initial_count - nrow(wide), "zero-variance markers.\n")
cat("Final markers retained:", nrow(wide), "\n")

# Re-sort to maintain clean order
setkey(wide, ID)
wide <- wide[intersect(variant_order, wide$ID)]

# 1. Write the clean genotype TSV file
cat("Writing genotype TSV to:", OUT_TSV, "\n")
fwrite(wide, OUT_TSV, sep="\t", quote=FALSE)

# -----------------------------------------------------------------------------
# BUILD GAPIT OBJECTS (GD, GM) AND SAVE RDATA
# -----------------------------------------------------------------------------

cat("Building GAPIT objects (GD and GM)...\n")

GM <- data.frame(
  Name = wide$ID,
  Chromosome = wide$CHR,
  Position = wide$START,
  stringsAsFactors = FALSE
)

gt_matrix <- as.matrix(wide[, ..sample_cols])
rownames(gt_matrix) <- wide$ID
GD_transposed <- t(gt_matrix)
storage.mode(GD_transposed) <- "numeric"

GD <- data.frame(
  Taxa = rownames(GD_transposed),
  GD_transposed,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Saving GAPIT.RData to:", OUT_RDATA, "\n")
save(GD, GM, file = OUT_RDATA)

cat("\nDONE\n")
cat("Genotype TSV :", OUT_TSV, "\n")
cat("GAPIT RData  :", OUT_RDATA, "\n")
cat("Markers      :", nrow(GM), "\n")
cat("Taxa (Samples):", nrow(GD), "\n")
