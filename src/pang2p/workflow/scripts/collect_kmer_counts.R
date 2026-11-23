# =============================================================================
# STEP 1 OF 2 — Collect k-mer counts across all accessions
# =============================================================================
#
# PURPOSE:
#   Read per-accession k-mer mapping files in parallel and produce two wide
#   TSV files — one for mapped_kmers, one for count — both in the format:
#       v.id | accession_1 | accession_2 | ... | accession_N
#   Rows are variants; every variant from the pan-genome reference appears in
#   every accession (missing mappings are filled with 0).
#
# INPUTS:
#   1. Pan_count_kmers.tsv
#        Two-column tab-separated file: variant_id, total_kmers
#        (produced by the pan k-mer counting step)
#
#   2. /path/to/04.genotyping/<ACCESSION>/count.tsv.gz
#        One file per accession, tab-separated, no header, three columns:
#            v.id  mapped_kmers  count
#        where:
#            mapped_kmers = number of distinct k-mers from this variant
#                           that were found in this accession's reads
#            count        = total number of read k-mers that matched
#                           (can exceed mapped_kmers if k-mers hit multiple times)
#
# OUTPUTS (written to OUTPUT_DIR):
#   wide_mapped_kmers.tsv  — wide matrix of mapped_kmers values
#   wide_count.tsv         — wide matrix of count values
#
# HOW TO RUN:
#   Rscript step1_collect_kmer_counts.R
#   Adjust the three path variables in the CONFIG section below as needed.
# =============================================================================

options(datatable.print.class = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
})

options(
  datatable.print.class = FALSE,
  scipen = 20
)

option_list <- list(
  make_option(
    c("-i", "--indir"),
    type = "character",
    help = "Directory with count.tsv.gz files (one per accession)"
  ),
  make_option(
    c("-p", "--pancounts"),
    type = "character",
    help = "File with pan kmers and their total counts (Pan_count_kmers.tsv)"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    help = "Output directory for wide TSV files"
  ),
  make_option(
    c("-c", "--cores"),
    type = "integer",
    default = 4,
    help = "Number of parallel workers (default: 4)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# if all three options are not given, print usage and exit
if (is.null(opt$indir) && is.null(opt$outdir) && is.null(opt$pancounts)) {
  print_help(OptionParser(option_list = option_list))
  quit()
}

# =============================================================================
# CONFIG — edit these paths before running
# =============================================================================

PAN_KMER_FILE  <- opt$pancounts  # path to Pan_count_kmers.tsv
GENOTYPING_DIR <- opt$indir      # path to 04.genotyping/ directory
OUTPUT_DIR     <- opt$outdir      # path to output directory for wide TSVs
CORES          <- opt$cores      # number of parallel workers

# =============================================================================
# STEP 1A: Load pan-genome variant reference
# =============================================================================
# gtPAV is the master list of all variants (v.id) and how many k-mers each
# variant was represented by in the pan-genome. This is used as the left-hand
# side of the join so that every variant appears for every accession, even if
# that accession had zero mappings to it.

cat("Loading pan k-mer reference...\n")
gtPAV <- fread(PAN_KMER_FILE, sep = "\t", header = FALSE)
setnames(gtPAV, c("v.id", "total_kmers"))
cat("  Loaded", nrow(gtPAV), "variants.\n")

# =============================================================================
# STEP 1B: Discover per-accession mapping files
# =============================================================================
# Files live at: <GENOTYPING_DIR>/<ACCESSION>/count.tsv.gz
# The accession name is the immediate parent directory of the file.

cat("Discovering accession files under", GENOTYPING_DIR, "...\n")
files <- fread(
  paste("find", GENOTYPING_DIR, "| grep 'count.tsv.gz$'"),
  header = FALSE, col.names = "file"
)
files[, accession := basename(dirname(file))]
cat("  Found", nrow(files), "accession files.\n")
cat("  Example accessions:", paste(head(files$accession, 5), collapse = ", "), "\n")

# =============================================================================
# STEP 1C: Read all files in parallel
# =============================================================================
# For each accession:
#   1. Decompress and read the three-column TSV (v.id, mapped_kmers, count).
#   2. Left-join against gtPAV so ALL variants appear, even unobserved ones.
#   3. Fill NA (= variant not seen in this accession) with 0 for both value cols.
#
# Returns a long-format data.table:
#   v.id | total_kmers | mapped_kmers | count | accession

cat("Reading files in parallel (cores =", CORES, ")...\n")

long_data <- rbindlist(
  mclapply(mc.cores = CORES, 1:nrow(files), function(i) {
    
    acc <- files[i, accession]
    f   <- files[i, file]
    
    tryCatch({
      
      z <- fread(cmd = paste("gzip -kdc", f), header = FALSE, sep = "\t")
      setnames(z, c("v.id", "mapped_kmers", "count"))
      z[, accession := acc]
      
      # Left-join: every row in gtPAV is kept; z rows that match fill in values
      merged <- z[gtPAV, on = "v.id"]
      merged[, accession    := acc]
      merged[is.na(mapped_kmers), mapped_kmers := 0L]
      merged[is.na(count),        count        := 0L]
      
      merged[, .(v.id, accession, mapped_kmers, count)]
      
    }, error = function(e) {
      message("  WARNING: skipping accession '", acc, "' — ", e$message)
      NULL
    })
  }),
  use.names = TRUE
)

# Report any failures
n_expected <- nrow(files)
n_got      <- uniqueN(long_data$accession)
if (n_got < n_expected) {
  missing <- setdiff(files$accession, unique(long_data$accession))
  message("WARNING: ", length(missing), " accession(s) failed and were skipped: ",
          paste(missing, collapse = ", "))
}
cat("  Successfully read", n_got, "of", n_expected, "accessions.\n")
cat("  Long table rows:", nrow(long_data), "\n")

# =============================================================================
# STEP 1D: Cast to wide format and write TSVs
# =============================================================================
# Two separate wide tables, one per value column.
# Shape: rows = variants (v.id), columns = accessions.

cat("Casting to wide format...\n")

# --- mapped_kmers table ---
wide_mapped <- dcast(long_data, v.id ~ accession, value.var = "mapped_kmers")
out_mapped  <- file.path(OUTPUT_DIR, "overall_mapped_kmers.tsv.gz")
fwrite(wide_mapped, file = out_mapped, sep = "\t", quote = FALSE)
cat("  Written:", out_mapped,
    "(", nrow(wide_mapped), "variants x", ncol(wide_mapped) - 1, "accessions )\n")

# --- count table ---
wide_count <- dcast(long_data, v.id ~ accession, value.var = "count")
out_count  <- file.path(OUTPUT_DIR, "overall_counts.tsv.gz")
fwrite(wide_count, file = out_count, sep = "\t", quote = FALSE)
cat("  Written:", out_count,
    "(", nrow(wide_count), "variants x", ncol(wide_count) - 1, "accessions )\n")

cat("\nStep 1 complete. Pass both TSV files to step2_normalize_genotype.R\n")