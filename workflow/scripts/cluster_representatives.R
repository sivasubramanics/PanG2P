#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(igraph)
})

options(
  datatable.print.class = FALSE,
  scipen = 20
)

option_list <- list(
  make_option(
    c("-c", "--clusters"),
    type = "character",
    help = "MMseqs cluster TSV",
    required = TRUE
  ),
  make_option(
    c("-f", "--fai"),
    type = "character",
    help = "FASTA index (.fai)",
    required = TRUE
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    help = "Output directory",
    required = TRUE
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

##
## Read cluster graph
##

mq <- fread(opt$clusters, sep = "\t", header = FALSE)
setnames(mq, c("seq1", "seq2"))

g <- graph_from_edgelist(
  as.matrix(mq[, .(seq1, seq2)]),
  directed = FALSE
)

##
## Read sequence metadata
##

info <- fread(
  opt$fai,
  select = 1:2,
  col.names = c("seq", "seq.length")
)

info[, genotype := sub("\\|.*", "", seq)]
info[, chr := sub(":.*", "", sub("^.*\\|", "", seq))]
info[, start := as.integer(sub("-.*", "", sub("^.*:", "", seq)))]
info[, end := as.integer(sub(".*-", "", seq))]

setkey(info, seq)

##
## Add singleton sequences
##

g <- g + vertex(info[!seq %in% c(mq$seq1, mq$seq2)]$seq)

##
## Cluster membership
##

mem <- data.table(
  seq = V(g)$name,
  components = paste0("cluster_", clusters(g)$membership)
)

mem <- info[mem, on = "seq"]

fwrite(
  mem,
  file.path(opt$outdir, "pan_cluster_all_member.txt"),
  sep = "\t",
  quote = FALSE
)

##
## Representative selection
##

setorder(mem, -seq.length)

mem[, cluster.idx := seq_len(.N), by = components]
mem[, rep := cluster.idx == 1]

cl <- mem[, .(size = .N), keyby = components]

cl <- mem[
  rep == TRUE,
  .(
    components,
    rep = seq,
    rep.length = seq.length
  )
][cl, on = "components"]

cl <- mem[
  ,
  .N,
  keyby = .(components, genotype)
][
  ,
  .(ngenotype = .N),
  keyby = components
][cl, on = "components"]

##
## Chromosome assignment
##

x <- mem[
  ,
  .N,
  keyby = .(chr, genotype, components)
][
  order(-N)
][
  ,
  c("idx", "pchr") := list(.GRP, N / sum(N)),
  by = .(genotype, components)
]

y <- x[!duplicated(idx), .(chr, pchr, idx, components, genotype)]

y <- x[
  chr != "chrUn",
  .(nchr = .N),
  keyby = idx
][
  y,
  on = "idx"
]

w <- mem[
  y[, .(chr, genotype, components)],
  on = c("chr", "genotype", "components")
]

yy <- y[
  w[
    ,
    .(
      start = min(start),
      end = max(end)
    ),
    keyby = .(components, genotype)
  ],
  on = c("components", "genotype")
]

yy[, length := end - start + 1]

yy <- mem[
  ,
  .(nseq = .N),
  keyby = .(genotype, components)
][
  yy,
  on = c("genotype", "components")
]

ww <- yy[
  ,
  .(
    max.length = max(length),
    min.length = min(length),
    mean.length = mean(length),
    sd.length = sd(length)
  ),
  keyby = components
]

ww[, cv.length := sd.length / mean.length]

ww <- y[
  chr != "chrUn",
  .N,
  keyby = .(components, chr)
][
  ,
  .(nchr = .N),
  keyby = components
][
  ww,
  on = "components"
]

cl <- ww[cl, on = "components"]

cl <- mem[
  ,
  .(
    rep = seq,
    rep.genotype = genotype
  )
][
  cl,
  on = "rep"
]

setkey(cl, components)

cl <- cl[order(rep.genotype)]

##
## BED output
##

cl[, c("i1", "i2") := tstrsplit(rep, ":", fixed = TRUE)]
cl[, c("start", "end") := tstrsplit(i2, "-", fixed = TRUE)]
cl[, c("gn", "chr") := tstrsplit(i1, "|", fixed = TRUE)]
cl[, c("i1", "i2", "gn") := NULL]

sc <- cl[
  ,
  .(
    chr,
    start,
    end,
    rep,
    rep.genotype
  )
]

##
## Write outputs
##

fwrite(
  cl,
  file.path(opt$outdir, "pan_cluster_mmseq_rep.txt"),
  sep = "\t",
  quote = FALSE
)

fwrite(
  sc,
  file.path(opt$outdir, "pan_cluster_bed_rep.txt"),
  sep = "\t",
  quote = FALSE
)
