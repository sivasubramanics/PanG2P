options(datatable.print.class = FALSE)
library(data.table)
library(parallel)
options(scipen = 999)


# Read kmer count for each v.ids
fread("/scratch/group/pgenomics/EXTERNAL/pecan/assemblies/Pan_count_kmers.tsv", sep = "\t") -> gtPAV
setnames(gtPAV, c("v.id", "total_kmers"))


# Now read kmer mapping results
fread("find /scratch/group/pgenomics/EXTERNAL/pecan/assemblies/genotype_result | grep 'count.tsv.gz$'",head=F, col.names="file")-> files
files[, accession := sub("_count.tsv.gz$", "", sub("^.*/genotype_result/", "", file))]


# Select accession that have phenotype
# format phenotype data
fread("/scratch/group/pgenomics/EXTERNAL/pecan/assemblies/phenotype.txt") -> ph

# Name the sample/accession name colum as taxa
setnames(ph, "sample", "taxa")

files[accession %in% ph$taxa] -> files
ph[taxa %in%files$accession] -> ph


# define function
read_kmer_count<-function(info, gt, cores=1){
  f <- info
  rbindlist(mclapply(mc.cores=cores, 1:nrow(files), function(i){
    n <- files[i, file]
    fread(cmd = paste("zcat", n), header = FALSE, sep = " ") ->  z 
    setnames(z, c("v.id", "mapped_kmers", "count"))
    z[, accession := files[i, accession]]
    z[gt, on="v.id"][, .(v.id, total_kmers, mapped_kmers, count, accession)] -> c
    c[is.na(mapped_kmers), accession := files[i, accession]]
    c[is.na(mapped_kmers), mapped_kmers := 0]
    c[is.na(count), count := 0]
  })) 
}

# Now read
read_kmer_count(info=files, gt=gtPAV, cores=20) -> PAV_gt


# change into wide table of mapped kmer count
dcast(PAV_gt, v.id + total_kmers ~ accession, value.var='mapped_kmers')->w


saveRDS(w, file="kmer_raw_count_table.Rds")

# Normalize mapped count to billion per sample
options(digits = 2)
setdiff(names(w), c("v.id", "total_kmers")) -> lst

for (i in lst) w[, (i) := w[[i]]/sum(w[[i]])*1000000000]

# Count row mean and st.deviation
w[, mean := rowMeans(.SD), .SDcols = lst] 
w[, sd := sd(.SD), by=v.id, .SDcols = lst]
w[, cutoff := (sd/mean)]

# Make chromsome info and check markers
w[, .(v.id)] -> vid
vid[, c("ac", "chid") := tstrsplit(v.id, "|", fixed=TRUE)][, c("chr","pos") := tstrsplit(chid, ":", fixed=TRUE)]
vid[, c("start", "end") := tstrsplit(pos, "-", fixed=TRUE)][, c("ac","chid", "pos") := c(NULL, NULL, NULL)]
vid[w, on="v.id"] -> w


####Fix dupliate position


w[, variant_id := paste(chr, start, sep = ":")]
w[, start := as.integer(start)]

# Function to resolve duplicates by incrementing positions

resolve_duplicates <- function(chr_vec, pos_vec) {
  seen <- new.env(hash = TRUE)
  resolved_pos <- integer(length(pos_vec))
  
  for (i in seq_along(pos_vec)) {
    chr <- chr_vec[i]
    pos <- pos_vec[i]
    key <- paste(chr, pos, sep = ":")
    
    while (!is.null(seen[[key]])) {
      pos <- pos + 1
      key <- paste(chr, pos, sep = ":")
    }
    
    seen[[key]] <- TRUE
    resolved_pos[i] <- pos
  }
  
  return(resolved_pos)
}


# Apply the function
w[, Position := as.integer(resolve_duplicates(chr, start))]

# modify chrs to numeric 
w[, Chromosome := as.integer(factor(chr, levels = sort(unique(chr))))]

# filter variants that show optimum variation for GWAS
w[cutoff >= 0.3] -> gw

#Check no. of markers per chromosome
w[, .N, by=Chromosome]
gw[, .N, by=Chromosome]


#saveRDS(gw, file="pan_WGS_GT_info_table.Rds")

gw[, c("mean", "sd", "cutoff", "total_kmers", "variant_id", "chr", "start", "end") := c(NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL)]

gw[, .(SNP=v.id, Chromosome, Position)] -> gm

# convert the data.table to long format
melt(gw, id="v.id", measure=setdiff(names(gw), c("v.id", "Chromosome", "Position")), variable.factor=F, variable.name="id", value.name="count")->z

#scale to interval 0:2 to make counts compatible with GAPIT
z[, x := 2*(count - min(count))/(max(count) - min(count)), by=v.id]



#cast to wide matrix for GAPIT: rows are samples, columns are markers
dcast(z, id ~ v.id, value.var='x')->q

#Save
saveRDS(q, file="Pan_WGS_GWAS_GAPIT_value.Rds")
