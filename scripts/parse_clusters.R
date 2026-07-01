source("pseudomolecule_construction.R")
library(igraph)
library(xlsx)
library(data.table)
library(parallel)
library(stringi)
library(zoo)
options(datatable.print.class = FALSE)
options(scipen=20) #needed to ensure proper formating of numbers

###
fread("/Users/selva001/projects/work/tamu/pangat/pangat/work/02.clustering/clusters.tsv", sep="\t", header=F) -> mq
setnames(mq, c("seq1", "seq2"))

#analyze graph
graph_from_edgelist(as.matrix(mq[, .(seq1, seq2)]), directed=F) -> g  # igraph 2.0.0.

# make info 
f<- '/Users/selva001/projects/work/tamu/pangat/pangat/work/01.single_copy/single_copy.fa.fai'
fread(f, sel=1:2, col.names=c('seq', 'seq.length'))->info
info[, genotype := sub("\\|.*", "", seq)]
info[, chr := sub(":.*", "", sub("^.*\\|", "", seq))]
info[, start := as.integer(sub("-.*", "", sub("^.*:", "", seq)))]
info[, end := as.integer(sub(".*-", "", seq))]
setkey(info, seq)

g + vertex(info[!c(mq$seq1, mq$seq2)]$seq) -> g

data.table(seq=V(g)$name, components=paste0("cluster_", clusters(g)$membership))->mem
info[mem, on='seq']->mem

# save the members list in each cluster
fwrite(mem, file="pan_cluster_all_member.txt", sep="\t", quote=F)



#the cluster representative is the longest sequence
setorder(mem, -seq.length)
mem[, cluster.idx := 1:.N, by=components]
mem[, rep := cluster.idx == 1]
mem[, .(size=.N), key=components] -> cl
mem[rep == T, .(components, rep=seq, rep.length=seq.length)][cl, on='components']->cl
mem[, .N, key=.(components, genotype)][, .(ngenotype=.N), key=components][cl, on='components']->cl

# assign clusters to chromosomes (majority rule as in TRITEX)
mem[, .N, key=.(chr, genotype, components)][order(-N)][, c("idx", "pchr") := list(.GRP, N/sum(N)), by=.(genotype, components)][]->x

x[!duplicated(idx)][, .(chr, pchr, idx, components, genotype)] -> y
x[chr != 'chrUn', .(nchr=.N), key=idx][y, on='idx']->y
mem[y[, .(chr, genotype, components)], on=c('chr', 'genotype', 'components')]->w
y[w[, .(start=min(start), end=max(end)), key=.(components, genotype)], on=c('components', 'genotype')]->yy
yy[, length := 1 + end - start]
yy[mem[, .(nseq=.N), key=.(genotype, components)], on=c('genotype', 'components')]->yy

yy[, .(max.length = max(length), min.length=min(length), mean.length=mean(length), 
       sd.length=sd(length)), key=.(components)]->ww
ww[, cv.length := sd.length / mean.length]
y[chr != 'chrUn', .N, key=.(components, chr)][, .(nchr=.N), key=.(components)][ww, on='components']->ww
ww[cl, on='components']->cl
mem[, .(rep=seq, rep.genotype=genotype)][cl, on='rep']->cl

setkey(cl, components) #information per cluster

cl[order(rep.genotype)] -> cl

#see reps per accessions
cl[, .N, by=rep.genotype]

cl[, c("i1","i2") := tstrsplit(rep, ":", fixed=TRUE)]
cl[, c("start","end") := tstrsplit(i2, "-", fixed=TRUE)]
cl[, c("gn","chr") := tstrsplit(i1, "|", fixed=TRUE)]
cl[, c("i1", "i2", "gn") := NULL]

cl[, .(chr,start,end,rep, rep.genotype)] -> sc

# save the results

# write list of cluster representatives as a master table
fwrite(cl, file="pan_cluster_mmseq_rep.txt", sep="\t", quote=F)

# write list of cluster representatives for downstream analysis
fwrite(sc, file="pan_cluster_bed_rep.txt", sep="\t", quote=F)

#exit
q()
