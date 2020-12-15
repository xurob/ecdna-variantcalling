library(GenomicRanges)
library(Rsamtools)
library(bamsignals)


MYCN <- GRanges(seqnames=Rle(c("2"), c(1)),
                ranges=IRanges(16080686:16087129),
                strand=rep(c("+"), 1))



MYCN

chrM <- GRanges(seqnames=Rle(c("MT"), c(1)),
                ranges=IRanges(30:16500))
chrM

files <- list.files(path = "/fast/users/xuro_c/scratch/vc_mito_mycn/mycn/DIR_rmdup_index", pattern=".bam", all.files=FALSE,
                    full.names=TRUE)
files
lapply(files, function(x) {
  y <- toString(x)
  t <- gsub(".bam","",y,ignore.case=".bam")
  t <- gsub("/fast/users/xuro_c/scratch/vc_mito_mycn/mycn/DIR_rmdup_index/","",y,ignore.case=".bam")
  t <- paste0(t, ".png")
  
  png(filename = file.path("/fast/users/xuro_c/scratch/vc_mito_mycn/mycn/DIR_rmdup_index/plots/",t), units="in", width=5, height=5, res=300)
  covSigs <- bamCoverage(x, MYCN, verbose=FALSE)
  xlab <- "offset from start of the region"
  ylab <- "reads per base pair"
  main <- paste0("read coverage MYCN  Chr", seqnames(MYCN),
                 ":", start(MYCN), "-", end(MYCN))
  plot(covSigs[1], ylim=c(0, max(covSigs[1])), ylab=ylab, xlab=xlab, main=main,
       type="l")
  dev.off()
  
})

