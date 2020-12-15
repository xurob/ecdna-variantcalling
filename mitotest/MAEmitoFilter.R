library(SummarizedExperiment)
library(BuenColors)
library(MultiAssayExperiment)
library(gplots)
library(BioQC)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(plotly)
library(RColorBrewer)


"%ni%" <- Negate("%in%")

# Import Caleb's genotype calls
raw <- readRDS("/Users/robin/Documents/projects/supervisedgenotypingpipeline/testdir/MAEopen.rds")

# Import Lauren's clonal definitions; make a vector for quick lookup
meta_df <- data.frame(fread("/Users/robin/Documents/projects/supervisedgenotypingpipeline/testdir/CellID.txt"))
clonevec <- as.character(meta_df$Clonotype.ID); names(clonevec) <- as.character(meta_df$Cell.Name)

# Collapse multi assay experiments into SummarizedExperiments
covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
#baq <- assays(allSE)[["BAQ"]] # don't need this for now.

# Pull out the allele frequencies and the coverages per position / allele
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], data.frame(rowRanges(allSE))[,c(6)], ">", data.frame(rowRanges(allSE))[,c(7)])

# Compute per-cell and per-variant properties
covVARmean <- Matrix::rowMeans(assays(covSE)[["coverage"]][start(rowRanges(allSE)),])
covVARmedian <- matrixStats::rowMedians(data.matrix(assays(covSE)[["coverage"]][start(rowRanges(allSE)),]))
covCell <-  Matrix::colMeans(assays(covSE)[["coverage"]])

# Filter variants and cells that don't have great coverages
minCellCoverage <- 1
af2 <- af[covVARmean> 1  & covVARmedian > 1 , covCell > minCellCoverage]
cov2 <- cov[covVARmean> 1 & covVARmedian > 1, covCell > minCellCoverage]
dim(af2)

# Filter for the minor allele
af3 <- af2[Matrix::rowMeans(af2) > 0.0001 & Matrix::rowMeans(af2) < 0.75, ]
cov3 <- cov2[Matrix::rowMeans(af2) > 0.0001 & Matrix::rowMeans(af2) < 0.75, ]

dim(af3)

# Adjust column names
colnames(af3) <- gsub(".mito", "", colnames(af3))
colnames(cov3) <- gsub(".mito", "", colnames(cov3))

# Remove samples where we don't have a clone ID... shouldn't be any
af3 <- af3[,colnames(af3) %in% names(clonevec)]
cov3 <- cov3[,colnames(cov3) %in% names(clonevec)]
dim(af3)

# Create a vector per cell with clone ID
clones <- clonevec[colnames(af3)]

# Identify all clones with at least 5 cells
data.frame(clones) %>% group_by(clones) %>%
  summarise(count = n()) %>% filter(count >=1) %>%
  pull(clones) %>% as.character() -> search_clones

# Loop over these clones to identify variants
lapply(search_clones, function(clone){
  
  # Do a Mann-Whitney test
  out <- sort(wmwTest(data.matrix(t(af3)), which(clones == clone),
                      valType = "abs.log10p.two.sided"), decreasing = TRUE)
  
  # Retain if reasonable p-value-- will threshold harder, later
  keep <- out[out> 3]
  df <- data.frame(Variant = names(keep), log10p = unname(keep), colony = clone)
  df
})  %>% rbindlist() %>% as.data.frame() -> variantScoresClones



# Try different threshold to see what is optimal
variantScoresClones %>% filter(log10p > 3) %>% dim()
keepVariants <- unique(variantScoresClones[variantScoresClones$log10p > 3, "Variant"])
length(keepVariants)
keepVariants <- keepVariants[-c(11,24)]
# Only keep variants that match this
boo <- rownames(af3) %in% keepVariants
af3 <- data.matrix(af3[boo,])
dim(af3)

# Threshold allele frequencies for a plot -- higher heteroplasmy values are definitely interesting
# but we will look at those case-by-case in box plots like in the existing manuscript
af4 <- af3
#af4[af4 > .2] <- .2

# Make an ordered meta-data file
#meta_df %>% filter(Cell.Name %in% colnames(af3)) %>% 
  #mutate(clone_color = Component.Size >= 0, Clonotype.ID, "small_clone") %>%
 # arrange(clone_color) -> annodf2


meta_df %>% filter(Cell.Name %in% colnames(af3)) %>% 
mutate(clone_color = c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3), Clonotype.ID, "small_clone") %>%
arrange(clone_color) -> annodf2

unique_clone_colors <- unique(as.character(annodf2$clone_color))

# Make a random vector of colors for annotation
colorsUniverse <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols <- sample(colorsUniverse, length(unique_clone_colors))
names(cols) <- unique_clone_colors

# That's it! Now make a heatmap
ha1 <- HeatmapAnnotation(df = annodf2[,c("clone_color"), drop = FALSE],
                         col = list(clone_color = cols)
)

pdf(file="../output/plots/heatmaps/P0322-heatmaps/P0322.heatmap.pdf", height = 9, width = 15)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap(sqrt(af4[,as.character(annodf2$Cell.Name)]),
        col=as.character(jdb_palette("brewer_red",type="continuous")),
        cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha1,
        name = "AF")
dev.off()
# New heatmap but with actual allele frequencies
pdf(file="../output/plots/heatmaps/P0322.small.heatmap.pdf", height = 5, width = 15)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap((af3[,as.character(annodf2$Cell.Name)])^(1/2),
        col=as.character(jdb_palette("brewer_red",type="continuous")),
        cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha1,
        name = "sqrt(AF)")
dev.off()


bild1 <- Heatmap((af3[,as.character(annodf2$Cell.Name)])^(1/3),clustering_distance_columns = function(x, y) 1 - cor(x, y),
                 col=as.character(jdb_palette("brewer_red",type="continuous")),
                 cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = TRUE,
                 row_names_gp = gpar(fontsize = 6),
                 top_annotation = ha1,
                 name = "AF^(1/3)")

bild1

 
 
af5 <- af3[,as.character(annodf2$Cell.Name)]
af6 <- dist(t(sqrt(af5)),method="euclidean") 
tr <- hclust(af6, method = "complete")
tr
plot(tr)
