#Script for Chip-seq data analysis using Bioconductor and ChIPSeeker
#source("https://bioconductor.org/biocLite.R")
#biocLite("ChIPseeker")
#devtools::install_github("GuangchuangYu/ChIPseeker")

library(GenomicFeatures)
library(biomaRt)

setwd("C:/cygwin64/home")

#txdb<-makeTranscriptDbFromBiomart(biomart ="fungal_mart" ,dataset="spombe_eg_gene" ,host="fungi.ensembl.org")
txdb <- makeTxDbFromGFF("genes.gtf",organism="Schizosaccharomyces pombe")

library(ChIPseeker)

peakC <- readPeakFile("MACS_C2_peaks.bed",as = "GRanges")
peakB <- readPeakFile("MACS_B_peaks.bed",as = "GRanges")
#head(peakC)
#Do wtkresu porownawczego heatmapy
peak_shuffled <- readPeakFile("peaks_shuffled.bed",as = "GRanges")

#covplot function calculates the coverage of peak regions over chromosomes and generate a figure to visualize
covplot(peakC, weightCol="V5",chrs=c("I", "II", "III"))
peaks=GenomicRanges::GRangesList(B=peakB,C=peakC)

#Profile of ChIP peaks binding to TSS regions
#for calculate the profile of ChIP peaks binding to TSS regions, we should prepare the TSS regions, which are defined as the flanking sequence of the TSS sites. Then align the peaks that are mapping to these regions, and generate the tagMatrix

promoter <- getPromoters(TxDb=txdb, upstream=500, downstream=500)
tagMatrixC <- getTagMatrix(peakC, windows=promoter)
tagMatrixB <- getTagMatrix(peakB, windows=promoter)

plotAvgProf(tagMatrix, xlim=c(-500, 500),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrixB, xlim=c(-500, 500), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrixC, xlim=c(-500, 500), conf = 0.95, resample = 1000)

#tagHeatmap(tagMatrixB, xlim=c(-500, 500), color="cyan3",title="Heatmap of ChIP binding to TSS regions - sample B")
#tagHeatmap(tagMatrixC, xlim=c(-500, 500), color="cyan3",title="Heatmap of ChIP binding to TSS regions - sample C")
peakHeatmap(peakB, TxDb=txdb, upstream=1000, downstream=1000, color="cyan3",title="Heatmap of ChIP binding to TSS regions - sample B")
peakHeatmap(peakC, TxDb=txdb, upstream=1000, downstream=1000, color="cyan3",title="Heatmap of ChIP binding to TSS regions - sample C")
peakHeatmap(peak_shuffled, TxDb=txdb, upstream=1000, downstream=1000, color="cyan3",title="Heatmap of ChIP binding to TSS regions - peaks shuffled")

#Adnotacja dla S pombe https://www.bioconductor.org/packages/devel/data/annotation/html/MeSH.Spo.972h.eg.db.html
biocLite("MeSH.Spo.972h.eg.db")
library(MeSH.Spo.972h.eg.db)
MeSH.Spo.972h.eg.db

peakAnnoC <- annotatePeak(peakC, tssRegion=c(-500, 500),TxDb=txdb, annoDb="MeSH.Spo.972h.eg.db")
peakAnnoB <- annotatePeak(peakB, tssRegion=c(-500, 500),TxDb=txdb, annoDb="MeSH.Spo.972h.eg.db")

plotAnnoPie(peakAnnoB)
plotAnnoPie(peakAnnoC)

plotAnnoBar(peakAnnoB)
plotAnnoBar(peakAnnoC)

upsetplot(peakAnnoB, vennpie=FALSE)
upsetplot(peakAnnoC, vennpie=FALSE)

#średnia odleglosc peakow od najblizszego genu
as.data.frame(peakAnnoB)["distanceToTSS"]
summary(as.data.frame(peakAnnoB)["distanceToTSS"])
summary(as.data.frame(peakAnnoC)["distanceToTSS"])

plotDistToTSS(peakAnnoB,title="Distribution of transcription factor-binding loci\nrelative to TSS - sample B")
plotDistToTSS(peakAnnoC,title="Distribution of transcription factor-binding loci\nrelative to TSS - sample C")

as.vector(as.data.frame(peakAnnoC)["distanceToTSS"])
hist(c(as.data.frame(peakAnnoB)["distanceToTSS"])$distanceToTSS)
b_dist <- as.data.frame(peakAnnoB)["distanceToTSS"]
c_dist <- as.data.frame(peakAnnoC)["distanceToTSS"]
b_dist$sample <- 'B'
c_dist$sample <- 'C'
distanceToTSS <- rbind(b_dist,c_dist)
ggplot(distanceToTSS, aes(distanceToTSS, fill = sample)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
ggplot(data=distanceToTSS, aes(distanceToTSS)) + geom_histogram()

qplot(distanceToTSS$distanceToTSS,
      geom="histogram",  
      main = "Histogram of distance to nearest gene and genomic region", 
      xlab = "Distance to nearest gene",  
      fill=I("blue"), 
      col=I("cyan3"), 
      alpha=I(.2))
#B
qplot(b_dist$distanceToTSS,
      geom="histogram",  
      main = "Histogram of distance to nearest gene and genomic region", 
      xlab = "Distance to nearest gene",  
      fill=I("cadetblue1"), 
      col=I("cyan"), 
      alpha=I(.4))
#C
qplot(c_dist$distanceToTSS,
      geom="histogram",  
      main = "Histogram of distance to nearest gene and genomic region", 
      xlab = "Distance to nearest gene",  
      fill=I("cadetblue1"), 
      col=I("cyan"), 
      alpha=I(.4))

#enrichment analysis
library("ReactomePA")
as.data.frame(peakAnnoB)$geneId
write(as.data.frame(peakAnnoC)$geneId,file="genesC.txt", sep="\n")


tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
tagHeatmap(tagMatrixList, xlim=c(-500, 500), color=NULL)
plotAvgProf(tagMatrixList, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb, tssRegion=c(-500, 500), verbose=FALSE)
plotDistToTSS(peakAnnoList)

library("clusterProfiler")
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))


compKEGG <- compareCluster(geneCluster   = genes, organism = "spo",
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
plot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
#AnnotationHub
biocLite("AnnotationHub")
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, "pombe")
# snapshotDate(): 2016-10-11 
# names(): AH10598
# $dataprovider: Inparanoid8
# $species: Schizosaccharomyces pombe
# $rdataclass: Inparanoid8Db
# $title: hom.Schizosaccharomyces_pombe.inp8.sqlite
# $description: Inparanoid 8 annotations about Schizosaccharomyces pombe
# $taxonomyid: 284812
# $genome: inparanoid8 genomes
# $sourcetype: Inparanoid
# $sourceurl: http://inparanoid.sbc.su.se/download/current/Orthologs/S.pombe
# $sourcelastmodifieddate: NA
# $sourcesize: NA
# $tags: c("Inparanoid", "Gene", "Homology", "Annotation") 
# retrieve record with 'object[["AH10598"]]
pombe <- hub[["AH10598"]]

keytypes(pombe)
keytypes(MeSH.Spo.972h.eg.db)
#nie ma takich id w obiekcie OrgDb
#gene.df <- bitr(genes, fromType = "SCHIZOSACCHAROMYCES_POMBE", toType = "ENSEMBL", "SYMBOL", OrgDb = pombe)

ego <- enrichGO(gene = genes, universe = names(genes), OrgDb = pombe, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

compGO <- compareCluster(geneCluster   = genes, OrgDb = "MeSH.Spo.972h.eg.db",
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")

#Sprawdzenie wersji pakietow
sessionInfo()
methods(class=class(peakAnnoB))
.libPaths()
listMarts(mart = NULL, host="fungi.ensembl.org", path="/biomart/martservice")

#Development mode
#install.packages('devtools')
#library(devtools)
#dev_mode(TRUE)
#devtools::install_github("GuangchuangYu/ChIPseeker")
#library(ChIPseeker)
#dev_mode(on=F)

#enrichment analisys
go_df = read.xls("go.xlsx", sheet = 1, header = TRUE)
data_go <- as.data.frame(genes$B)
colnames(data_go) <- c("idB")
key_go <- go_df
#data$term <- key[match(data_go$idB, key_go$GO_name), 'systematic_id']
merge(data_go, key_go, by.x = "idB", by.y = "systematic_id")
#dla C
data_goC <- as.data.frame(genes$C)
colnames(data_goC) <- c("idC")
goC <- merge(data_goC, key_go, by.x = "idC", by.y = "systematic_id")

#Histogram szerokości peaków
hist(width(peakB))
hist(width(peakC))

b_width <- as.data.frame(width(peakB))
colnames(b_width) <- c("width")
c_width <- as.data.frame(width(peakC))
colnames(c_width) <- c("width")
b_width$sample <- 'B'
c_width$sample <- 'C'
widthAll <- rbind(b_width,c_width)
ggplot(width, aes(width, fill = sample)) + geom_histogram(alpha = 0.5, position = 'identity')

qplot(c_width$width,
      geom="histogram",  
      main = "Histogram of peak widths",  
      xlab = "Width",
      fill=I("cadetblue1"), 
      col=I("cyan"), 
      alpha=I(.4))


select(txdb,keys=genes$B,columns=c("TXSTART","TXEND","TXCHROM"),keytype="GENEID") -> genesB_coord
select(txdb,keys=genes$C,columns=c("TXSTART","TXEND","TXCHROM"),keytype="GENEID") -> genesC_coord

write.csv(genesC_coord,file="genesC_coord.csv")
write.csv(genesB_coord,file="genesB_coord.csv")


#Różnice w pokryciu genów przypisanych do peaków
pokB <- read.table(file="pokrycieB.txt")
colnames(pokB) <- c("gene","coverage")
pokB_input <- read.table(file="pokrycieB_input.txt")
colnames(pokB_input) <- c("gene","coverage")
pokB$coverage_input <- pokB_input$coverage
pokB$diff <- pokB$coverage - pokB$coverage_input
write.csv(pokB,file="cov_diffB.csv")

#wykres fpkm dla wszystkich genów
library("ggplot2", lib.loc="~/R/win-library/3.3")
cuff_table <- read.table(file="genes_fpkm_tracking.txt",header=TRUE)
head(cuff_table$FPKM)
df_all_fpkm <- data.frame(fpkm = cuff_table$FPKM, Sample = 'all')
#dla ilu z B są wyniki z cufflings
sum(cuff_table$tracking_id %in% genes$B)
#27
sum(cuff_table$tracking_id %in% genes$C)
#187
b_fpkm <- cuff_table$FPKM[cuff_table$tracking_id %in% genes$B]
c_fpkm <- cuff_table$FPKM[cuff_table$tracking_id %in% genes$C]
df_b_fpkm <- data.frame(fpkm = b_fpkm, Sample = 'B')
df_c_fpkm <- data.frame(fpkm = c_fpkm, Sample = 'C')
df_fpkm <- rbind(df_all_fpkm,df_b_fpkm,df_c_fpkm)
qplot(fpkm, colour=factor(Sample), data=df_fpkm, geom="density") + xlim(0,500)
df_allc <- rbind(df_all_fpkm,df_c_fpkm)
qplot(fpkm, colour=factor(Sample), data=df_allc, geom="density") + xlim(0,500)
#box plot wartosci fpkm
ggplot(data=df_fpkm, aes(x=Sample, y=fpkm, fill=Sample)) + geom_boxplot()
df_bc <- rbind(df_b_fpkm,df_c_fpkm)
ggplot(data=df_bc, aes(x=Sample, y=fpkm, fill=Sample)) + geom_boxplot() + ylim(0,1000)
qplot(fpkm, colour=factor(Sample), data=df_bc, geom="density") + xlim(0,500)
############
head(cuff_table$tracking_id)
startsWith(as.character(cuff_table$tracking_id[1]), 'CUFF')
#novel transcripts
sum(startsWith(as.character(cuff_table$tracking_id), 'CUFF'))
#3037
length(cuff_table$tracking_id)
#5324