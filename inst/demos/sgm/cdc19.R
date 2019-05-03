library(MotifDb)
library(trena)
library(trenaSGM)
library(org.Sc.sgd.db)
library(TrenaProjectScerevisiae)

genome <- "SacCer3"
targetGene <- "CDC19"
tp <- TrenaProjectScerevisiae()
setTargetGene(tp, targetGene)
tbl.transcript <- getTranscriptsTable(tp)

tbl.regions <- getGeneRegulatoryRegions(tp)
tbl.regions$chrom <- sprintf("chr%s", tbl.regions$chrom)

pfms <- as.list(jaspar.yeast.pfms <- query(MotifDb, c("cerevisiae", "jaspar2018")))

motifMatcher <- MotifMatcher(genomeName="SacCer3", pfms)
tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=95L)
dim(tbl.motifs)
checkTrue(nrow(tbl.motifs) > 80 & nrow(tbl.motifs) < 120)
checkEquals(ncol(tbl.motifs), 13)
tf.geneSymbols <- mcols(MotifDb[tbl.motifs$motifName])$geneSymbol
tbl.motifs$targetGene <- targetGene
tbl.motifs$tf <- tf.geneSymbols

candidate.tfs <- sort(mcols(MotifDb[unique(tbl.motifs$motifName)])$geneSymbol)
length(candidate.tfs)
getExpressionMatrixNames(tp)
mtx <- getExpressionMatrix(tp, "GSE115556.asinh")
dim(mtx)

build.spec <- list(title="CDC19.noDNA.allTFs",
                   type="noDNA.tfsSupplied",
                   matrix=mtx,
                   tfPool=candidate.tfs,
                   tfPrefilterCorrelation=0.0,
                   annotationDbFile=dbfile(org.Sc.sgd.db),
                   orderModelByColumn="rfScore",
                   solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                   quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
x <- build(builder)
fivenum(x$model$rfScore)
tbl.model <- subset(x$mode, rfScore > fivenum(x$model$rfScore)[4])
dim(tbl.model)
tbl.regRegions <- subset(tbl.motifs, tf %in% tbl.model$gene)
dim(tbl.regRegions)

library(igvR)
igv <- igvR()
getSupportedGenomes(igv)
setGenome(igv, "sacCer3")
loc <- with(tbl.regions, sprintf("%s:%d-%d", chrom, start, end))
showGenomicRegion(igv, loc)

library (RColorBrewer)
totalColorCount <- 10
colors <- brewer.pal(totalColorCount, "Paired")
currentColorNumber <- 0

tfs <- sort(unique(tbl.regRegions$tf))
for(tf.geneSymbol in tfs){
   tbl.sub <- subset(tbl.regRegions, tf==tf.geneSymbol)[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
   currentColorNumber <- (currentColorNumber %% totalColorCount) + 1
   color <- colors[currentColorNumber]
   track <- DataFrameQuantitativeTrack(tf.geneSymbol, tbl.sub, autoscale=TRUE, color=color)
   displayTrack(igv, track)
   } # for tf
