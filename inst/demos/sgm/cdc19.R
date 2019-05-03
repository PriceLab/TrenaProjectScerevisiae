library(MotifDb)
library(trena)
library(trenaSGM)
library(org.Sc.sgd.db)
library(TrenaProjectScerevisiae)
library(RUnit)

#------------------------------------------------------------------------------------------------------------------------
trimModel <- function(tbl.model, tbl.reg, tf.keepers=c(), votesNeeded=3)
{
   matched.keeper.rows <- unlist(lapply(tf.keepers, function(tf) grep(tf, tbl.model$gene)))
   pvals <- -log10(tbl.model$lassoPValue)
   good.lassoPval <- which(pvals > 5)
   good.betaLasso <- which(abs(tbl.model$betaLasso) > 0.1)
   good.betaRidge <- which(abs(tbl.model$betaRidge) > 0.1)
   spearman.cutoff <- fivenum(abs(tbl.model$spearmanCoeff))[4]
   good.spearmanCoeff <- which(abs(tbl.model$spearmanCoeff) >= spearman.cutoff)
   randomForest.cutoff <- fivenum(tbl.model$rfScore)[4]
   forest.tfs <- subset(tbl.model, rfScore >= randomForest.cutoff)$gene
   good.rfScore <- unlist(lapply(forest.tfs, function(tf) grep(tf, tbl.model$gene)))
   all.counts <- c(good.lassoPval, good.betaLasso, good.betaRidge, good.spearmanCoeff, good.rfScore)
   tbl.freq <- as.data.frame(table(all.counts), stringsAsFactors=FALSE)
   colnames(tbl.freq) <- c("rowNumber", "count")
   tbl.freq <- tbl.freq[order(tbl.freq$count, decreasing=TRUE),]
   tbl.freq$rowNumber <- as.integer(tbl.freq$rowNumber)
   good.tf.rows <- subset(tbl.freq, count >= votesNeeded)$rowNumber
   not.yet.included <- setdiff(matched.keeper.rows, good.tf.rows)
   if(length(not.yet.included) > 0)
      good.tf.rows <- c(good.tf.rows, not.yet.included)
   tbl.model <- tbl.model[good.tf.rows,]
   new.order <- order(abs(tbl.model$spearmanCoeff), decreasing=TRUE)
   tbl.model <- tbl.model[new.order,]
   tbl.reg <- subset(tbl.reg, tf %in% tbl.model$gene)
   return(list(model=tbl.model, regulatoryRegions=tbl.reg))

} # trimModel
#------------------------------------------------------------------------------------------------------------------------

genome <- "SacCer3"
targetGene <- "CDC19"
tp <- TrenaProjectScerevisiae()
setTargetGene(tp, targetGene)
tbl.transcript <- getTranscriptsTable(tp)

tbl.regions <- getClassicalGenePromoter(tp, upstream=2000, downstream=200)
tbl.regions$chrom <- sprintf("chr%s", tbl.regions$chrom)

pfms <- as.list(jaspar.yeast.pfms <- query(MotifDb, c("cerevisiae", "jaspar2018")))

motifMatcher <- MotifMatcher(genomeName="SacCer3", pfms)
tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=85L)
dim(tbl.motifs)
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
#fivenum(x$model$rfScore)
#tbl.model <- subset(x$mode, rfScore > fivenum(x$model$rfScore)[4])
x.trimmed <- trimModel(x$model, tbl.motifs)
tbl.model <- x.trimmed$model
dim(tbl.model)
tbl.regRegions <- x.trimmed$regulatoryRegions
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
