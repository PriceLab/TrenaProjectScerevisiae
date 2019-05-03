library(TrenaProjectScerevisiae)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tProj")) {
   message(sprintf("--- creating instance of TrenaProjectScerevisiae"))
   tProj <- TrenaProjectScerevisiae();
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_supportedGenes()
   test_variants()
   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()
   test_getGeneRegulatoryRegions()

   test_buildSingleGeneModel()


} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectScerevisiae", "TrenaProject") %in% is(tProj)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("CDC19")
   checkTrue(all(subset.expected %in% getSupportedGenes(tProj)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(getVariantDatasetNames(tProj), character(0))

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   checkTrue(is.na(getFootprintDatabaseNames(tProj)))
   checkTrue(is.na(getFootprintDatabaseHost(tProj)))

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("GSE115556.asinh")
   checkTrue(all(expected %in% getExpressionMatrixNames(tProj)))

   mtx <- getExpressionMatrix(tProj, expected[1])
   checkEquals(dim(mtx), c(6692, 295))

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tProj, "CDC19")
   checkEquals(getTargetGene(tProj), "CDC19")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tProj)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(tbl.transcripts$chr, "I")
   checkEquals(tbl.transcripts$start, 71786)
   checkEquals(tbl.transcripts$end , 73288)
   checkEquals(tbl.transcripts$tss, 71786)
   checkEquals(tbl.transcripts$strand, "+")

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tProj, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "I:71786-73288")


   #message(sprintf("    geneGeneEnhancersRegion"))
   #region <- getGeneEnhancersRegion(tProj, flankingPercent=0)
   #checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   #checkEquals(region$chromLocString, "chr4:52901678-53289004")

   #message(sprintf("    encode DHS"))
   #tbl.dhs <- getEncodeDHS(tProj)
   #checkEquals(nrow(tbl.dhs), 0)

   #message(sprintf("    ChIP-seq"))
   #tbl.chipSeq <- getChipSeq(tProj, chrom=chromosome, start=start, end=end, tfs=NA)
   #checkEquals(nrow(tbl.chipSeq), 0)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
test_getGeneRegulatoryRegions <- function()
{
   printf("--- test_getGeneRegulatoryRegions")

   goi <- "CDC19"  # plus strand
   tbl <- getGeneRegulatoryRegions(tProj, goi)   # default is -2000, +500

   checkEquals(tbl$chrom, "I")
   checkEquals(tbl$start, 69786)
   checkEquals(tbl$end, 72286)
   checkTrue(tbl$start < tbl$end)   # since + strand
   checkEquals(with(tbl, end-start), 2500)
   checkEquals(tbl$geneSymbol, goi)
   checkEquals(tbl$type, "Promoter")

   goi <- "HPA2"  # on the minus strand
   tbl <- getGeneRegulatoryRegions(tProj, goi)
   checkEquals(tbl$chrom, "XVI")
   checkEquals(tbl$start, 925379)
   checkEquals(tbl$end,   922879)
   checkTrue(tbl$start > tbl$end)   # since + strand
   checkEquals(with(tbl, start-end), 2500)
   checkEquals(tbl$geneSymbol, goi)
   checkEquals(tbl$type, "Promoter")

   tbl <- getClassicalGenePromoter(tProj, goi, 0, 0)
   checkEquals(tbl$start, tbl$end)
   tss <- getTranscriptsTable(tProj, goi)$tss
   checkEquals(tbl$start, tss)

   tbl.goi <- getClassicalGenePromoter(tProj, goi, 100, 10)
   with(tbl.goi, {
           checkEquals(start, 923479);
           checkEquals(end,   923369);
           checkEquals(chrom, "XVI")
           })

   tbl.targetGene <- getClassicalGenePromoter(tProj, targetGene=NA, 10, 10)
   with(tbl.targetGene, {
           checkEquals(start, 71776);
           checkEquals(end,   71796);
           checkEquals(chrom, "I")
           })

} # test_getGeneRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel <- function()
{
   printf("--- test_buildSingleGeneModel")

   genome <- "SacCer3"
   targetGene <- "CDC19"
   tp <- TrenaProjectScerevisiae()
   setTargetGene(tp, targetGene)
   tbl.transcript <- getTranscriptsTable(tp)

   tbl.regions <- getGeneRegulatoryRegions(tp)
   tbl.regions$chrom <- sprintf("chr%s", tbl.regions$chrom)

   library(MotifDb)
   library(trena)
   library(trenaSGM)
   library(org.Sc.sgd.db)
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
   getExpressionMatrixNames(tProj)
   mtx <- getExpressionMatrix(tProj, "GSE115556.asinh")
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

} # test_buildSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
