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
   checkEquals(region$chromLocString, "chr4:53030787-53159895")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tProj)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 0)

     # no enhancers, DHS or ChIP-seq for mouse

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
if(!interactive())
   runTests()
