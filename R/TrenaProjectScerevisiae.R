#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectScerevisiae-class
#'
#' @name TrenaProjectScerevisiae-class
#' @rdname TrenaProjectScerevisiae-class
#' @aliases TrenaProjectScerevisiae
#' @exportClass TrenaProjectScerevisiae
#'

.TrenaProjectScerevisiae <- setClass("TrenaProjectScerevisiae",
                                  contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectScerevisiae
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectScerevisiae-class
#'
#' @export
#'
#' @return An object of the TrenaProjectScerevisiae class
#'

TrenaProjectScerevisiae <- function(quiet=TRUE)

{
   genomeName <- "mm10"

   directory <- system.file(package="TrenaProjectScerevisiae", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()
   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE)
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   geneInfoTable.path <- system.file(package="TrenaProjectScerevisiae", "extdata", "geneInfoTable.RData")

   footprintDatabaseNames <- NA_character_;
   expressionDirectory <- system.file(package="TrenaProjectScerevisiae", "extdata", "expression")
   stopifnot(file.exists(expressionDirectory))

   variantsDirectory <- system.file(package="TrenaProjectScerevisiae", "extdata", "variants")

   covariatesFile <- NA_character_;

   .TrenaProjectScerevisiae(TrenaProject(projectName="TrenaProjectScerevisiae",
                                         supportedGenes=geneSets[[1]],
                                         genomeName="SacCer3",
                                         geneInfoTable.path=geneInfoTable.path,
                                         footprintDatabaseHost=NA_character_,
                                         footprintDatabaseNames=NA_character_,
                                         footprintDatabasePort=NA_integer_,
                                         expressionDirectory=expressionDirectory,
                                         variantsDirectory=variantsDirectory,
                                         covariatesFile=covariatesFile,
                                         quiet=quiet
                                         ))

} # TrenaProjectScerevisiae, the constructor
#----------------------------------------------------------------------------------------------------
