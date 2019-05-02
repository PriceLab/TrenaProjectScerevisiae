library(AnnotationHub)
library(org.Sc.sgd.db)


hub = AnnotationHub()
snapshotDate(): 2019-04-29
length(query(hub, c("ensembl", "gtf", "release-96")))  # [1] 471
length(query(hub, c("fasta", "release-96", "twobit"))) # [1] 799

as.list(query(hub, c("ensembl", "gtf", "cerevisiae", "R64")))
       # retrieve records with 'object[["AH64833"]]'
x <- hub[["AH64833"]]   # length: 41606
tbl.x <- as.data.frame(x)
dim(tbl.x)  # [1] 41606    18
tbl.geneInfo <- subset(tbl.x, type=="gene") # [1] 7036   18
dim(tbl.geneInfo)
orfs <- tbl.geneInfo$gene_id
length(orfs)          #  7036
length(unique(orfs))  #  7036

# keytypes(org.Sc.sgd.db)
# demonstrate
tbl.names <- select(org.Sc.sgd.db, keys=head(unique(tbl.x$gene_id)), columns="GENENAME", keytype="ORF")
  #             ORF        SGD GENENAME
  #    1   YDL248W S000002407     COS7
  #    2 YDL247W-A S000007602     <NA>

tbl.names <- select(org.Sc.sgd.db, keys=orfs, columns="GENENAME", keytype="ORF")
dim(tbl.names)
unmapped <- which(is.na(tbl.names$GENENAME))
length(unmapped)
tbl.names$GENENAME[unmapped] <- tbl.names$ORF[unmapped]

tbl.merged <- merge(tbl.geneInfo, tbl.names, by.x="gene_id", by.y="ORF")
tbl.merged.2 <- tbl.merged[, c("seqnames", "start", "end", "strand", "width", "GENENAME", "type", "gene_biotype", "gene_id")]
tbl.merged.2$tss <- tbl.merged.2$start
reverse.strand <- which(tbl.merged.2$strand == "-")
length(reverse.strand)
tbl.merged.2$tss[reverse.strand] <- tbl.merged.2$end[reverse.strand]

colnames(tbl.merged.2) <- c("chrom", "start", "end", "strand", "width", "geneSymbol", "type", "gene_biotype", "geneID", "tss")
tbl.geneInfo <- tbl.merged.2[, c("chrom", "start", "end", "strand", "tss", "geneID", "geneSymbol", "type", "gene_biotype", "width")]
tbl.geneInfo$chrom <- as.character(tbl.geneInfo$chrom)
tbl.geneInfo$strand <- as.character(tbl.geneInfo$strand)
tbl.geneInfo$type <- as.character(tbl.geneInfo$type)
save(tbl.geneInfo, file="../../../inst/extdata/geneInfoTable.RData")
