file <- "GSE115556_normalized_counts_matrix.csv"
tbl <- read.table(file, sep=",", as.is=TRUE, nrow=-1, header=TRUE)
rownames(tbl) <- tbl[, 1]
mtx <- as.matrix(tbl[,-1])
dim(mtx)           # 6692 295
mtx[1:10, 1:10]
goi <- "CDC19"
goi %in% rownames(mtx)

plot(mtx[goi,])
fivenum(mtx)          # 0.00000     72.39078    283.97727    700.83245   489829.29230
mtx <- asinh(mtx)     # 0.000000     4.975274     6.342044     7.245417      13.794959

save(mtx, file="../../../inst/extdata/expression/GSE115556.asinh.RData")

