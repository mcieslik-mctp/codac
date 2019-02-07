readCounts <- function(dir, ann) {
    fn <- grep("_(gene|full)-", dir$cts.fn, value=TRUE)
    tbl <- data.table::fread(cmd=sprintf("zcat %s", fn), skip=1, drop=c(2,3,4,5), col.names = c("gene_id", "length", "count"))
    ## add locus_id
    setkey(tbl, gene_id)
    setkey(ann$gene.ovr, gene_id)
    tbl <- ann$gene.ovr[tbl]
    tbl <- tbl[, .(count=sum(count)),by=locus_id]
    setkey(tbl, locus_id)
    tbl <- tbl[ann$loci$locus_id]
    tbl[,cpm:=count/(sum(count, na.rm=TRUE)/1e6),]
    setkey(tbl, locus_id)
    return(tbl)
}
