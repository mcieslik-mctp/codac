.readSJ <- function(dir, ann) {
    tmp <- fread(paste0(system2("zcat", dir$sj.fn, stdout=TRUE), collapse="\n", showProgress=FALSE),
                 showProgress=FALSE)
    tmp <- tmp[(V7 > 1)] # unique reads supporting junction w/ motif
    sj <- GRanges(seqnames=tmp$V1, IRanges(start=tmp$V2, end=tmp$V3),
                  strand=ifelse(tmp$V4==1,"+",ifelse(tmp$V4==2, "-", "*")),
                  nfrag=tmp$V7, motif=(tmp$V5 %in% c(1, 2)))
    seqlevels(sj, pruning.mode="coarse") <- seqlevels(ann$seqi)
    seqinfo(sj) <- ann$seqi
    return(sj)
}

#' @export
readSplices <- function(dir, ann) {
    ## locsjexp
    sjs <- .readSJ(dir, ann)
    hits <- suppressWarnings(data.table(data.frame(findOverlaps(ann$loci, sjs))))
    hits$nfrag <- mcols(sjs[hits$subjectHits])$nfrag
    tmp <- hits[,list(nsj=.N, nfrag=sum(nfrag)),by=queryHits]
    tmp$locus_id <- mcols(ann$loci[tmp$queryHits])$locus_id
    sjexp <- data.table(
      sid=dir$sid,
      locus_id=ann$loci$locus_id,
                     nsfrag = 0L,
                        nsj = 0L)
    setkey(tmp, "locus_id")
    setkey(sjexp, "locus_id")
    sjexp[tmp,":="(nsfrag=i.nfrag,nsj=i.nsj)]
    sjexp[,cspm:=nsfrag / (sum(nsfrag) / 1e6)]
    setkeyv(sjexp, c("sid", "locus_id"))
    spl <- list(
      sjs=sjs,
      locus.sjexp=sjexp
    )
    return(spl)
}
