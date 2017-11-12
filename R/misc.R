.naTo0 <- function(DT, omit=NULL, zero=0) {
    for (i in setdiff(names(DT), omit))
        DT[is.na(get(i)),(i):=zero]
    return(DT)
}

.naToBlank <- function(DT, omit=NULL) {
    for (i in setdiff(names(DT), omit))
        DT[is.na(get(i)),(i):=""]
    return(DT)
}

.omitBlank <- function(x) {
    return(x[!(x=="")])
}

.dustySplit <- function(x, batch=100000) {
     unlist(unname(lapply(split(x, factor(ceiling(seq_along(x) / batch))), dustyScore)))
}

.l2g.full <- function(ann) {
    tmp <- data.table(as.data.frame(mcols(ann$genes)))
    setkey(tmp, gene_id)
    setkey(ann$gene.ovr, gene_id)
    l2g <- ann$gene.ovr[tmp]
    setkey(ann$gene.ovr, locus_id)
    return(l2g)
}

.l2g <- function(ann, only.prot=FALSE) {
    sel <- if(only.prot) quote(protein) else quote(TRUE)
    l2g <- .l2g.full(ann)
    l2g <- l2g[,list(
      gene_ids=paste0(gene_id[eval(sel)], collapse=":"),
      gene_names=paste0(gene_name[eval(sel)], collapse=":")), by=locus_id]
    setkey(l2g, "locus_id")
    return(l2g)
}

.lgt <- function(ann) {
    go <- ann$gene.ovr
    setkey(go, gene_id)
    ge <- as.data.table(mcols(ann$genes)[,c("gene_id", "gene_name")])
    setkey(ge, gene_id)
    go <- ge[go]
    gt <- data.table(as.data.frame(mcols(ann$transcripts)))
    setkey(go, gene_id)
    setkey(gt, gene_id)
    lgt <- go[gt]
    lgt[,tags:=factor(
        ifelse(grepl("appris_principal_1", tags), "ap1",
        ifelse(grepl("appris_principal_2", tags), "ap2",
        ifelse(grepl("appris_principal_3", tags), "ap3",
        ifelse(grepl("CCDS", tags), "rev",
        ifelse(cdsN!=0, "cds", "nil"))))),
        levels=c("ap1", "ap2", "ap3", "rev", "cds", "nil"))]
    setkey(lgt, locus_id)
    return(lgt)
}
