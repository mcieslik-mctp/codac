#' @export
bsFormat <- function(tbl) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        breakpoints=sum.bpt,
        distance=tx.dst,
        spliced=d2a,
        "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
        "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
        "repetitive(5';3)'" = paste(art.5, art.3, sep=";"),
        chain=bs.chain
    )][order(chain)]
    return(rep)
}

#' @export
slFormat <- function(tbl) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        breakpoints=sum.bpt,
        spliced=ifelse(d2a, "✓", ""),
        HQ=ifelse(hq.bpt, "✓", ""),
        inframe=ifelse(orf, "✓", ""),
        distance=dst,
        topology=topo,
        "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
        chain=sl.chain
    )][order(chain)]
    return(rep)
}

#' @export
svFormat <- function(tbl) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        encompassing.reads=sum.enc,
        breakpoints=sum.bpt,
        spliced=ifelse(d2a, "✓", ""),
        HQ=ifelse(hq.bpt, "✓", ""),
        HR=ifelse(hr.bpt, "✓", ""),
        inframe=ifelse(orf, "✓", ""),
        distance=dst,
        topology=topo,
        "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
        "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
        "recurrent(5';3')" = paste(unq.rec.5, unq.rec.3, sep=";"),
        "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
        chain=sv.chain,
        read.5=seq.5,
        read.3=seq.3
    )][order(chain)]
    return(rep)
}

#' @export
tsFormat <- function(tbl) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        breakpoints=sum.bpt,
        spliced=ifelse(d2a, "✓", ""),
        HQ=ifelse(hq.bpt, "✓", ""),
        HR=ifelse(hr.bpt, "✓", ""),
        inframe=ifelse(orf, "✓", ""),
        distance=dst,
        topology=topo,
        "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
        "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
        "recurrent(5';3')" = paste(unq.rec.5, unq.rec.3, sep=";"),
        "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
        chain=ts.chain
    )][order(chain)]
    return(rep)
}

#' @export
qcFormat <- function(stat, ann) {
    qc <- stat[sapply(stat, is.numeric)]
    qc.tbl <- data.table(var=names(qc), val=round(unlist(qc), 2))
    setkey(qc.tbl, var)
    setkey(ann$cutoffs, var)
    qc.tbl <- ann$cutoffs[qc.tbl]
    qc.tbl[,status:=
                ifelse(val >= error.lo & val <= error.hi, "fail",
                ifelse(val >= warn.lo  & val <=  warn.hi, "warn", "pass"))
    ]
    qc.tbl[is.na(status), status:="pass"]
    qc.tbl <- qc.tbl[,.(var, val, status)]
    return(qc.tbl)
}

#' @export
readFormat <- function(fn) {
    tbl <- fread(fn, colClasses=c(
                  "gene.5"="character", "gene.3"="character",
                  "spliced"="character", "HQ"="character",
                  "HR"="character", "inframe"="character"))
    tmp <- str_match(basename(fn), "(.*)-parc-(sv|ts|bs|sl)-rep.csv")
    tbl$sid <- tmp[,2]
    tbl$type <- tmp[,3]
    return(tbl)
}
