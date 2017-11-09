.addNames <- function(tbl, ann, only.prot) {
    ## locus_id to gene_name
    l2g <- .l2g(ann, only.prot)
    setkey(tbl, "locus_id.3.2")
    tbl <- l2g[tbl]
    setnames(tbl, c("locus_id", "gene_names", "gene_ids"), c("locus_id.3.2", "gene_names.3.2", "gene_ids.3.2"))
    setkey(tbl, "locus_id.5.2")
    tbl <- l2g[tbl]
    setnames(tbl, c("locus_id", "gene_names", "gene_ids"), c("locus_id.5.2", "gene_names.5.2", "gene_ids.5.2"))
    setkey(tbl, "locus_id.3.1")
    tbl <- l2g[tbl]
    setnames(tbl, c("locus_id", "gene_names", "gene_ids"), c("locus_id.3.1", "gene_names.3.1", "gene_ids.3.1"))
    setkey(tbl, "locus_id.5.1")
    tbl <- l2g[tbl]
    setnames(tbl, c("locus_id", "gene_names", "gene_ids"), c("locus_id.5.1", "gene_names.5.1", "gene_ids.5.1"))
    tbl <- .naToBlank(tbl, omit=c("tx.dst", "gx.dst", "tx.pos.5", "tx.pos.3"))
    return(tbl)
}

.addCyto <- function(tbl, ann) {
    ## locus_id to cytoband
    l2c <- ann$cytoband.ovr
    setkey(l2c, locus_id)
    setkey(tbl, "locus_id.3.2")
    tbl <- l2c[tbl]
    setnames(tbl, c("locus_id", "cytoband_id"), c("locus_id.3.2", "cyt.3.2"))
    setkey(tbl, "locus_id.5.2")
    tbl <- l2c[tbl]
    setnames(tbl, c("locus_id", "cytoband_id"), c("locus_id.5.2", "cyt.5.2"))
    setkey(tbl, "locus_id.3.1")
    tbl <- l2c[tbl]
    setnames(tbl, c("locus_id", "cytoband_id"), c("locus_id.3.1", "cyt.3.1"))
    setkey(tbl, "locus_id.5.1")
    tbl <- l2c[tbl]
    setnames(tbl, c("locus_id", "cytoband_id"), c("locus_id.5.1", "cyt.5.1"))
    return(tbl)
}

.addExpressionSJ <- function(tbl, spl) {
    sjs <- data.table(as.data.frame(spl$sjs))
    sjs <- sjs[,.(chr.5=seqnames, pos.5=ifelse(strand=="+", start, end), str.5=strand,
                  chr.3=seqnames, pos.3=ifelse(strand=="+", end, start), str.3=strand,
                  nfrag)]
    sjs[, tlj.5:=sum(nfrag), by=.(chr.5, pos.5, str.5)]
    sjs[, tlj.3:=sum(nfrag), by=.(chr.3, pos.3, str.3)]
    ##
    tbl[,tot.sp.jnc.5:=0L]
    tbl[,tot.sp.jnc.3:=0L]
    setkey(tbl, chr.5, pos.5, str.5)
    setkey(sjs, chr.5, pos.5, str.5)
    tbl[sjs, tot.sp.jnc.5:=tlj.5]
    setkey(tbl, chr.3, pos.3, str.3)
    setkey(sjs, chr.3, pos.3, str.3)
    tbl[sjs, tot.sp.jnc.3:=tlj.3]
    return(tbl)
}

.addExpressionGeneLoci <- function(tbl, cts) {
    setkey(tbl, locus_id.5.1)
    tbl[cts, cpm.5:=cpm]
    setkey(tbl, locus_id.3.1)
    tbl[cts, cpm.3:=cpm]
    return(tbl)
}

.reportSV <- function(bun, spl, cts, ann, only.prot) {
    bpt <- copy(bun$bpt)
    ## total reads
    bpt[, tot.enc:=NA_integer_]
    bpt[, tot.enc:=sum(sum.jnc * (l1=="enc")), by=CHM.KEY]
    bpt[, tot.jnc:=NA_integer_]
    bpt[, tot.jnc:=sum(sum.jnc * (l1=="spn")), by=CHM.KEY]
    ## order and filter breakpoints
    bpt <- bpt[order((gmap.valid | mm2.valid), orf, d2a, l1=="spn", sum.jnc, decreasing=TRUE)]
    if (nrow(bpt) > 0) {
        bpt[, first:=(1:.N), by=CHM.KEY]
        if (ann$par$only.spn.bpt) {
            bpt <- bpt[(l1=="spn")]
        }
        if (ann$par$only.hx.bpt) {
            bpt <- bpt[(hq.bpt | (hi.bpt & first==1))]
        }
        bpt[, first:=NULL]
    }
    ## select best JNC per BPT
    jnc <- bun$jnc
    hqj.pen <- ifelse(jnc$hq.jnc, 1, 1e4)
    ovr.min <- pmin(jnc$ovr.5, jnc$ovr.3)
    ovr.len <- ifelse(jnc$ovr.5 < jnc$ovr.3, jnc$seg.5.len, jnc$seg.3.len)
    ovr.pen <- abs(ovr.min / ovr.len - 0.5)
    ovr.rnk <- order((hqj.pen * ovr.pen))
    jnc <- jnc[ovr.rnk]
    jnc <- jnc[, .SD[1], by=BPT.KEY]
    ## merge CHM, BPT, JNC
    setkeyv(bpt, BPT.KEY)
    setkeyv(jnc, BPT.KEY)
    tbl <- jnc[bpt]
    setkeyv(tbl, CHM.KEY)
    tbl <- .addNames(tbl, ann, only.prot)
    tbl <- .addCyto(tbl, ann)
    tbl <- .addExpressionSJ(tbl, spl)
    suppressWarnings({
        tbl[, tot.jnc.5:=sum(sum.jnc), by=.(chr.5, pos.5, str.5, l1=="spn")]
        tbl[, tot.jnc.3:=sum(sum.jnc), by=.(chr.3, pos.3, str.3, l1=="spn")]
    })
    tbl <- .addExpressionGeneLoci(tbl, cts)
    return(tbl)
}

.reportMini <- function(bun, spl, cts, ann, only.d2a, only.prot) {
    bpt <- copy(bun$bpt)
    if (only.d2a) {
        bpt <- bpt[(d2a)]
    }
    if (ann$par$only.spn.bpt) {
        bpt <- bpt[(type > -1)]
    }
    setkeyv(bpt, BPT.KEY)
    bpt <- .addNames(bpt, ann, only.prot)
    bpt <- .addCyto(bpt, ann)
    bpt <- .addExpressionSJ(bpt, spl)
    suppressWarnings({
        bpt[, tot.jnc.5:=sum(sum.jnc), by=.(chr.5, pos.5, str.5, l1=="spn")]
        bpt[, tot.jnc.3:=sum(sum.jnc), by=.(chr.3, pos.3, str.3, l1=="spn")]
    })
    bpt <- .addExpressionGeneLoci(bpt, cts)
    return(bpt)
}
