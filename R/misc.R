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

.prepareTxs <- function(bun, lgt) {
    lids <- unique(c(
        as.character(bun$bpt$locus_id.5.1),
        as.character(bun$bpt$locus_id.3.1),
        as.character(bun$bpt$locus_id.5.2),
        as.character(bun$bpt$locus_id.3.2)
    ))
    txs <- lgt[lids,nomatch=0]
    txs[,mrna:=list(as.list(DNAStringSet(tx.seq)))]
    txs[,prot:=list(as.list(suppressWarnings(translate(subseq(DNAStringSet(tx.seq), cds1, cdsN)))))]
    txs <-  txs[,.(locus_id, transcript_id, mrna, prot, cds1, cdsN, tags)]
    return(txs)
}

.splice.rate <- function(bw.fn, sj.fn, ann) {
    ## SJs in protein coding transcripts
    prot.tx.id <- ann$transcripts[ann$transcripts$cdsN>0]$transcript_id
    tmp <- as.data.table(ann$features[ann$features$type=="I"])
    setkey(tmp, split_id)
    sj.ref <- unique(tmp[prot.tx.id,.(chr=seqnames, beg=start, end=end, str=strand), nomatch=0])
    ## SJs observed in sample
    sj.tbl <- setNames(fread(paste(system2("zcat", sj.fn, stdout=TRUE), collapse="\n"), showProgress=FALSE),
                       c("chr", "beg", "end", "str", "mot", "known", "n.unq", "n.mult", "ovr"))
    sj.tbl[,str:=ifelse(str==1, "+", "-")]
    setkey(sj.tbl, chr, beg, end, str)
    setkey(sj.ref, chr, beg, end, str)
    sj.tbl.flt <- sj.tbl[sj.ref, nomatch=0]
    sj.gr.flt.beg0 <- with(sj.tbl.flt, GRanges(chr, IRanges(beg, beg), "*"))
    sj.gr.flt.end0 <- with(sj.tbl.flt, GRanges(chr, IRanges(end, end), "*"))
    sj.gr.flt.0 <- c(unique(sj.gr.flt.beg0), unique(sj.gr.flt.end0))
    sj.gr.flt.beg1 <- with(sj.tbl.flt, GRanges(chr, IRanges(beg-1, beg-1), "*"))
    sj.gr.flt.end1 <- with(sj.tbl.flt, GRanges(chr, IRanges(end+1, end+1), "*"))
    sj.gr.flt.1 <- c(unique(sj.gr.flt.beg1), unique(sj.gr.flt.end1))
    tmp0 <- unname(unlist(import(bw.fn, which = sj.gr.flt.0, as="NumericList")))
    tmp1 <- unname(unlist(import(bw.fn, which = sj.gr.flt.1, as="NumericList")))
    res <- 1 - sum(tmp0)/(sum(tmp0)+sum(tmp1))
    return(res)
}

.coverage.tbl <- function(bw.fn, ann, n=10) {
    vseq <- Vectorize(function(from, to, lo) floor(seq.int(from, to, length.out=lo)), vectorize.args=c("from", "to"))
    cbt <- split(ann$features[ann$features$type=="C"], ann$features[ann$features$type=="C"]$split_id)
    tmp <- data.table(transcript_id=names(cbt), t(vseq(1, sum(width(cbt)), lo=n)))
    tmp <- melt(tmp, id.vars=c("transcript_id"))
    tmp[,variable:=as.integer(str_sub(variable, -1))] # remove V
    setkey(tmp, transcript_id, variable)
    hit <- unlist(pmapFromTranscripts(IRanges(tmp$value, width=1), cbt[tmp$transcript_id]))
    hit <- hit[hit$hit]
    cov <- unname(unlist(import.bw(bw.fn, which=hit, as="NumericList")))
    tbl <- tmp[,.(transcript_id, pos=value, npos=variable, cov)]
    tbl[,avg.cov:=mean(cov), by=transcript_id]
    tbl[,med.cov:=median(cov), by=transcript_id]
    tbl[,cds.len:=max(pos), by=transcript_id]
    tbl[,cut.cov:=.bincode(avg.cov, quantile(avg.cov, seq(0,1,length.out=n+1)), include.lowest=TRUE)-1]
    tbl[,cut.len:=.bincode(cds.len, quantile(cds.len, seq(0,1,length.out=n+1)), include.lowest=TRUE)-1]
    return(tbl)
}

.alignment.stat <- function(bun, ann) {
    logs <- bun$dir$aln.log.fn
    nreads <- sum(sapply(logs, function(log) as.integer(str_match(readLines(log)[6], "\t(.*)")[,2])))
    if (length(bun$dir$bw.fn)>0) {
        cov.tbl <- .coverage.tbl(bun$dir$bw.fn, ann)
        cov <- cov.tbl[,.(avg.cov=mean(cov)),.(cut.cov, cut.len, npos)]
        spl.rte <- .splice.rate(bun$dir$bw.fn, bun$dir$sj.fn, ann)
        lng.rte <- cov.tbl[cut.len>=8,sum(cov)]/1e6
    } else {
        cov.tbl <- NULL
        cov <- NULL
        spl.rte <- NA_real_
        lng.rte <- NA_real_
    }
    ## alignment rate
    aln.rte <- mean(sapply(logs, function(log) as.numeric(str_sub(str_match(readLines(log)[10], "\t(.*)")[,2], end=-2))) / 100)
    nreads.eff <- (nreads * aln.rte) + nrow(bun$jnc)
    ## on-target rate
    tgt.stat <- lapply(bun$dir$cts.sum.fn, function(fn) {x <- fread(fn); unname(unlist(x[1,2,with=FALSE] / sum(x[,2,with=FALSE])))})
    names(tgt.stat) <- paste0(str_match(bun$dir$cts.sum.fn, "_nsort_(.*)-cts.cts.summary")[,2], ".rte")
    aln.stat <- c(
        list(
            sum.frg = nreads,
            eff.frg = nreads.eff,
            sum.jnc = nrow(bun$jnc),
            sum.bpt = nrow(bun$bpt),
            sum.spl = sum(bun$spl$locus.sjexp$nsfrag),
            sum.nsj = sum(bun$spl$locus.sjexp$nsj),
            cov = cov,
            aln.rte = aln.rte,
            spl.rte = spl.rte,
            lng.rte = lng.rte
        ),
        tgt.stat
    )
    return(aln.stat)
}

.breakpoint.stat <- function(bun, eff.frg) {
    ## dt selectors
    sl.sel <- quote(l3=="sl")
    ts.sel <- quote(l3=="ts")
    art.sel <- quote(art==TRUE)
    dup.sel <- quote(art.5=="nil" & art.3=="nil" & l2=="dist" & l3=="nd" & !d2a)
    spn.dist.nd.sel <- quote(l1=="spn" & l2=="dist" & l3=="nd")
    spn.prox.nd.sel <- quote(l1=="spn" & l2=="prox" & l3=="nd")
    dist.nd.sel <- quote(!art & l2=="dist" & l3=="nd")
    prox.nd.sel <- quote(!art & l2=="prox" & l3=="nd")
    dist.lig.nd.sel <- quote(l1=="spn" & l2=="dist" & l3=="nd" & !d2a)
    prox.lig.nd.sel <- quote(l1=="spn" & l2=="prox" & l3=="nd" & !d2a)
    spn.dist.sel <- quote(l1=="spn" & l2=="dist")
    enc.dist.sel <- quote(l1=="enc" & l2=="dist")
    spn.dist.gtag.sel <- quote(l1=="spn" & l2=="dist" & ((str_sub(mot.5,7,8)=="GT") & (str_sub(mot.3,7,8)=="AG")))
    enc.dist.gtag.sel <- quote(l1=="enc" & l2=="dist" & ((str_sub(mot.5,7,8)=="GT") & (str_sub(mot.3,7,8)=="AG")))
    spn.prox.d2a.dup.sel <- quote(l1=="spn" & l2=="prox" & d2a & topo=="dup")
    ## artifact rate
    art.rte <- bun$bpt[eval(art.sel), sum(sum.jnc)] / eff.frg
    ## duplication rate
    dup.rte <- bun$bpt[eval(dup.sel), 1-sum(unq.sum.jnc)/sum(sum.jnc)]
    ## ligation rate
    spn.dist.rte <- bun$bpt[eval(spn.dist.nd.sel), sum(sum.jnc)] / bun$bpt[eval(dist.nd.sel), sum(sum.jnc)]
    spn.prox.rte <- bun$bpt[eval(spn.prox.nd.sel), sum(sum.jnc)] / bun$bpt[eval(prox.nd.sel), sum(sum.jnc)]
    spn.dist.lig <- bun$bpt[eval(dist.lig.nd.sel), sum(sum.jnc)]
    spn.prox.lig <- bun$bpt[eval(prox.lig.nd.sel), sum(sum.jnc)]
    dist.lig <- spn.dist.lig / spn.dist.rte
    prox.lig <- spn.prox.lig / spn.prox.rte
    lig.rte <- (dist.lig+prox.lig) / eff.frg
    ## trans-splice log fold-change
    spn.dist <- bun$bpt[eval(spn.dist.sel), sum(sum.jnc)]
    enc.dist <- bun$bpt[eval(enc.dist.sel), sum(sum.jnc)]
    spn.dist.gtag <- bun$bpt[eval(spn.dist.gtag.sel), sum(sum.jnc)]
    enc.dist.gtag <- bun$bpt[eval(enc.dist.gtag.sel), sum(sum.jnc)]
    gtag.exp.rte <- (enc.dist.gtag/enc.dist)
    gtag.obs.rte <- (spn.dist.gtag/spn.dist)
    ts.lfc <- log2(gtag.obs.rte / gtag.exp.rte)
    ## trans-splice rate
    gtag.exp.cnt <- spn.dist * gtag.exp.rte
    gtag.obs.cnt <- spn.dist.gtag
    ts.rte.lo <- (bun$bpt[eval(ts.sel), sum(sum.jnc)] / spn.dist.rte) / eff.frg
    ts.rte.hi <- (max(gtag.obs.cnt - gtag.exp.cnt, 0) / spn.dist.rte) / eff.frg
    ## back-splice rate
    spn.prox.d2a <- bun$bpt[eval(spn.prox.d2a.dup.sel), sum(sum.jnc)]
    prox.d2a <- spn.prox.d2a / spn.prox.rte
    bs.rte <- prox.d2a / eff.frg
    ## stem-loop rate
    sl.rte <- bun$bpt[eval(sl.sel), sum(sum.jnc)] / eff.frg
    bpt.surv <- bun$bpt[,.(
        n=.N,
        hq.sum.jnc=mean(hq.sum.jnc), unq.sum.jnc=mean(unq.sum.jnc), sum.jnc=mean(sum.jnc),
        avg.err.5=mean(avg.err.5), avg.err.3=mean(avg.err.3),
        avg.low.5=mean(avg.low.5), avg.low.3=mean(avg.low.3),
        max.ovr.5=mean(max.ovr.5), max.ovr.3=mean(max.ovr.3),
        unq.ovr.5=mean(unq.ovr.5), unq.ovr.3=mean(unq.ovr.3),
        orf=mean(orf), spn=mean(type>-1),
        hq.bpt=mean(hq.bpt), hi.bpt=mean(hi.bpt)
    ), by=.(
           l1, l2, l3, topo, dst, d2a,
           gene.5=(cls.5 != "."), gene.3=(cls.3 != "."),
           ok.5=(art.5=="nil"), ok.3=(art.3=="nil")
       )
    ]
    res <- list(art.rte=art.rte, dup.rte=dup.rte, lig.rte=lig.rte,
                ts.rte.lo=ts.rte.lo, ts.rte.hi=ts.rte.hi, ts.lfc=ts.lfc,
                bs.rte=bs.rte, sl.rte=sl.rte,
                s2e.prox.rte=spn.prox.rte, s2e.dist.rte=spn.dist.rte,
                bpt.surv=bpt.surv)
    return(res)
}
