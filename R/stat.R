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

.alignment.stat <- function(jnc, bpt, spl, dir, ann) {
    logs <- dir$aln.log.fn
    nreads <- sum(sapply(logs, function(log) as.integer(str_match(readLines(log)[6], "\t(.*)")[,2])))
    if (length(dir$bw.fn)>0) {
        cov.tbl <- .coverage.tbl(dir$bw.fn, ann)
        cov <- cov.tbl[,.(avg.cov=mean(cov)),.(cut.cov, cut.len, npos)]
        spl.rte <- .splice.rate(dir$bw.fn, dir$sj.fn, ann)
        lng.rte <- cov.tbl[cut.len>=8,sum(cov)]/1e6
    } else {
        cov.tbl <- NULL
        cov <- NULL
        spl.rte <- NA_real_
        lng.rte <- NA_real_
    }
    ## alignment rate
    aln.rte <- mean(sapply(logs, function(log) as.numeric(str_sub(str_match(readLines(log)[10], "\t(.*)")[,2], end=-2))) / 100)
    nreads.eff <- (nreads * aln.rte) + nrow(jnc)
    ## on-target rate
    tgt.stat <- lapply(dir$cts.sum.fn, function(fn) {x <- fread(fn); unname(unlist(x[1,2,with=FALSE] / sum(x[,2,with=FALSE])))})
    names(tgt.stat) <- paste0(str_match(dir$cts.sum.fn, "_nsort_(.*)-cts.cts.summary")[,2], ".rte")
    aln.stat <- c(
        list(
            sum.frg = nreads,
            eff.frg = nreads.eff,
            sum.jnc = nrow(jnc),
            sum.bpt = nrow(bpt),
            sum.spl = sum(spl$locus.sjexp$nsfrag),
            sum.nsj = sum(spl$locus.sjexp$nsj),
            cov = cov,
            aln.rte = aln.rte,
            spl.rte = spl.rte,
            lng.rte = lng.rte
        ),
        tgt.stat
    )
    return(aln.stat)
}

.breakpoint.stat <- function(bpt, eff.frg) {
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
    art.rte <- bpt[eval(art.sel), sum(sum.jnc)] / eff.frg
    ## duplication rate
    dup.rte <- bpt[eval(dup.sel), 1-sum(unq.sum.jnc)/sum(sum.jnc)]
    ## ligation rate
    spn.dist.rte <- bpt[eval(spn.dist.nd.sel), sum(sum.jnc)] / bpt[eval(dist.nd.sel), sum(sum.jnc)]
    spn.prox.rte <- bpt[eval(spn.prox.nd.sel), sum(sum.jnc)] / bpt[eval(prox.nd.sel), sum(sum.jnc)]
    spn.dist.lig <- bpt[eval(dist.lig.nd.sel), sum(sum.jnc)]
    spn.prox.lig <- bpt[eval(prox.lig.nd.sel), sum(sum.jnc)]
    dist.lig <- spn.dist.lig / spn.dist.rte
    prox.lig <- spn.prox.lig / spn.prox.rte
    lig.rte <- (dist.lig+prox.lig) / eff.frg
    ## trans-splice log fold-change
    spn.dist <- bpt[eval(spn.dist.sel), sum(sum.jnc)]
    enc.dist <- bpt[eval(enc.dist.sel), sum(sum.jnc)]
    spn.dist.gtag <- bpt[eval(spn.dist.gtag.sel), sum(sum.jnc)]
    enc.dist.gtag <- bpt[eval(enc.dist.gtag.sel), sum(sum.jnc)]
    gtag.exp.rte <- (enc.dist.gtag/enc.dist)
    gtag.obs.rte <- (spn.dist.gtag/spn.dist)
    ts.lfc <- log2(gtag.obs.rte / gtag.exp.rte)
    ## trans-splice rate
    gtag.exp.cnt <- spn.dist * gtag.exp.rte
    gtag.obs.cnt <- spn.dist.gtag
    ts.rte.lo <- (bpt[eval(ts.sel), sum(sum.jnc)] / spn.dist.rte) / eff.frg
    ts.rte.hi <- (max(gtag.obs.cnt - gtag.exp.cnt, 0) / spn.dist.rte) / eff.frg
    ## back-splice rate
    spn.prox.d2a <- bpt[eval(spn.prox.d2a.dup.sel), sum(sum.jnc)]
    prox.d2a <- spn.prox.d2a / spn.prox.rte
    bs.rte <- prox.d2a / eff.frg
    ## stem-loop rate
    sl.rte <- bpt[eval(sl.sel), sum(sum.jnc)] / eff.frg
    bpt.surv <- bpt[,.(
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

.locus.stat <- function(bpt, ann) {
    tmp.x <- data.table(
        locus_id=factor(),
        enc_dist_lig=integer(),
        enc_prox_lig=integer(),
        spn_dist_lig=integer(),
        spn_dist_mot=integer(),
        spn_prox_lig=integer(),
        spn_prox_mot=integer()
    )
    tmp.5 <- bpt[,.(N=sum(sum.jnc)), by=.(locus_id=locus_id.5.1, l1, l2, l3, d2a)]
    tmp.5 <- dcast.data.table(tmp.5, locus_id~l1+l2+ifelse(d2a, "mot", "lig"), value.var = "N", fun.aggregate = sum)
    tmp.5 <- rbind(tmp.5, tmp.x, fill=TRUE)[,colnames(tmp.x),with=FALSE]
    setnames(tmp.5, c("locus_id", paste0(colnames(tmp.5)[-1], ".5")))
    tmp.3 <- bpt[,.(N=sum(sum.jnc)), by=.(locus_id=locus_id.3.1, l1, l2, l3, d2a)]
    tmp.3 <- dcast.data.table(tmp.3, locus_id~l1+l2+ifelse(d2a, "mot", "lig"), value.var = "N", fun.aggregate = sum)
    tmp.3 <- rbind(tmp.3, tmp.x, fill=TRUE)[,colnames(tmp.x),with=FALSE]
    setnames(tmp.3, c("locus_id", paste0(colnames(tmp.3)[-1], ".3")))
    setkey(tmp.5, locus_id)
    setkey(tmp.3, locus_id)
    tmp.53 <- merge(tmp.5, tmp.3, all=TRUE)
    tmp.53 <- tmp.53[ann$loci$locus_id]
    tmp.53 <- .naTo0(tmp.53)
    setkey(tmp.53, locus_id)
    return(list(loc.surv=tmp.53))
}

#' @export
makeStat <- function(jnc, bpt, spl, dir, ann) {
    aln.stat <- .alignment.stat(jnc, bpt, spl, dir, ann)
    bpt.stat <- .breakpoint.stat(bpt, aln.stat$eff.frg)
    loc.stat <- .locus.stat(bpt, ann)
    log.stat <- list(version=packageVersion("codac"), beg.time=NULL, end.time=NULL, gc=NULL)
    stat <- c(aln.stat, bpt.stat, loc.stat)
    stat$log <- log.stat
    stat$dir <- dir
    return(stat)
}
