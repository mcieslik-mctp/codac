.collapseBreakpoints <- function(jnc) {
    bpt <- jnc[, .(
      sum.jnc=.N,
      hq.sum.jnc=sum(hq.jnc, na.rm=TRUE),
      unq.sum.jnc=sum(unq.jnc, na.rm=TRUE),
      avg.err.5=mean(err.5, na.rm=TRUE),
      avg.err.3=mean(err.3, na.rm=TRUE),
      avg.low.5=mean(low.5, na.rm=TRUE),
      avg.low.3=mean(low.3, na.rm=TRUE),
      max.ovr.5=max(ovr.5, na.rm=TRUE),
      max.ovr.3=max(ovr.3, na.rm=TRUE),
      unq.ovr.5=length(.Internal(unique(ovr.5[!is.na(ovr.5)], FALSE, FALSE, NA))),
      unq.ovr.3=length(.Internal(unique(ovr.3[!is.na(ovr.3)], FALSE, FALSE, NA)))
    ), by=BPT.KEY]
    setkeyv(bpt, BPT.KEY)
    return(bpt)
}

.classifyBreakpoints <- function(bpt, ann) {
    ## breakpoint class
    cls <- ann$features[ann$features$type %in% c("3", "5", "C", "E", "I", ".")]
    bpt.5s <- split(with(bpt, GRanges(chr.5, IRanges(pos.5-2, pos.5+2), str.5)), factor(ceiling(seq_len(nrow(bpt)) / 100000)))
    bpt.3s <- split(with(bpt, GRanges(chr.3, IRanges(pos.3-2, pos.3+2), str.3)), factor(ceiling(seq_len(nrow(bpt)) / 100000)))
    bpt[,cls.5:=unlist(unname(lapply(bpt.5s, function(b5) {cls[findOverlaps(b5, cls, select="first")]$type})))]
    bpt[,cls.3:=unlist(unname(lapply(bpt.3s, function(b3) {cls[findOverlaps(b3, cls, select="first")]$type})))]
    ss <- ann$features[ann$features$type %in% c("D", "A", ".")]
    ss.5o <- unlist(GRangesList(unname(lapply(bpt.5s, function(b5) {ss[findOverlaps(b5, ss, select="first")]}))))
    ss.3o <- unlist(GRangesList(unname(lapply(bpt.3s, function(b3) {ss[findOverlaps(b3, ss, select="first")]}))))
    bpt[,ss.5:=ss.5o$type]
    bpt[,ss.3:=ss.3o$type]
    bpt[,pos.5.fix:=ifelse(width(ss.5o)==1, start(ss.5o), bpt$pos.5)]
    bpt[,pos.3.fix:=ifelse(width(ss.3o)==1, start(ss.3o), bpt$pos.3)]
    rm(ss, ss.5o, ss.3o); gc() ## memory optimization
    ## splice class
    tmp <- cls[cls$type==".", NULL]
    tmp$frame <- 3
    tmp$type <- "."
    frm <- c(ann$frames$frame.da, tmp)
    bpt[,frm.5:=unlist(unname(lapply(bpt.5s, function(b5) {frm[findOverlaps(b5, frm, select="first")]$frame})))]
    bpt[,frm.3:=unlist(unname(lapply(bpt.3s, function(b3) {frm[findOverlaps(b3, frm, select="first")]$frame})))]
    rm(cls, tmp, frm, bpt.5s, bpt.3s); gc() ## memory optimization
    ## frame prediction
    bpt[, orf:=
          (
            (type > -1) & (frm.5 == frm.3) & (frm.5 < 3) & (frm.3 < 3)
          ) |
          (
            (type > -1) & (cls.5 %in% c("5", "E")) & (cls.3 %in% c("5", "C"))
          )   
        ]
    bpt[, d2a:=
          (
            (type == 1)
            |
            (type == 0) & (ss.5 == "D" & ss.3 == "A")
          )
        ]
    bpt[,pos.5.fix:=NULL]
    bpt[,pos.3.fix:=NULL]
    return(bpt)
}

.motifBreakpoints <- function(bpt) {
    fl <- 7
    rng.5 <- GRanges(bpt$chr.5, IRanges(bpt$pos.5-fl, bpt$pos.5+fl), bpt$str.5)
    rng.3 <- GRanges(bpt$chr.3, IRanges(bpt$pos.3-fl, bpt$pos.3+fl), bpt$str.3)
    bpt[,mot.5:=str_sub(getSeq(BSgenome.Hsapiens.UCSC.hg38, rng.5, as.character=TRUE), 1+1, 2*fl+1)]
    bpt[,mot.3:=str_sub(getSeq(BSgenome.Hsapiens.UCSC.hg38, rng.3, as.character=TRUE), 1,   2*fl  )]
    return(bpt)
}

.distanceBreakpoints <- function(bpt, ann) {
    batch <- 100000
    ebt <- split(ann$features[ann$features$type=="E"], ann$features[ann$features$type=="E"]$split_id)
    batch.ids <- factor(ceiling(seq_len(nrow(bpt)) / batch))
    bpt.5s <- split(with(bpt, GRanges(chr.5, IRanges(pos.5 + ifelse(str.5=="+", -1, +1), width=1), str.5)), batch.ids)
    bpt.3s <- split(with(bpt, GRanges(chr.3, IRanges(pos.3 + ifelse(str.3=="+", +1, -1), width=1), str.3)), batch.ids)
    tmp.5 <- unlist(GRangesList(unname(lapply(names(bpt.5s), function(b5n) {
        b5 <- bpt.5s[[b5n]]
        b5m <- mapToTranscripts(b5, ebt, ignore.strand=TRUE)
        b5m$xHits <- b5m$xHits + ((as.integer(b5n) - 1) * batch)
        return(b5m)
    }))))
    ## bpt.5 <- GRanges(bpt$chr.5, IRanges(bpt$pos.5, width=1), bpt$str.5)
    ## tmp.5 <- mapToTranscripts(bpt.5, ebt, ignore.strand=TRUE)
    mcols(tmp.5)$idx.5 <- 1:length(tmp.5)
    tmp.3 <- unlist(GRangesList(unname(lapply(names(bpt.3s), function(b3n) {
        b3 <- bpt.3s[[b3n]]
        b3m <- mapToTranscripts(b3, ebt, ignore.strand=TRUE)
        b3m$xHits <- b3m$xHits + ((as.integer(b3n) - 1) * batch)
        return(b3m)
    }))))
    ## bpt.3 <- GRanges(bpt$chr.3, IRanges(bpt$pos.3, width=1), bpt$str.3)
    ## tmp.3 <- mapToTranscripts(bpt.3, ebt, ignore.strand=TRUE)
    mcols(tmp.3)$idx.3 <- 1:length(tmp.3)
    tmp.5.dt <- data.table(as.data.frame(mcols(tmp.5)))
    setkeyv(tmp.5.dt, c("xHits", "transcriptsHits"))
    tmp.3.dt <- data.table(as.data.frame(mcols(tmp.3)))
    setkeyv(tmp.3.dt, c("xHits", "transcriptsHits"))
    hits <- merge(tmp.5.dt, tmp.3.dt)
    hits$tx.pos.5 <- start(tmp.5[hits$idx.5])
    hits$tx.pos.3 <- start(tmp.3[hits$idx.3])
    hits[, distance:=tx.pos.5-tx.pos.3]
    hits[, whichdst:=which.min(abs(distance)), by=xHits]
    ## minimum possible transcriptome distance
    tmp <- hits[,.(
      tx.dst=distance[whichdst],
      tx.pos.5=tx.pos.5[whichdst],
      tx.pos.3=tx.pos.3[whichdst]
    ), by=xHits]
    bpt$tx.dst <- NA_integer_
    bpt$tx.pos.5 <- NA_integer_
    bpt$tx.pos.3 <- NA_integer_
    bpt[tmp$xHits, ":="(
      tx.dst=tmp$tx.dst * ifelse(str.5 == "+" & str.3 == "+", -1L, 1L),
      tx.pos.5=tmp$tx.pos.5,
      tx.pos.3=tmp$tx.pos.3
    )]
    ## genomic distance
    bpt$gx.dst <- NA_integer_
    bpt[chr.5==chr.3, gx.dst:=pos.5-pos.3]
    bpt[, bp.dst:=pmin(abs(ifelse(is.na(tx.dst), MAXINT, tx.dst)),
                       abs(ifelse(is.na(gx.dst), MAXINT, gx.dst)))]
    return(bpt)
}

.associateBreakpoints <- function(bpt, ann) {
    ## associate
    bpt.5 <- GRanges(bpt$chr.5, IRanges(bpt$pos.5, width=1),bpt$str.5)
    bpt.3 <- GRanges(bpt$chr.3, IRanges(bpt$pos.3, width=1),bpt$str.3)
    hits.5 <- suppressWarnings(findOverlaps(bpt.5, ann$loci))
    hits.3 <- suppressWarnings(findOverlaps(bpt.3, ann$loci))
    rm(bpt.5, bpt.3); gc() ## memory optimization
    locus_id.5.1 <- mcols(ann$loci[subjectHits(hits.5)])$locus_id
    locus_id.3.1 <- mcols(ann$loci[subjectHits(hits.3)])$locus_id
    rm(hits.5, hits.3); gc() ## memory optimization
    ## associate reverse and complement
    bpt.5 <- GRanges(bpt$chr.3, IRanges(bpt$pos.3, width=1), ifelse(bpt$str.3=="+","-","+"))
    bpt.3 <- GRanges(bpt$chr.5, IRanges(bpt$pos.5, width=1), ifelse(bpt$str.5=="+","-","+"))
    hits.5 <- suppressWarnings(findOverlaps(bpt.5, ann$loci))
    hits.3 <- suppressWarnings(findOverlaps(bpt.3, ann$loci))
    rm(bpt.5, bpt.3); gc() ## memory optimization
    locus_id.5.2 <- mcols(ann$loci[subjectHits(hits.5)])$locus_id
    locus_id.3.2 <- mcols(ann$loci[subjectHits(hits.3)])$locus_id
    rm(hits.5, hits.3); gc() ## memory optimization
    ## undo reverse and assign
    need.flip <- (locus_id.5.1 == locus_id.5.2) & (locus_id.3.1 == locus_id.3.2)
    bpt$locus_id.5.1 <- factor(locus_id.5.1, levels=ann$loci$locus_id)
    bpt$locus_id.3.1 <- factor(locus_id.3.1, levels=ann$loci$locus_id)
    bpt$locus_id.5.2 <- factor(ifelse(need.flip, locus_id.3.2, locus_id.5.2), levels=ann$loci$locus_id)
    bpt$locus_id.3.2 <- factor(ifelse(need.flip, locus_id.5.2, locus_id.3.2), levels=ann$loci$locus_id)
    return(bpt)
}

.topologyBreakpoints <- function(bpt, ann) {
    ## number of spanning breakpoints
    bpt[, ":="(sum.bpt=sum(type>-1)), by=CHM.KEY]
    ## cytoband
    setkeyv(ann$locus.ovr, c("locus_id.5", "locus_id.3"))
    setkeyv(bpt, c("locus_id.5.1", "locus_id.3.1"))
    loc.5 <- ann$loci[as.integer(bpt$locus_id.5.1)]
    str.5 <- strand(loc.5)
    loc.3 <- ann$loci[as.integer(bpt$locus_id.3.1)]
    str.3 <- strand(loc.3)
    cb_id.5 <- ann$cytoband.ovr[as.integer(bpt$locus_id.5.1)]$cytoband_id
    cb_id.3 <- ann$cytoband.ovr[as.integer(bpt$locus_id.3.1)]$cytoband_id
    cb_idx.5 <- ann$cytoband2idx[cb_id.5]$idx
    cb_idx.3 <- ann$cytoband2idx[cb_id.3]$idx
    has.cb <- (!is.na(cb_idx.5)) & (!is.na(cb_idx.3))
    ## distances
    same.chr <- as.logical(seqnames(loc.5) == seqnames(loc.3))
    same.arm <- same.chr & has.cb & (mcols(ann$cytobands)$arm[cb_idx.5] ==
                                     mcols(ann$cytobands)$arm[cb_idx.3])
    same.cyt <- has.cb & (cb_id.5 == cb_id.3)
    same.frg <- bpt$bp.dst < ann$par$bpt.dst.frg
    ## update bpt
    bpt$dst <- factor("gme", levels=c("frg", "loc", "adj", "cyt", "arm", "chr", "gme"), ordered=TRUE)         
    bpt[same.chr, dst:="chr"]
    bpt[same.arm, dst:="arm"]
    bpt[same.cyt, dst:="cyt"]
    bpt[ann$locus.adj, dst:="adj"]
    setkeyv(bpt, c("locus_id.5.1", "locus_id.5.2"))
    bpt[ann$locus.adj, dst:="adj"]
    setkeyv(bpt, c("locus_id.5.1", "locus_id.3.1"))
    bpt[(locus_id.5.1 == locus_id.3.1) | (locus_id.5.1 == locus_id.5.2) |
        (locus_id.3.1 == locus_id.3.2) | ((locus_id.5.1 == locus_id.3.2) &
                                          (locus_id.3.1 == locus_id.5.2)), dst:="loc"]
    bpt[same.frg, dst:="frg"]
    ## inversions, duplications, translocations, deletions
    sme.str <- as.logical(str.5 == str.3)
    pls.str <- as.logical(same.chr & sme.str & (str.5 == "+"))
    mns.str <- as.logical(same.chr & sme.str & (str.5 == "-"))
    inv <- !sme.str
    dup <- (pls.str & loc.5 > loc.3) | (mns.str & loc.5 < loc.3) | (sme.str & loc.5==loc.3)
    bpt[, topo:=factor(ifelse(dst!="gme",
                       ifelse(inv, "inv",
                       ifelse(dup, "dup", "del")), "tloc"),
            levels=c("inv", "dup", "del", "tloc"), ordered=TRUE)]
    return(bpt)
}

.defineHQbreakpoints <- function(bpt, ann) {
    par <- ann$par
    pen <- par$bpt.art.pen
    bpt[,hq.bpt:=FALSE]
    bpt[
        ( ## not known alignment artifact
            art.5 != "aln" & art.3 != "aln"
        )
        &
        ( ## read support
          ( ## HQ reads
            (((art.5 == "nil" & art.3 == "nil") | d2a) & (hq.sum.jnc >= par$bpt.sum.hq.jnc))
            |
            (hq.sum.jnc >= par$bpt.sum.hq.jnc * pen)
          )
          &
          ( ## unique reads
            (((art.5 == "nil" & art.3 == "nil") | d2a) & (unq.sum.jnc >= par$bpt.sum.unq.jnc))
            |
            (unq.sum.jnc >= par$bpt.sum.unq.jnc * pen)
          )
          &
          ( ## total reads
            (((art.5 == "nil" & art.3 == "nil") | d2a) & (sum.jnc >= par$bpt.sum.jnc))
            |
            (sum.jnc >= par$bpt.sum.jnc * pen)
          )
        )
        &
        ( ## breakpoint topology and distance
            (topo %in% c("inv", "tloc")) | (bp.dst > par$bpt.dst.spl)
        )
        &
        ( ## breakpoint quality
          (((art.5 == "nil" | d2a) & (unq.ovr.5 >= par$bpt.unq.ovr)) | (unq.ovr.5 >= par$bpt.unq.ovr * pen))
          &
          (((art.3 == "nil" | d2a) & (unq.ovr.3 >= par$bpt.unq.ovr)) | (unq.ovr.3 >= par$bpt.unq.ovr * pen))
          &
          (((art.5 == "nil" | d2a) & (max.ovr.5 >= par$bpt.max.ovr)) | (max.ovr.5 >= par$bpt.max.ovr * pen))
          &
          (((art.3 == "nil" | d2a) & (max.ovr.3 >= par$bpt.max.ovr)) | (max.ovr.3 >= par$bpt.max.ovr * pen))
          &
          (((art.5 == "nil" | d2a) & (avg.err.5 <= par$bpt.err)) | (avg.err.5 <= par$bpt.err / pen))
          &
          (((art.3 == "nil" | d2a) & (avg.err.3 <= par$bpt.err)) | (avg.err.3 <= par$bpt.err / pen))
          &
          (((art.5 == "nil" | d2a) & (avg.low.5 <= par$bpt.low)) | (avg.low.5 <= par$bpt.low / pen))
          &
          (((art.3 == "nil" | d2a) & (avg.low.3 <= par$bpt.low)) | (avg.low.3 <= par$bpt.low / pen))
        ),
        hq.bpt:=TRUE
        ]
    return(bpt)
}

.defineHIbreakpoints <- function(bpt, ann) {
    hi.art <- c("aln", "gtb", "cmt")
    hi.lid.goi <- .l2g.full(ann)[(goi)]$locus_id
    hi.lid.loi <- ann$loci[ann$loci %over% ann$loi]$locus_id
    hi.lid <- unique(c(hi.lid.goi, hi.lid.loi))
    bpt[,hi.bpt:=FALSE]
    bpt[
        ( ## not known alignment artifact or blacklisted gene
            !(art.5 %in% hi.art) & !(art.3 %in% hi.art)
        )
        &
        ( ## high-importance locus
            (locus_id.5.1 %in% hi.lid) |
            (locus_id.3.1 %in% hi.lid) |
            (locus_id.5.2 %in% hi.lid) |
            (locus_id.3.2 %in% hi.lid)
        ),
        hi.bpt:=
            (
                (sum(unq.sum.jnc) > ann$par$chm.sum.hi.jnc) &
                (all(dst %in% c("adj", "cyt", "arm", "chr", "gme"))) &
                (any(type>-1 & d2a)) # at least one spliced spanning read
            ), by=CHM.KEY]
    return(bpt)
}

.defineHRbreakpoints <- function(bpt, ann) {
    ## recurrent d2a breakpoints
    hr.art <- c("rpx", "ret", "cmt")
    bpt[,hr.bpt:=FALSE]
    bpt[l1=="spn" & l2=="dist" & d2a & !(art.5 %in% hr.art) & !(art.3 %in% hr.art), ":="(
        rec.5=.N,
        unq.rec.5=1L
    ), by=.(chr.5, pos.5, str.5)]
    bpt[rec.5>1, unq.rec.5:=uniqueN(.SD), by=.(chr.5, pos.5, str.5), .SDcols=c("locus_id.3.1", "locus_id.5.2")]
    bpt[l1=="spn" & l2=="dist" & d2a & !(art.5 %in% hr.art) & !(art.3 %in% hr.art), ":="(
        rec.3=.N,
        unq.rec.3=1L
    ), by=.(chr.3, pos.3, str.3)]
    bpt[rec.3>1, unq.rec.3:=uniqueN(.SD), by=.(chr.3, pos.3, str.3), .SDcols=c("locus_id.5.1", "locus_id.3.2")]
    rec <- ann$par$bpt.rec.ts
    tmp <- bpt[(unq.rec.5 >= rec | unq.rec.3 >= rec), .(l1, chr.5, pos.5, str.5, chr.3, pos.3, str.3)]
    setkey(tmp, l1, chr.5, pos.5, str.5)
    setkey(bpt, l1, chr.5, pos.5, str.5)
    bpt[unique(tmp, by=key(tmp)), hr.bpt:=TRUE]
    setkey(tmp, l1, chr.3, pos.3, str.3)
    setkey(bpt, l1, chr.3, pos.3, str.3)
    bpt[unique(tmp, by=key(tmp)), hr.bpt:=TRUE]
    return(bpt)
}

.defineHCbreakpoints <- function(bpt, ann) {
    bpt[,hc.bpt:=FALSE]
    bpt[(l1=="spn" & l2=="prox" & (topo == "inv" | !d2a)), hc.bpt:=(.N >= ann$par$bpt.rec.chm),
        by=.(locus_id.5.1, locus_id.3.1, locus_id.5.2, locus_id.3.2)]
    return(bpt)
}

.defineLevels <- function(bpt, ann) {
    ## Level 1 (artifact, spanning, encompassing)
    bpt[,l1:=ifelse(type>-1, "spn", "enc")]
    ## Level 2 (proximal, distal)
    bpt[,l2:=ifelse(
    (
        ## inversion
        ((topo == "inv") & (dst > "loc")) |
        ## duplication
        ((topo == "dup") & (dst > "adj")) |
        ## deletion
        ((topo == "del") & (dst > "adj")) |
        ##
        ((topo == "tloc"))
    ), "dist", "prox")]
    ## Level 3
    bpt[,l3:=factor("nd", levels=c("nd", "sl", "bs", "ts", "sv"), ordered=TRUE)]
    return(bpt)
}

.defineHXBreakpoints <- function(bpt, ann) {
    bpt <- .defineHQbreakpoints(bpt, ann) # high-quality SV
    bpt <- .defineHIbreakpoints(bpt, ann) # high-importance SV
    bpt <- .defineHRbreakpoints(bpt, ann) # recurrent TS
    bpt <- .defineHCbreakpoints(bpt, ann) # clustered SL
    return(bpt)
}

.linkBreakpoints <- function(bpt, ann) {
    ## not adjacent read-through, not adjacent tandem duplication, not within locus
    hr.art <- c("rpx", "ret", "cmt")
    types=list(
        "sl"=quote((hc.bpt & (topo == "inv" | !d2a) & !art ) & l1=="spn" & l2=="prox"),
        "bs"=quote((         (topo == "dup" &  d2a) & !art ) & l1=="spn" & l2=="prox"),
        "ts"=quote((hr.bpt                                 ) & l1=="spn" & l2=="dist"),
        "sv"=quote((hi.bpt |                (hq.bpt & !art)) & l1=="spn" & l2=="dist")
    )
    for (name in names(types)) {
        sel <- types[[name]]
        sel.link <- paste0(name, ".link")
        bpt[,(sel.link):=0L]
        hx.bpt <- bpt[eval(sel), BPT.KEY, with=FALSE]
        if (nrow(hx.bpt) > 0) {
            bpt.edg <- data.table(
                bpt.5=paste(hx.bpt$chr.5, hx.bpt$pos.5, hx.bpt$str.5),
                bpt.3=paste(hx.bpt$chr.3, hx.bpt$pos.3, hx.bpt$str.3)
            )
            setkeyv(bpt.edg, c("bpt.5", "bpt.3"))
            bpt.edg <- unique(bpt.edg)
            g <- graph.data.frame(bpt.edg, directed=FALSE)
            bpt.x <- V(g)$name
            link <- as.integer(clusters(g)$membership)
            tmp <- data.table(bpt.x=bpt.x,
                              chn=link
                              )
            bpt.chn <- data.table(matrix(unlist(str_split(tmp$bpt.x, " ")), ncol=3, byrow=TRUE), tmp$chn)
            setnames(bpt.chn, c("chr", "pos", "str", "chn"))
            bpt.chn[,pos:=as.integer(pos)]
            setkeyv(bpt.chn, c("chr", "pos", "str"))
            setkeyv(hx.bpt, c("chr.5", "pos.5", "str.5"))
            hx.bpt$hx.link <- 0L
            hx.bpt[bpt.chn, hx.link:=chn]
        }
        setkeyv(bpt, BPT.KEY)
        setkeyv(hx.bpt, BPT.KEY)
        bpt[hx.bpt, (sel.link):=hx.link]
        bpt[get(sel.link)>0, l3:=name]
    }
    return(bpt)
}

.chainBreakpoints <- function(bpt, ann) {
    setkeyv(bpt, c("locus_id.5.1", "locus_id.3.1", "locus_id.5.2", "locus_id.3.2"))
    for (l3.sel in c("sl", "bs", "ts", "sv")) {
        sel.chain <- paste0(l3.sel, ".chain")
        bpt[,(sel.chain):=0L]
        chm <- bpt[l3==l3.sel,.N,by=CHM.KEY]
        loc.obs <- data.table(locus_id.1=c(as.character(chm$locus_id.5.1)),
                              locus_id.2=c(as.character(chm$locus_id.3.1)))
        loc.lnk <- data.table(locus_id.1=c(as.character(chm$locus_id.5.1),
                                           as.character(chm$locus_id.3.1),
                                           as.character(chm$locus_id.5.2)),
                              locus_id.2=c(as.character(chm$locus_id.3.2),
                                           as.character(chm$locus_id.5.2),
                                           as.character(chm$locus_id.3.2)))
        loc.edg <- rbind(loc.obs, loc.lnk)
        setkeyv(loc.edg, c("locus_id.1", "locus_id.2"))
        loc.edg <- unique(loc.edg)
        g <- graph.data.frame(loc.edg, directed=FALSE)
        locus_id <- V(g)$name
        chn <- as.integer(clusters(g)$membership)
        if (is.null(locus_id)) {
            locus_id <- character(0)
        }
        loc.chn <- data.table(locus_id=locus_id, chn=chn)
        setkey(loc.chn, "locus_id")
        ## 
        loc.obs$chain <- 0L
        setkeyv(loc.obs, "locus_id.1")
        loc.obs[loc.chn, chain:=chn]
        setkeyv(loc.obs, "locus_id.2")
        loc.obs[loc.chn, chain:=chn]
        ##
        setkeyv(loc.obs, c("locus_id.1", "locus_id.2"))
        bpt[loc.obs, (sel.chain):=chain]
    }
    return(bpt)
}

#' @export
collapseBreakpoints <- function(jnc, ann) {
    bpt <- .collapseBreakpoints(jnc)
    bpt <- .classifyBreakpoints(bpt, ann)
    bpt <- .motifBreakpoints(bpt)
    bpt <- .distanceBreakpoints(bpt, ann)
    bpt <- .defineBreakpointArtifacts(bpt, ann)
    bpt <- .associateBreakpoints(bpt, ann)
    bpt <- .topologyBreakpoints(bpt, ann)
    bpt <- .defineLevels(bpt, ann)
    bpt <- .defineHXBreakpoints(bpt, ann)
    bpt <- .linkBreakpoints(bpt, ann)
    bpt <- .chainBreakpoints(bpt, ann)
    bpt <- .naTo0(bpt, c("tx.dst", "gx.dst", "tx.pos.5", "tx.pos.3"))
    setkeyv(bpt, BPT.KEY)
    return(bpt)
}
