.defineBreakpointArtifacts <- function(xpt, ann) {
    par <- ann$par
    bpt.5 <- with(xpt, GRanges(chr.5, IRanges(ifelse(str.5=="+", pos.5-par$jnc.ovr, pos.5-ceiling(par$jnc.ovr/3)),
                                              ifelse(str.5=="+", pos.5+ceiling(par$jnc.ovr/3), pos.5+par$jnc.ovr)), str.5))
    bpt.3 <- with(xpt, GRanges(chr.3, IRanges(ifelse(str.3=="+", pos.3-ceiling(par$jnc.ovr/3), pos.3-par$jnc.ovr),
                                              ifelse(str.3=="+", pos.3+par$jnc.ovr, pos.3+ceiling(par$jnc.ovr/3))), str.3))

    ## retros
    ret.5 <- suppressWarnings(bpt.5 %over% ann$retros)
    ret.3 <- suppressWarnings(bpt.3 %over% ann$retros)

    ## SVs
    svs.5 <- suppressWarnings(bpt.5 %over% ann$svs[width(ann$svs) > 20000])
    svs.3 <- suppressWarnings(bpt.3 %over% ann$svs[width(ann$svs) > 20000])

    ## LOW
    low.5 <- suppressWarnings(bpt.5 %over% ann$low[width(ann$low) >= par$jnc.ovr])
    low.3 <- suppressWarnings(bpt.3 %over% ann$low[width(ann$low) >= par$jnc.ovr])

    ## XLW
    xlw.5 <- xpt$avg.low.5 > par$jnc.low*4
    xlw.3 <- xpt$avg.low.3 > par$jnc.low*4
    
    ## repeats
    ## implementation is a workaround around new findOverlaps implementation NCList
    tmp <- ann$repeats[width(ann$repeats) < 2000]
    rep <- split(tmp, factor(ceiling(seq_along(tmp) / 100000)))
    rep.5 <- rowSums(as.data.table(lapply(rep, function(r) {suppressWarnings(bpt.5 %over% r)}))) > 0
    rep.3 <- rowSums(as.data.table(lapply(rep, function(r) {suppressWarnings(bpt.3 %over% r)}))) > 0
    rm(tmp, rep); gc()
    
    ## medium repeats
    mtr <- ann$repeats[grep("\\([AGCT]{3,}\\)n", ann$repeats$name)]
    mtr <- mtr[(mtr$score / width(mtr) > 0.90)]
    mtr.5 <- suppressWarnings(bpt.5 %over% mtr)
    mtr.3 <- suppressWarnings(bpt.3 %over% mtr)

    ## simple repeats
    sim <- ann$repeats[grep("\\([AGCT]{1,2}\\)n", ann$repeats$name)]
    sim.5 <- suppressWarnings(bpt.5 %over% sim)
    sim.3 <- suppressWarnings(bpt.3 %over% sim)
    
    ## Alu repeats
    alu <- ann$repeats[grep("Alu", ann$repeats$name)]
    alu.5 <- suppressWarnings(bpt.5 %over% alu)
    alu.3 <- suppressWarnings(bpt.3 %over% alu)

    ## satellite repeats
    msr <- ann$repeats[grep("MSR1", ann$repeats$name)]
    msr.5 <- suppressWarnings(bpt.5 %over% msr)
    msr.3 <- suppressWarnings(bpt.3 %over% msr)

    ## L1 repeats
    l1 <- ann$repeats[grep("^L1|HAL1", ann$repeats$name)]
    l1.5 <- suppressWarnings(bpt.5 %over% l1)
    l1.3 <- suppressWarnings(bpt.3 %over% l1)
    
    ## L2 repeats
    l2 <- ann$repeats[grep("^L2", ann$repeats$name)]
    l2.5 <- suppressWarnings(bpt.5 %over% l2)
    l2.3 <- suppressWarnings(bpt.3 %over% l2)

    ## L3 repeats
    l3 <- ann$repeats[grep("^L3", ann$repeats$name)]
    l3.5 <- suppressWarnings(bpt.5 %over% l3)
    l3.3 <- suppressWarnings(bpt.3 %over% l3)

    ## MLT repeats
    mlt <- ann$repeats[grep("^MLT", ann$repeats$name)]
    mlt.5 <- suppressWarnings(bpt.5 %over% mlt)
    mlt.3 <- suppressWarnings(bpt.3 %over% mlt)

    ## MER repeats
    mer <- ann$repeats[grep("^MER", ann$repeats$name)]
    mer.5 <- suppressWarnings(bpt.5 %over% mer)
    mer.3 <- suppressWarnings(bpt.3 %over% mer)

    ## ERV repeats
    erv <- ann$repeats[grep("ERV", ann$repeats$name)]
    erv.5 <- suppressWarnings(bpt.5 %over% erv)
    erv.3 <- suppressWarnings(bpt.3 %over% erv)

    ## LTR repeats
    ltr <- ann$repeats[grep("LTR", ann$repeats$name)]
    ltr.5 <- suppressWarnings(bpt.5 %over% ltr)
    ltr.3 <- suppressWarnings(bpt.3 %over% ltr)
    
    ## ALTLOCs 
    alt.5 <- suppressWarnings(bpt.5 %over% ann$altlocs)
    alt.3 <- suppressWarnings(bpt.3 %over% ann$altlocs)
    
    ## segmental duplications
    seg.5 <- with(ann$segdups, GRanges(chr.5, IRanges(start.5, end.5), strand="*"))
    seg.3 <- with(ann$segdups, GRanges(chr.3, IRanges(start.3, end.3), strand="*"))
    b2s.5 <- data.table(as.data.frame(findOverlaps(bpt.5, seg.5)))
    setkeyv(b2s.5, c("queryHits", "subjectHits"))
    b2s.3 <- data.table(as.data.frame(findOverlaps(bpt.3, seg.3)))
    setkeyv(b2s.3, c("queryHits", "subjectHits"))
    seg.53 <- unique(merge(b2s.5, b2s.3, allow.cartesian=TRUE)$queryHits)

    ## known artifacts
    aln.5 <- with(ann$artifacts, GRanges(chr.5, IRanges(start.5, end.5), strand="*"))
    aln.3 <- with(ann$artifacts, GRanges(chr.3, IRanges(start.3, end.3), strand="*"))
    b2a.5 <- data.table(as.data.frame(findOverlaps(bpt.5, aln.5)))
    setkeyv(b2a.5, c("queryHits", "subjectHits"))
    b2a.3 <- data.table(as.data.frame(findOverlaps(bpt.3, aln.3)))
    setkeyv(b2a.3, c("queryHits", "subjectHits"))
    aln.53 <- unique(merge(b2a.5, b2a.3, allow.cartesian=TRUE)$queryHits)
    
    ## rn7s
    rn7s <- ann$genes[grep("^RN7S", ann$genes$gene_name)]
    strand(rn7s) <- "*"
    rn7s.5 <- suppressWarnings(bpt.5 %over% rn7s)
    rn7s.3 <- suppressWarnings(bpt.3 %over% rn7s)
    
    ## RP[SL]
    rpls <- ann$genes[grep("^M{0,1}RP[LS][0-9AP]", ann$genes$gene_name)]
    strand(rpls) <- "*"
    rpls.5 <- suppressWarnings(bpt.5 %over% rpls)
    rpls.3 <- suppressWarnings(bpt.3 %over% rpls)
    
    ## SNOR[A-Z]
    sno <- ann$snornas
    strand(sno) <- "*"
    sno.5 <- suppressWarnings(bpt.5 %over% sno)
    sno.3 <- suppressWarnings(bpt.3 %over% sno)

    ## PSEUDO
    psi <- ann$genes[grep("pseudo", ann$genes$gene_type, perl=TRUE)]
    strand(psi) <- "*"
    psi.5 <- suppressWarnings(bpt.5 %over% psi)
    psi.3 <- suppressWarnings(bpt.3 %over% psi)
    
    ## RNU
    ug <- ann$genes[grep("^RNU", ann$genes$gene_name), NULL]
    ur <- ann$repeats[grep("^U[0-9]", ann$repeats$name), NULL]
    rnu <- reduce(c(ug, ur))
    rnu.5 <- suppressWarnings(bpt.5 %over% rnu)
    rnu.3 <- suppressWarnings(bpt.3 %over% rnu)
    
    ## IGX
    igx.5 <- bpt.5 %over% ann$igx
    igx.3 <- bpt.3 %over% ann$igx

    ## GTB
    gtb.5 <- bpt.5 %over% ann$gtb
    gtb.3 <- bpt.3 %over% ann$gtb

    ## ChrM
    cmt <- ann$loci[seqnames(ann$loci)=="chrM"]
    gmt <- ann$svs[ann$svs$meth=="NUMT_umich"]
    cmt.5 <- suppressWarnings(bpt.5 %over% cmt | bpt.5 %over% gmt)
    cmt.3 <- suppressWarnings(bpt.3 %over% cmt | bpt.3 %over% gmt)
    
    ##
    xpt$art.5 <- "nil"
    xpt$art.3 <- "nil"
    xpt[rep.5, art.5:="rep"]
    xpt[rep.3, art.3:="rep"]
    xpt[low.5, art.5:="low"]
    xpt[low.3, art.3:="low"]
    ## not so bad
    xpt[svs.5, art.5:="sv"]
    xpt[svs.3, art.3:="sv"]
    xpt[l1.5, art.5:="l1"]
    xpt[l1.3, art.3:="l1"]
    xpt[l2.5, art.5:="l2"]
    xpt[l2.3, art.3:="l2"]
    xpt[l3.5, art.5:="l3"]
    xpt[l3.3, art.3:="l3"]
    xpt[mlt.5, art.5:="mlt"]
    xpt[mlt.3, art.3:="mlt"]
    xpt[mer.5, art.5:="mer"]
    xpt[mer.3, art.3:="mer"]
    xpt[erv.5, art.5:="erv"]
    xpt[erv.3, art.3:="erv"]
    xpt[ltr.5, art.5:="ltr"]
    xpt[ltr.3, art.3:="ltr"]
    xpt[igx.5, art.5:="igx"]
    xpt[igx.3, art.3:="igx"]
    ## expression
    xpt[rn7s.5, art.5:="rn7"]
    xpt[rn7s.3, art.3:="rn7"]
    xpt[rpls.5, art.5:="rpx"]
    xpt[rpls.3, art.3:="rpx"]
    xpt[sno.5, art.5:="sno"]
    xpt[sno.3, art.3:="sno"]
    xpt[rnu.5, art.5:="rnu"]
    xpt[rnu.3, art.3:="rnu"]
    ## "structural" repeats
    xpt[alu.5, art.5:="alu"]
    xpt[alu.3, art.3:="alu"]
    xpt[psi.5, art.5:="psi"]
    xpt[psi.3, art.3:="psi"]
    xpt[mtr.5, art.5:="mtr"]
    xpt[mtr.3, art.3:="mtr"]
    xpt[alt.5, art.5:="alt"]
    xpt[alt.3, art.3:="alt"]
    ## segmental duplications
    xpt[seg.53, ":="(art.5="seg", art.3="seg")]
    ## alignment artifacts
    xpt[aln.53, ":="(art.5="aln", art.3="aln")]
    ## bad repeats
    xpt[xlw.5, art.5:="xlw"]
    xpt[xlw.3, art.3:="xlw"]
    xpt[ret.5, art.5:="ret"]
    xpt[ret.3, art.3:="ret"]
    xpt[msr.5, art.5:="msr"]
    xpt[msr.3, art.3:="msr"]
    xpt[sim.5, art.5:="sim"]
    xpt[sim.3, art.3:="sim"]
    xpt[gtb.5, art.5:="gtb"]
    xpt[gtb.3, art.3:="gtb"]
    ## chromosome M
    xpt[cmt.5, art.5:="cmt"]
    xpt[cmt.3, art.3:="cmt"]
    ## levels
    art.levels <- sort(unique(c(xpt$art.5, xpt$art.3)))
    xpt$art.5 <- factor(xpt$art.5, levels=art.levels)
    xpt$art.3 <- factor(xpt$art.3, levels=art.levels)
    setkeyv(xpt, BPT.KEY)
    ####
    xpt$art <- TRUE
    xpt[(ss.5 == "D" & ss.3 == "A"), art:=FALSE] # known transcript splice
    xpt[
        ## cannot be same "long" repeat
        (!(
           art.5 %in% c("l1", "l2", "l3", "mlt", "mer", "erv", "ltr", "igx", "sv")
           &
           art.3 %in% c("l1", "l2", "l3", "mlt", "mer", "erv", "ltr", "igx", "sv")
           &
           art.5 == art.3
         )
        )
        ## cannot be expression artifact on both ends
        &
        (!(
           art.5 %in% c("rn7", "rpx", "sno", "rnu")
           &
           art.3 %in% c("rn7", "rpx", "sno", "rnu")
         )
        )
        ## cannot be bad aritfact 
        &
        (!(
           (art.5 %in% c("alu", "psi", "mtr", "alt")) & (art.3 != "nil")
           |
           (art.3 %in% c("alu", "psi", "mtr", "alt")) & (art.5 != "nil")
          )
        )
        ## cannot be ugly artifact on either end
        &
        (!(
           art.5 %in% c("ret", "cmt", "msr", "sim", "seg", "xlw", "gtb", "aln")
           |
           art.3 %in% c("ret", "cmt", "msr", "sim", "seg", "xlw", "gtb", "aln")
         )
        )
        ,
        art:=FALSE
        ]
    return(xpt)
}
