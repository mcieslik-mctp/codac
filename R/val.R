.breakpointValidate <- function(aln, maxdst, minmapq) {
    tmp <- aln[order(contig_id, first_base_q)][n.seg==2]
    if (nrow(tmp)>0) {
        tmp <- cbind(
            tmp[seq(1, nrow(tmp), 2), ..BPT.KEY],
            tmp[seq(1, nrow(tmp), 2), .(contig_id, x.flag=flag)],
            tmp[seq(1, nrow(tmp), 2), .(x.chr.5=rname, x.pos.5=ifelse(strand=="+", last_base_r, first_base_r),
                                        x.str.5=strand, x.mapq.5=mapq)],
            tmp[seq(2, nrow(tmp), 2), .(x.chr.3=rname, x.pos.3=ifelse(strand=="+", first_base_r, last_base_r),
                                        x.str.3=strand, x.mapq.3=mapq)]
        )
        tmp[,x.val:=FALSE]
        tmp[chr.5==x.chr.5 & chr.3==x.chr.3 &
            str.5==x.str.5 & str.3==x.str.3 &
            abs(pos.5-x.pos.5) <= maxdst &
            abs(pos.3-x.pos.3) <= maxdst &
            x.mapq.5 >= minmapq &
            x.mapq.3 >= minmapq, x.val:=TRUE]
        bpt.val.x <- tmp[,.(x.val=any(x.val)),by=BPT.KEY]
    } else {
        bpt.val.x <- tmp[,..BPT.KEY]
        bpt.val.x[,x.val:=logical()]
    }
    return(bpt.val.x)
}

#' @export
validateBreakpoints <- function(bun, ann) {
    val.mm2 <- .breakpointValidate(bun$aln.mm2, ann$par$asm.mm2.maxdst, ann$par$asm.mm2.minmapq)
    setnames(val.mm2, "x.val", "mm2.val")
    setkeyv(val.mm2, BPT.KEY)
    val.gmap <- .breakpointValidate(bun$aln.gmap, ann$par$asm.gmap.maxdst, ann$par$asm.gmap.minmapq)
    setnames(val.gmap, "x.val", "gmap.val")
    setkeyv(val.gmap, BPT.KEY)
    bun$bpt[,mm2.valid:=FALSE]
    bun$bpt[val.mm2, mm2.valid:=mm2.val]
    bun$bpt[,gmap.valid:=FALSE]
    bun$bpt[val.gmap, gmap.valid:=gmap.val]
    return(bun)    
}
