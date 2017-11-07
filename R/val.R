.breakpointValidateMinimap2 <- function(aln.mm2, ann) {
    tmp <- aln.mm2[order(contig_id, first_base_q)][n.seg==2]
    if (nrow(tmp)>0) {
        tmp <- cbind(
            tmp[seq(1, nrow(tmp), 2), ..BPT.KEY],
            tmp[seq(1, nrow(tmp), 2), .(contig_id, mm.flag=flag)],
            tmp[seq(1, nrow(tmp), 2), .(mm.chr.5=rname, mm.pos.5=ifelse(strand=="+", last_base_r, first_base_r),
                                        mm.str.5=strand, mm.mapq.5=mapq)],
            tmp[seq(2, nrow(tmp), 2), .(mm.chr.3=rname, mm.pos.3=ifelse(strand=="+", first_base_r, last_base_r),
                                        mm.str.3=strand, mm.mapq.3=mapq)]
        )
        tmp[,mm2.val:=FALSE]
        tmp[chr.5==mm.chr.5 & chr.3==mm.chr.3 &
            str.5==mm.str.5 & str.3==mm.str.3 &
            abs(pos.5-mm.pos.5) <= ann$par$asm.mm2.maxdst &
            abs(pos.3-mm.pos.3) <= ann$par$asm.mm2.maxdst &
            mm.mapq.5 >= ann$par$asm.mm2.minmapq &
            mm.mapq.3 >= ann$par$asm.mm2.minmapq, mm2.val:=TRUE]
        bpt.val.mm2 <- tmp[,.(mm2.val=any(mm2.val)),by=BPT.KEY]
    } else {
        bpt.val.mm2 <- tmp[,..BPT.KEY]
        bpt.val.mm2[,mm2.val:=logical()]
    }
    setkeyv(bpt.val.mm2, BPT.KEY)
    return(bpt.val.mm2)
}

.breakpointValidateGmap <- function(aln.gmap, ann) {
    bpt.val.gmap <- aln.gmap[,.(gmap.val=any(XO=="UT")), by=BPT.KEY]
    setkeyv(bpt.val.gmap, BPT.KEY)
    return(bpt.val.gmap)
}

#' @export
validateBreakpoints <- function(bun, ann) {
    val.gmap <- .breakpointValidateGmap(bun$aln.gmap, ann)
    val.mm2 <- .breakpointValidateMinimap2(bun$aln.mm2, ann)
    bun$bpt[,mm2.valid:=FALSE]
    bun$bpt[val.mm2, mm2.valid:=mm2.val]
    bun$bpt[,gmap.valid:=FALSE]
    bun$bpt[val.gmap, gmap.valid:=gmap.val]
    return(bun)    
}
