.extractChimeric <- function(sv.sew, ann) {
    ctg.hq <- sv.asm$ctg[ctg.cov > ann$par$asm.ctg.cov*3 & ctg.len > ann$par$asm.ctg.len*2]
    setkey(ctg.hq, sv.chain, contig_id)
    chm.hq <- sv.asm$chm[(x.rna.len.3 > -1 & x.rna.len.5 == -1) & x.prot.pos.3 > 1]
    setkey(chm.hq, sv.chain, contig_id)
    mrg.hq <- ctg.hq[chm.hq, nomatch=0]
    setkey(mrg.hq, transcript_id, contig_id, sv.chain)
    fus.hq <- sv.asm$fus[(!homology)]
    setkey(fus.hq, transcript_id.5, contig_id, sv.chain)
    neo.hq <- fus.hq[mrg.hq]
    setkey(neo.hq, transcript_id.5)
    neo.hq[lgt, ":="(gene_id=gene_id, tags=tags, cds1=cds1)]
    neo.hq.pri <- neo.hq[order(-inframe.3, -x.rna.len.3, tags)]
    neo.hq.pri.sel <- neo.hq.pri[order(-inframe.3, -x.rna.len.3, tags),.SD[1], by=.(sv.chain, contig_id, gene_id)]
    setkey(neo.hq.pri.sel, sv.chain, gene_id, contig_id)
    return(neo.hq.pri.sel)
}
