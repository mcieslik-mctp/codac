#' @export
neoReport <- function(asm, ann) {
    ## prepare annotation
    lgt <- .lgt(ann)
    setkey(lgt, transcript_id)
    ## merge
    ctg.map <- asm$ctg[!is.na(locus_id.5.1) & !is.na(locus_id.3.1)]
    setkey(ctg.map, sv.chain, contig_id)
    chm.3 <- asm$chm[(x.rna.len.3 > -1 & x.rna.len.5 == -1) & x.prot.pos.3 > 1]
    setkey(chm.3, sv.chain, contig_id)
    chm.map.3 <- ctg.map[chm.3, nomatch=0]
    setkey(chm.map.3, transcript_id, contig_id, sv.chain)
    setkey(asm$fus, transcript_id.5, contig_id, sv.chain)
    chm.map.3 <- asm$fus[chm.map.3, allow.cartesian=TRUE]
    setkey(chm.map.3, transcript_id.5)
    chm.map.3[lgt, ":="(gene_id.5=gene_id, gene_name.5=gene_name, tags=tags, cds1=cds1)]
    setkey(chm.map.3, transcript_id.3)
    chm.map.3[lgt, ":="(gene_id.3=gene_id, gene_name.3=gene_name)]
    ## select 1 per contig/chain/gene
    chm.map.3 <- chm.map.3[order(-inframe.3, -x.rna.len.3, tags)]
    neo <- chm.map.3[,.SD[1], by=.(sv.chain, contig_id)]
    setkey(neo, sv.chain, contig_id)
    return(neo)
}
