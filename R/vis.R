#' @export
svSig <- function(sv.bun) {
    if (nrow(sv.bun$bpt)>0) {
        ## 
        all.topo <- factor(levels(sv.bun$bpt$topo), levels(sv.bun$bpt$topo), ordered=TRUE)
        all.dst <- factor(levels(sv.bun$bpt$dst), levels(sv.bun$bpt$dst), ordered=TRUE)
        sv.grid <- data.table(expand.grid(topo=all.topo, dst=all.dst))
        sv.grid <- sv.grid[!(dst %in% c("frg", "loc"))]
        ##
        chm.topo.dst <- sv.bun$bpt[, .(topo=topo[1], dst=dst[1]), CHM.KEY]
        sv.sig <- chm.topo.dst[, .(event.n=.N), by=.(topo, dst)]
        sv.sig[,event.frq:=event.n/sum(event.n)]
        ##
        setkey(sv.sig, topo, dst)
        setkey(sv.grid, topo, dst)
        sv.sig <- .naTo0(sv.sig[sv.grid])
    } else {
        sv.sig <- NULL
    }
    return(sv.sig)
}

#' @export
svSigPlot <- function(sv.sig) {
    plt <- ggplot2::ggplot(sv.sig) +
           ggplot2::aes(x=dst, y=topo, fill=event.frq) +
           ggplot2::geom_tile() +
           ggplot2::scale_fill_gradient(low="white", high="red", name="freq.", labels=scales::percent) +
           ggplot2::xlab("distance") +
           ggplot2::ylab("topology") +
           ggplot2::coord_fixed() +
           ggplot2::theme_bw(base_size=16)
    return(plt)
}

#' @export
svChainCircos <- function(sv.rep, ann) {
    circlize::circos.clear()
    cyt <- system.file(package="codac", "extdata", "hg38.cyt")
    chr <- seqnames(ann$seqi)
    circlize::circos.initializeWithIdeogram(cyt, chromosome.index=chr)
    ## select unique breakpoints
    sv.rep.unq <- sv.rep[!duplicated(sv.rep[,eval(CHM.KEY),with=FALSE])]
    sv.rep.unq[, sv.chain.n:=.N, by=sv.chain]
    sv.rep.unq[sv.chain.n==1, sv.chain:=0L]
    sv.rep.unq[,":="(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), NA_character_)),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), NA_character_))
    )]
    ## draw links
    sv.by.chain <- split(sv.rep.unq, sv.rep.unq$sv.chain)
    sv.by.chain <- sv.by.chain[sort(names(sv.by.chain))]
    for (idx in seq_along(sv.by.chain)) {
        sv.rep.sel <- sv.by.chain[[idx]]
        bed.sel.5 <- sv.rep.sel[,.(chr=chr.5, start=pos.5, end=pos.5)]
        bed.sel.3 <- sv.rep.sel[,.(chr=chr.3, start=pos.3, end=pos.3)]
        circlize::circos.genomicLink(bed.sel.5, bed.sel.3, lwd=4, col=c("darkgray", COLS)[idx])
    }
    ## draw labels
    tbed <- rbind(sv.rep.unq[,.(chr=chr.5, start=pos.5, end=pos.5, label=gene.5)],
                  sv.rep.unq[,.(chr=chr.3, start=pos.3, end=pos.3, label=gene.3)])
    tbed <- tbed[!is.na(label) & !duplicated(label)][order(chr, start)]
    circlize::circos.genomicTrackPlotRegion(tbed, ylim = c(0, 10), track.height = 0.05, bg.border = NA,
        panel.fun = function(region, value, ...) {
            circlize::circos.genomicText(region, value, y = 0, labels.column = 1,
                               facing = "clockwise", adj = c(1, 0.5),
                               posTransform = circlize::posTransform.text, cex = 0.75, padding=1)
        })
    current.track.index = circlize::get.cell.meta.data("track.index")
    circlize::circos.genomicPosTransformLines(tbed, track.index = current.track.index, direction = "inside",
        posTransform = function(region, value) {
            circlize::posTransform.text(region, y = 0, labels = value[[1]],
                                        track.index = current.track.index,
                                        cex = 0.75, 
                                        padding = 1
                                        )
        })
}

statPlot <- function(stat, ann) {
    var <- c("eff.frg", "aln.rte", "spl.rte", "sum.nsj", "lng.rte", "gene.rte",
                 "art.rte", "dup.rte", "lig.rte", "bs.rte", "ts.lfc", "sl.rte")
    writeLines(var, "blah")


    
    stat$aln.rte
}
