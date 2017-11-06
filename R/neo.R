ASM.EMPTY <- list(
    ctg = data.table(contig_id = character(), ctg.seq=list(), ctg.cov=numeric(), ctg.len=integer(), sv.chain=integer()),
    aln = data.table(contig_id = character(), transcript_id = character(), aln=list(), aln.len=integer(), sv.chain=integer()),
    chm = data.table(contig_id = character(), transcript_id = character(), x.rna.5 = list(), x.rna.3 = list(),
                     x.prot.5 = list(), x.prot.3 = list(), x.rna.len.5 = numeric(), x.rna.len.3 = numeric(),
                     x.prot.pos.5 = numeric(), x.prot.pos.3 = numeric(), sv.chain=integer()),
    fus = data.table(contig_id = character(), transcript_id.5 = character(), transcript_id.3 = character(),
                     x.seq.rna = list(), x.seq.prot = list(), inframe.5 = logical(), inframe.3 = logical(),
                     homology = logical(), splinter = logical(),
                     warn.5 = logical(), warn.3 = logical(), sv.chain=integer())
)


.filterTranscripts <- function(txs, ann) {
    txs.flt <- txs[order(tags)]
    txs.flt <- txs.flt[,.SD[1:min(ann$par$asm.tx.max, nrow(.SD))], by=locus_id]
    return(txs.flt)
}

.alignContigs <- function(ctgs, txs, ann) {
    alns <- data.table(expand.grid(contig_id=ctgs$contig_id, transcript_id=txs$transcript_id, stringsAsFactors=FALSE))
    setkey(alns, contig_id)
    setkey(ctgs, contig_id)
    alns <- alns[ctgs]
    setkey(alns, transcript_id)
    setkey(txs, transcript_id)
    alns <- alns[txs, nomatch=0]
    vpa <- Vectorize(pairwiseAlignment, vectorize.args=c("pattern", "subject"))
    alns[, aln:=list(vpa(mrna, ctg.seq, type="local", gapOpening=MAXINT/2, gapExtension=MAXINT/2))]
    alns[, aln.len:=as.integer(lapply(aln, nchar))]
    alns <- alns[(aln.len >= ann$par$asm.aln.len), .(contig_id, transcript_id, aln, aln.len)]
    setkey(alns, contig_id, transcript_id)
    return(alns)
}

.assembleChimeras <- function(ctgs, alns, txs, ann) {
    setkey(ctgs, contig_id)
    setkey(alns, contig_id)
    tmp <- ctgs[alns]
    setkey(tmp, transcript_id)
    setkey(txs, transcript_id)
    inps <- tmp[txs, nomatch=0, allow.cartesian=TRUE]
    ##
    ##                  cg5---------cg3
    ##     ____________ |    all    | ________________
    ##                 \_____________/ 
    ##  -----M----------------M----------M----*---------- (txr)
    ##       |          |           |         |
    ##       cds1       tx5---------tx3       cdsN
    ##
    outs <- lapply(seq_len(nrow(inps)), function(i) {
        with(inps[i], {
            out <- list(
                x.rna.5=list(DNAString()), x.rna.3=list(DNAString()),
                x.prot.5=list(AAString()), x.prot.3=list(AAString()),
                x.rna.pos.5=-1, x.rna.pos.3=-1,
                x.rna.len.5=-1, x.rna.len.3=-1,
                x.prot.pos.5=-1, x.prot.pos.3=-1
            )
            aln <- aln[[1]]
            mrna <- mrna[[1]]
            ctg.seq <- ctg.seq[[1]]
            txl <- nchar(mrna)
            tx5 <- start(pattern(aln))
            tx3 <- (tx5 - 1) + aln.len
            txr <- nchar(mrna) - tx3
            cgl <- nchar(ctg.seq)
            cg5 <- start(subject(aln))
            cg3 <- (cg5 - 1) + aln.len
            cgr <- nchar(ctg.seq) - cg3
            if (cgr > 0) { # mrna->ctg
                seq.5 <- subseq(mrna, end=tx5-1)
                seq.3 <- subseq(ctg.seq, start=cg5)
                out$x.rna.3 <- as.list(DNAStringSet(xscat(seq.5, seq.3)))
                out$x.rna.pos.3 <- as.integer(tx3+1) # first chimeric base
                out$x.rna.len.3 <- cgr
                if ((cdsN > 0) && ((cds1 + 2) <= tx3) && (cdsN > tx3)) {
                    out$x.prot.3 <- suppressWarnings(
                        as.list(AAStringSet(translate(subseq(out$x.rna.3[[1]], start=cds1)))))
                    out$x.prot.pos.3 <- as.integer(((tx3 - cds1 + 1) %/% 3) + 1)
                }
            }
            if (cg5 > 1) { # ctg->mrna
                seq.5 <- subseq(ctg.seq, end=cg3)
                seq.3 <- subseq(mrna, start=tx3+1)
                out$x.rna.5 <- as.list(DNAStringSet(xscat(seq.5, seq.3)))
                out$x.rna.pos.5 <- as.integer(cg5-1) # last chimeric base
                out$x.rna.len.5 <- cg5-1
                if ((cdsN > 0) && (cds1 >= tx5)) {
                    cds1.x <- cds1 - tx5 + cg5
                    out$x.prot.5 <- suppressWarnings(
                        as.list(AAStringSet(translate(subseq(out$x.rna.5[[1]], start=cds1.x)))))
                    out$x.prot.pos.5 <- 0L
                }
            }
            out <- out[colnames(ASM.EMPTY$chm)[-c(1,2,11)]]
            return(out)
        })
    })
    outs0 <- ASM.EMPTY$chm[,-c(1,2,11)]
    outs <- cbind(inps[,.(contig_id, transcript_id)], rbind(outs0, rbindlist(outs)))
    chms <- outs[((x.rna.len.5>-1) | (x.rna.len.3>-1))]
    setkey(chms, contig_id, transcript_id)
    return(chms)
}

.assembleFusions <- function(ctgs, alns, chms, txs, ann) {
    setkey(chms, contig_id)
    setkey(ctgs, contig_id)
    tbls <- ctgs[chms]
    setkey(txs, transcript_id)
    setkey(tbls, transcript_id)
    tbls <- txs[tbls]
    setkey(tbls, contig_id, transcript_id)
    setkey(alns, contig_id, transcript_id)
    tbls <- alns[tbls]
    a5 <- tbls[x.rna.len.3 > -1] # mrna->ctg
    a3 <- tbls[x.rna.len.5 > -1] # ctg->mrna
    setkey(a5, contig_id)
    setkey(a3, contig_id)
    inps <- merge(
        a5[,.(contig_id, ctg.seq, ctg.len,
              locus_id.5=locus_id, transcript_id.5=transcript_id, aln.5=aln, aln.len.5=aln.len,
              mrna.5=mrna, prot.5=prot, cds1.5=cds1, cdsN.5=cdsN)],
        a3[,.(contig_id,
              locus_id.3=locus_id, transcript_id.3=transcript_id, aln.3=aln, aln.len.3=aln.len,
              mrna.3=mrna, prot.3=prot, cds1.3=cds1, cdsN.3=cdsN)], allow.cartesian=TRUE)
    inps <- inps[locus_id.5!=locus_id.3]
    ##               cds1.3       tx5.3    tx3.3     cdsN.3
    ##               |            |            |     |
    ##  -------------M------------------M------------*------- mrna.3
    ##                           /¯¯¯¯¯¯¯¯¯¯¯¯¯¯\    ?
    ##                          / |            | ¯¯¯¯cgr.3
    ##                         ?  cg5.3    cg3.3 
    ##         cg5.5   cg3.5  ?
    ##    ?___ |           | /
    ##        \_____________/ 
    ##  -M-----------------M--------------------*------------ mrna.5
    ##   |     |           |                    |
    ##  cds1.5 tx5.5   tx3.5                    cdsN.5
    outs <- lapply(seq_len(nrow(inps)), function(i) {
        with(inps[i], {
            aln.5 <- aln.5[[1]]
            aln.3 <- aln.3[[1]]
            mrna.5 <- mrna.5[[1]]
            mrna.3 <- mrna.3[[1]]
            prot.5 <- prot.5[[1]]
            prot.3 <- prot.3[[1]]
            ctg.seq <- ctg.seq[[1]]
            txl.5 <- nchar(mrna.5)
            txl.3 <- nchar(mrna.3)
            tx5.5 <- start(pattern(aln.5))
            tx5.3 <- start(pattern(aln.3))
            tx3.5 <- (tx5.5 - 1) + aln.len.5
            tx3.3 <- (tx5.3 - 1) + aln.len.3
            txr.5 <- nchar(mrna.5) - tx3.5
            txr.3 <- nchar(mrna.3) - tx3.3
            cg5.5 <- start(subject(aln.5))
            cg5.3 <- start(subject(aln.3))
            cg3.5 <- (cg5.5 - 1) + aln.len.5
            cg3.3 <- (cg5.3 - 1) + aln.len.3
            cgr.3 <- ctg.len - cg3.3
            x.seq.rna <- DNAString()
            x.seq.cds <- DNAString()
            if ((cg5.5 < cg5.3) && (cg3.5 < cg3.3)) {
                seq.5 <- subseq(mrna.5, end=tx5.5-1)
                seq.x <- subseq(ctg.seq, start=cg5.5, end=cg3.3)
                seq.3 <- subseq(mrna.3, start=tx3.3+1)
                x.seq.rna <- DNAString(xscat(seq.5, seq.x, seq.3))
                if ((cdsN.5 > 0) && ((cds1.5+2) <= tx3.5)) { ## translate from 5' gene
                    x.seq.cds <- subseq(x.seq.rna, start=cds1.5)
                } else if ((cdsN.3 > 0) && (cds1.3 >= tx5.3)) { ## translate from 3' gene
                    cds1.3.rel <- (cds1.3-tx5.3) + ((tx5.5-1)+(cg5.3-cg5.5+1))
                    x.seq.cds <- subseq(x.seq.rna, start=cds1.3.rel)
                }
            }
            x.seq.prot <- AAString(str_match(suppressWarnings(translate(x.seq.cds)), "[^\\*]*")[,1]) # truncate till *
            aln.prot.5 <- pairwiseAlignment(x.seq.prot, prot.5, gapOpening=MAXINT/2, gapExtension=MAXINT/2, type="local")
            aln.prot.3 <- pairwiseAlignment(x.seq.prot, prot.3, gapOpening=MAXINT/2, gapExtension=MAXINT/2, type="local")
            aln.mrna.53 <- pairwiseAlignment(mrna.5, mrna.3, type="local")
            aln.prot.53 <- pairwiseAlignment(prot.5, prot.3, type="local")
            list(
                x.seq.rna=as.list(DNAStringSet(x.seq.rna)),
                x.seq.prot=as.list(AAStringSet(x.seq.prot)),
                inframe.5=nchar(aln.prot.5)>ann$par$asm.hom.len %/% 3,
                inframe.3=nchar(aln.prot.3)>ann$par$asm.hom.len %/% 3,
                homology=(nchar(aln.prot.53) > ann$par$asm.hom.len %/% 3) ||
                         (nchar(aln.mrna.53) > ann$par$asm.hom.len),
                splinter=(cg5.3-cg3.5-1)>0,
                warn.5=cg5.5>1,
                warn.3=cgr.3>0
            )
        })
    })
    outs0 <- ASM.EMPTY$fus[,-c(1,2,3,12)]
    outs <- cbind(inps[,.(contig_id, transcript_id.5, transcript_id.3)], rbind(outs0, rbindlist(outs)))
    return(outs)
}

.filterTS <- function(bun, ann) {
    bpt <- copy(bun$bpt)
    bpt[,max.rec.5:=suppressWarnings(max(rec.5)),by=eval(CHM.KEY[c(2,4)])]
    bpt[,max.rec.3:=suppressWarnings(max(rec.3)),by=eval(CHM.KEY[c(3,5)])]
    filterBundle(bpt[max.rec.5<=ann$par$asm.rec.max & max.rec.3<=ann$par$asm.rec.max], bun, ann)
}

.assembleChain <- function(bun, ann) {
    ctgs <- .runInchworm(bun, ann)
    if (nrow(ctgs)>0) {
        maps <- .runGmap(ctgs, ann)
        setkey(maps, contig_id)
        setkey(ctgs, contig_id)
        ctgs[maps, ":="(gmap.chr.5=chr.5, gmap.pos.5=pos.5, gmap.str.5=str.5,
                                gmap.chr.3=chr.3, gmap.pos.3=pos.3, gmap.str.3=str.3,
                                gmap.nm=NM, gmap.mapq=mapq, gmap.status=XO)]
        ctgs <- cbind(ctgs)
    } else {
        ctgs <- cbind(ctgs, GMAP.EMPTY)
    }
    ctgs[, ":="(gmap.chr.5=factor(gmap.chr.5, levels=seqlevels(ann$seqi), ordered=TRUE),
                    gmap.chr.3=factor(gmap.chr.3, levels=seqlevels(ann$seqi), ordered=TRUE),
                    gmap.str.5=factor(gmap.str.5, levels=c("+", "-"), ordered=TRUE),
                    gmap.str.3=factor(gmap.str.3, levels=c("+", "-"), ordered=TRUE)
                    )]
    return(ctgs)
}

.stitchChain <- function(bun, ann) {
    ## prepare
    lgt <- .lgt(ann)
    txs <- .prepareTxs(bun, lgt)
    txs <- .filterTranscripts(txs, ann)
    ## assemble & stitch
    ctgs <- .assembleChain(bun, ann)
    alns <- .alignContigs(ctgs, txs, ann)
    chms <- .assembleChimeras(ctgs, alns, txs, ann)
    fuss <- .assembleFusions(ctgs, alns, chms, txs, ann)
    stitch <- list(ctg=ctgs, aln=alns, chm=chms, fus=fuss)
    return(stitch)
}


.stitchChains <- function(bun, ann) {
    chain.type <- paste0(bun$type, ".chain")
    bun <- .filterTS(bun, ann)
    buns <- splitBundle(bun, chain.type)
    stitches <- mclapply(names(buns), function(chain) {
        sel.bun <- buns[[chain]]
        stitch <- .stitchChain(sel.bun, ann)
        stitch <- lapply(stitch, function(tbl) tbl[,(chain.type):=as.integer(chain)])
        return(stitch)
    })
    if (length(stitches)>0) {
        nasm <- c("ctg", "aln", "chm", "fus")
        asm <- lapply(nasm, function(x) rbindlist(lapply(stitches, "[[", x)))
        names(asm) <- nasm
    } else {
        asm <- ASM.EMPTY
    }
    return(asm)
}

.extractChimeric <- function(sv.asm, ann) {
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
