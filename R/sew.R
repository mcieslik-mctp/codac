ASM.EMPTY <- list(
    ctg = data.table(contig_id = character(), locus_id.5.1=character(), locus_id.3.1=character(),
                     sv.chain=integer(), ctg.seq=list(), ctg.cov=numeric(), ctg.len=integer()),
    aln = data.table(contig_id = character(), transcript_id = character(), aln=list(), aln.len=integer(), sv.chain=integer()),
    chm = data.table(contig_id = character(), transcript_id = character(), x.rna.5 = list(), x.rna.3 = list(),
                     x.prot.5 = list(), x.prot.3 = list(), x.rna.len.5 = numeric(), x.rna.len.3 = numeric(),
                     x.prot.pos.5 = numeric(), x.prot.pos.3 = numeric(), sv.chain=integer()),
    fus = data.table(contig_id = character(), transcript_id.5 = character(), transcript_id.3 = character(),
                     x.seq.rna = list(), x.seq.prot = list(), inframe.5 = logical(), inframe.3 = logical(),
                     homology = logical(), splinter = logical(),
                     warn.5 = logical(), warn.3 = logical(), sv.chain=integer())
)

.alignContigsToTranscripts <- function(ctgs, txs, ann) {
    alns <- data.table(expand.grid(contig_id=ctgs$contig_id, transcript_id=txs$transcript_id, stringsAsFactors=FALSE))
    setkey(alns, contig_id)
    setkey(ctgs, contig_id)
    alns <- alns[ctgs,allow.cartesian=TRUE]
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

.alignWrapperChain <- function(f, ctg, ann) {
    bam <- f(ctg, ann)
    setkey(ctg, contig_id)
    setkey(bam, contig_id)
    aln <- bam[ctg]
    CTG.KEY <- colnames(CTG.EMPTY)
    BAM.KEY <- setdiff(intersect(colnames(BAM.EMPTY.MINIMAP2), colnames(BAM.EMPTY.GMAP)), CTG.KEY)
    UNQ.KEY <- setdiff(colnames(aln), c(CTG.KEY, BAM.KEY))
    aln <- aln[,c(CTG.KEY, BAM.KEY, UNQ.KEY),with=FALSE]
    return(aln)
}

.pairAlignmentBreakpoints <- function(aln, ann) {
    tmp <- aln[order(contig_id, first_base_q)][n.seg==2]
    if (nrow(tmp)>0) {
        aln.bpt <- cbind(
            tmp[seq(1, nrow(tmp), 2), c("contig_id")],
            tmp[seq(1, nrow(tmp), 2), .(chr.5=rname, pos.5=ifelse(strand=="+", last_base_r, first_base_r),
                                        str.5=strand, mapq.5=mapq)],
            tmp[seq(2, nrow(tmp), 2), .(chr.3=rname, pos.3=ifelse(strand=="+", first_base_r, last_base_r),
                                        str.3=strand, mapq.3=mapq)]
        )
    } else {
        aln.bpt <- data.table(
            contig_id=character(),
            chr.5=character(), pos.5=integer(), str.5=character(), mapq.5=integer(),
            chr.3=character(), pos.3=integer(), str.3=character(), mapq.3=integer()
        )
    }
    aln.bpt <- aln.bpt[chr.5 %in% seqnames(ann$seqi) & chr.3 %in% seqnames(ann$seqi)]
    return(aln.bpt)
}

.filterContigTranscripts <- function(ctg, lgt, max.size=MAXINT, max.n=MAXINT) {
    lids <- unique(c(
        as.character(ctg$locus_id.5.1),
        as.character(ctg$locus_id.3.1)
    ))
    txs <- lgt[lids, nomatch=0]
    txs <- txs[nchar(tx.seq) <= max.size]
    txs[,mrna:=list(as.list(DNAStringSet(tx.seq)))]
    txs[,prot:=list(as.list(suppressWarnings(translate(subseq(DNAStringSet(tx.seq), cds1, cdsN)))))]
    txs <- txs[,.(locus_id, transcript_id, mrna, prot, cds1, cdsN, tags)]
    txs <- txs[order(tags)]
    txs <- txs[,.SD[1:min(max.n, nrow(.SD))], by=locus_id]
    return(txs)
}

.alignmentsToBreakpoints <- function(alns, ann) {
    alns <- .positionAlignments(alns)
    cbpt <- .pairAlignmentBreakpoints(alns, ann)
    cbpt <- .associateBreakpoints(cbpt, ann)
    return(cbpt)
}

.realignContigs <- function(ctgs, ann) {
    ## realign
    alns.mm2 <- .alignWrapperChain(.runMinimap2, ctgs, ann)
    cbpt.mm2 <- .alignmentsToBreakpoints(alns.mm2, ann)
    alns.gmap <- .alignWrapperChain(.runGmap, ctgs, ann)
    cbpt.gmap <- .alignmentsToBreakpoints(alns.gmap, ann)
    cbpt <- unique( rbind(cbpt.mm2, cbpt.gmap)[,.(contig_id, locus_id.5.1, locus_id.3.1)])
    setkey(ctgs, contig_id)
    setkey(cbpt, contig_id)
    ctgs <- cbpt[ctgs]
    return(ctgs)
}

.sewContigs <- function(ctgs, ann) {
    lgt <- .lgt(ann)
    txs <- .filterContigTranscripts(ctgs, lgt)
    alns <- .alignContigsToTranscripts(ctgs, txs, ann)
    chms <- .assembleChimeras(ctgs, alns, txs, ann)
    fuss <- .assembleFusions(ctgs, alns, chms, txs, ann)
    sew <- list(ctg=ctgs, aln=alns, chm=chms, fus=fuss)
    return(sew)
}

#' @export
sewBreakpoints <- function(bun, ann) {
    ## assemble
    chain.type <- paste0(bun$type, ".chain")
    buns <- splitBundle(bun, chain.type)
    ctgs <- rbindlist(mclapply(buns, function(bun) {
        ctg <- .runInchworm(bun$bpt, bun$jnc, ann)
        if (nrow(ctg)>0) {
            ctg <- cbind(bun$bpt[1,"sv.chain"], ctg)
        } else {
            ctg <- cbind(bun$bpt[0,"sv.chain"], ctg)
        }
        return(ctg)
    }))
    ## realign
    ctgs[,contig_id:=paste0(sv.chain, contig_id)]
    ctgs <- .realignContigs(ctgs, ann)
    ctgss <- split(ctgs, ctgs$sv.chain)
    ## sew
    sews <- mclapply(names(ctgss), function(chain) {
        ctgs <- ctgss[[chain]]
        sew <- .sewContigs(ctgs, ann)
        sew <- lapply(sew, function(tbl) tbl[,(chain.type):=as.integer(chain)])
        return(sew)
    })
    if (length(sews)>0) {
        nasm <- c("ctg", "aln", "chm", "fus")
        asm <- lapply(nasm, function(x) rbindlist(lapply(sews, "[[", x)))
        names(asm) <- nasm
    } else {
        asm <- ASM.EMPTY
    }
    return(asm)
}
