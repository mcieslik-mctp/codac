CTG.EMPTY = data.table(contig_id = character(), ctg.seq=list(), ctg.cov=numeric(), ctg.len=integer())

.writePairFastas <- function(bpt, jnc) {
    setkeyv(bpt, BPT.KEY)
    setkeyv(jnc, BPT.KEY)
    tmp <- jnc[bpt, .(qname, seq.5, seq.3)]
    fa5 <- DNAStringSet(tmp$seq.5)
    fa3 <- reverseComplement(DNAStringSet(tmp$seq.3))
    fa53 <- c(fa5, fa3)
    names(fa53) <- c(paste0(tmp$qname, "_5"), paste0(tmp$qname, "_3"))
    fn <- tempfile()
    writeXStringSet(fa53, fn)
    return(fn)
}

.runInchworm <- function(bpt, jnc, ann) {
    cfn <- .writePairFastas(bpt, jnc)
    ofn <- paste0(cfn, ".iworm")
    err <- system2("inchworm", c("--reads", cfn, "--run_inchworm"), stdout=ofn, stderr=NULL)
    out <- readDNAStringSet(ofn)
    if (err==0 && length(out)>0) {
        ctg <- data.table(
            contig_id = str_match(names(out), "(.*);")[,2],
            ctg.seq = as.list(out),
            ctg.cov = as.numeric(str_match(names(out), ";([^ ]*)")[,2]),
            ctg.len = lengths(out))
        ctg <- ctg[(ctg.cov >= ann$par$asm.ctg.cov) | (ctg.len >= ann$par$asm.ctg.len)]
    } else {
        ctg <- CTG.EMPTY
    }
    unlink(ofn)
    unlink(cfn)
    setkey(ctg, contig_id)
    return(ctg)
}

.assembleBreakpoint <- function(bun, ann) {
    ctg <- .runInchworm(bun$bpt, bun$jnc, ann)
    if (nrow(ctg)>0) {
        ctg <- cbind(bun$bpt[1,..BPT.KEY], ctg)
    } else {
        ctg <- cbind(bun$bpt[0,..BPT.KEY], ctg)
    }
    return(ctg)
}

#' @export
assembleBreakpoints <- function(bun, ann) {
    buns <- splitBundle(bun)
    ctg <- rbindlist(mclapply(names(buns), function(idx) {
        sel.bun <- buns[[idx]]
        sel.ctg <- .assembleBreakpoint(sel.bun, ann)
        return(sel.ctg)
    }))
    ctg[,contig_id:=paste0(.GRP, contig_id), by=BPT.KEY]
    bun$ctg <- ctg
    return(bun)
}
