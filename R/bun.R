#' @export
makeBundle <- function(ann, pth) {
    dir <- makeDirectory(pth)
    spl <- readSplices(dir, ann)
    jnc <- readJunctions(dir, ann)
    bpt <- collapseBreakpoints(jnc, ann)
    bun <- list(bpt=bpt, jnc=jnc, spl=spl, dir=dir, type="full")
    return(bun)
}

#' @export
statBundle <- function(bun, ann) {
    aln.stat <- .alignment.stat(bun, ann)
    bpt.stat <- .breakpoint.stat(bun, aln.stat$eff.frg)
    loc.stat <- .locus.stat(bun, ann)
    log.stat <- list(version=packageVersion("codac"), beg.time=NULL, end.time=NULL, gc=NULL)
    bun.stat <- c(aln.stat, bpt.stat, loc.stat)
    bun.stat$log <- log.stat
    return(bun.stat)
}

#' @export
filterBundle <- function(sel.bpt, bun, ann) {
    ## get junctions
    setkeyv(sel.bpt, BPT.KEY)
    setkeyv(bun$jnc, BPT.KEY)
    sel.jnc <- bun$jnc[sel.bpt[, BPT.KEY, with=FALSE], nomatch=0]
    setkeyv(bun$jnc, JNC.KEY)
    setkeyv(sel.jnc, JNC.KEY)
    ## output bundle
    sel.bun <- list(bpt=sel.bpt, jnc=sel.jnc, spl=NULL, cts=NULL, dir=bun$dir, type=bun$type)
    return(sel.bun)
}

#' @export
splitBundle <- function(bun, col) {
    split.bpt <- split(bun$bpt, bun$bpt[[col]])
    split.bun <- lapply(split.bpt, filterBundle, bun=bun, ann=ann)
    return(split.bun)
}

#' @export
importBundles <- function(smp.pth, all.sfx=c("bs", "sl", "sv", "ts", "sv-asm")) {
    all.sfx <- c("spl", "stat", all.sfx)
    rds.fns <- list.files(smp.pth, "*.rds", full.names=TRUE)
    rds.pfx <- str_match(basename(rds.fns), "(.*)-codac")[,2]
    rds.fns.split <- split(rds.fns, rds.pfx)
    buns <- lapply(rds.fns.split, function(fns) {
        sfx <- str_match(basename(fns), "-codac-([^.]*)")[,2]
        names(fns) <- sfx
        names(all.sfx) <- all.sfx
        res <- lapply(all.sfx, function(sfx) if (is.na(fns[sfx])) NULL else readRDS(fns[sfx]))
        return(res)
    })
    return(buns)
}
