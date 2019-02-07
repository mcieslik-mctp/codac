#' @export
makeBundle <- function(bpt, jnc, type) {
    bun <- list(bpt=bpt, jnc=jnc, type="full")
    return(bun)
}

#' @export
filterBundle <- function(sel.bpt, bun, ann) {
    setkeyv(sel.bpt, BPT.KEY)
    sel.bun <- lapply(bun, function(x) NULL)
    ## preserve required
    sel.bun$bpt <- sel.bpt
    sel.bun$type <- bun$type
    ## filter junctions, contigs, and alignments ...
    for (key in setdiff(names(bun), c("type", "bpt"))) {
        setkeyv(bun[[key]], BPT.KEY)
        sel.bun[[key]] <- bun[[key]][sel.bpt[, BPT.KEY, with=FALSE], nomatch=0]
    }
    return(sel.bun)
}

#' @export
splitBundle <- function(bun, split.col=NULL) {
    if (is.null(split.col)) {
        split.val <- bun$bpt[,.I]
    } else {
        split.val <- bun$bpt[[split.col]]
    }
    split.bpt <- split(bun$bpt, split.val)
    split.bun <- lapply(split.bpt, filterBundle, bun=bun, ann=ann)
    return(split.bun)
}

#' @export
mergeBundles <- function(buns) {
    mrg.bun <- lapply(buns[[1]], function(x) NULL)
    ## preserve type
    mrg.bun$type <- buns[[1]]$type
    ## merge breakpoints, junctions, contigs, and alignments
    for (key in setdiff(names(mrg.bun), c("type"))) {
            mrg.key <- rbindlist(lapply(buns, "[[", key))
            setkeyv(mrg.key, BPT.KEY)
            mrg.bun[[key]] <- mrg.key
    }
    return(mrg.bun)
}

#' @export
importBundles <- function(run.pth, bun.sfx=c("bs", "sl", "sv", "ts")) {
    rds.fns <- list.files(run.pth, "*.rds", full.names=TRUE)
    names(rds.fns) <- str_match(basename(rds.fns), "-codac-([^.]*)")[,2]
    buns <- lapply(bun.sfx, function(sfx) {
        tryCatch(readRDS(rds.fns[sfx]), warning=function(w) NULL)
    })
    names(buns) <- bun.sfx
    return(buns)
}
