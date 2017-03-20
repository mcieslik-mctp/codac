#' @export
bsBundle <- function(bun, ann) {
    bs.bpt <- bun$bpt[(bs.chain > 0)]
    bs.bun <- filterBundle(bs.bpt, bun, ann)
    bs.bun$type <- "bs"
    return(bs.bun)
}

#' @export
bsReport <- function(bs.bun, spl, ann, only.prot) {
    bs.rep <- .reportMini(bs.bun, spl, ann, TRUE, only.prot)
    return(bs.rep)
}

#' @export
slBundle <- function(bun, ann) {
    sl.bpt <- bun$bpt[(sl.chain>0)]
    sl.bun <- filterBundle(sl.bpt, bun, ann)
    sl.bun$type <- "sl"
    return(sl.bun)
}

#' @export
slReport <- function(sl.bun, spl, ann, only.prot) {
    sl.rep <- .reportMini(sl.bun, spl, ann, FALSE, only.prot)
    return(sl.rep)
}

#' @export
svBundle <- function(bun, ann) {
    sv.bpt <- bun$bpt[(sv.chain>0)]
    sv.bun <- filterBundle(sv.bpt, bun, ann)
    sv.bun$type <- "sv"
    return(sv.bun)
}

#' @export
svReport <- function(sv.bun, spl, ann, only.prot) {
    sv.rep <- .reportSV(sv.bun, spl, ann, only.prot)
    return(sv.rep)
}

#' @export
tsBundle <- function(bun, ann) {
    ts.bpt <- bun$bpt[(ts.chain>0)]
    ts.bun <- filterBundle(ts.bpt, bun, ann)
    ts.bun$type <- "ts"
    return(ts.bun)
}

#' @export
tsReport <- function(ts.bun, spl, ann, only.prot) {
    ts.rep <- .reportMini(ts.bun, spl, ann, TRUE, only.prot)
    return(ts.rep)
}
