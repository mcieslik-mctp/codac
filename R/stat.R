#' @export
chimeric.sig <- function(sel.bun) {
    tmp.1 <- data.table(expand.grid(dst=levels(sel.bun$bpt$dst), topo=levels(sel.bun$bpt$topo)))
    tmp.2 <- sel.bun$bpt[,.(n = uniqueN(.SD,by=CHM.KEY)), by=.(dst,topo)]
    setkey(tmp.1, dst, topo)
    setkey(tmp.2, dst, topo)
    tmp.2 <- .naTo0(tmp.2[tmp.1])
    return(tmp.2)
}
