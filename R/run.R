#' @export
runDetect <- function(inp.pth, cfg.pth) {
    beg.time <- Sys.time()
    ##
    ann <- readRDS(cfg.pth)
    ## annotation / configuration
    if (packageVersion("codac")[,1:2] != ann$par$version[,1:2]) {
        write("Configuration and program version mismatch.\n", stderr())
        quit("no", 1)
    }
    dir <- makeDirectory(inp.pth)
    cts <- readCounts(dir, ann)
    spl <- readSplices(dir, ann)
    jnc <- readJunctions(dir, ann)
    bpt <- collapseBreakpoints(jnc, ann)
    stat <- makeStat(jnc, bpt, spl, dir, ann)
    bun <- makeBundle(bpt, jnc, "full")
    end.time <- Sys.time()
    stat$log$beg.time <- beg.time
    stat$log$end.time <- end.time
    stat$log$gc <- gc()
    run <- list(ann=ann, dir=dir, cts=cts, spl=spl, stat=stat, bun=bun)
    return(run)
}
