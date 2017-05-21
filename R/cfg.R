#' @export
makeDirectory <- function(dir.fn) {
    dir <- list(
      jnc.pe.fn = list.files(dir.fn, pattern="*pe.jnc.gz$", full.names=TRUE),
      jnc.se.fn = list.files(dir.fn, pattern="*se.jnc.gz$", full.names=TRUE),
      bam.pe.fn = list.files(dir.fn, pattern="*chim_pe_csort.cram$", full.names=TRUE),
      bam.se.fn = list.files(dir.fn, pattern="*chim_se_csort.cram$", full.names=TRUE),
      sj.fn = list.files(dir.fn, pattern="*_alig.sj.gz", full.names=TRUE),
      aln.log.fn = list.files(dir.fn, pattern="*_alig.log$", full.names=TRUE),
      bw.fn = list.files(dir.fn, "*alig_csort.bw$", full.names=TRUE),
      cts.sum.fn = list.files(dir.fn, "*cts.cts.summary$", full.names=TRUE),
      sid = basename(dir.fn)
    )
    return(dir)
}

#' @export
makeParams <- function(gtf.fn, genome="hg38", config.file="built-in", preset="longread.balanced",
                       stranded=TRUE, lib.type="poly", only.spn.bpt=TRUE, only.hx.bpt=TRUE, opts=list()) {
    ## presets
    if (config.file=="built-in") {
        config.file <- system.file("extdata", "presets.conf", package="codac")
    }
    tmp <- fread(config.file)
    PRESETS <- lapply(tmp[,-1,with=FALSE], function(col){names(col) <- tmp$param; as.list(col)})
    pres <- PRESETS[[preset]]
    pres$stranded <- stranded
    pres$lib.type <- lib.type
    pres$only.spn.bpt <- only.spn.bpt
    pres$only.hx.bpt <- only.hx.bpt
    pres$read.length <- str_split(preset, "\\.")[[1]][1]
    ## genome parameters
    par38 <- list(
      gme = "hg38",
      gtf.fn = gtf.fn,
      cut.fn = system.file("extdata", sprintf("qc-cutoffs-%s.conf", pres$lib.type), package="codac"),
      goi.fn = system.file("extdata", "hg38.genes-of-interest.txt", package="codac"),
      loi.fn = system.file("extdata", "hg38.loci-of-interest.txt", package="codac"),
      gtb.fn = system.file("extdata", "hg38.genes-to-blacklist.txt", package="codac"),
      art.fn = system.file("extdata", sprintf("hg38.artifacts-%s.csv.gz", pres$read.length), package="codac"),
      rep.fn = system.file("extdata", "hg38.rm.bed.gz", package="codac"),
      ret.fn = system.file("extdata", "hg38.retro.bed.gz", package="codac"),
      cyt.fn = system.file("extdata", "hg38.cyt", package="codac"),
      seg.fn = system.file("extdata", "hg38.segdup.gz", package="codac"),
      alt.fn = system.file("extdata", "hg38.scaffold.placement.txt", package="codac"),
      sv.fn  = system.file("extdata", "hg38.sv.1000g.v2.vcf.gz", package="codac"),
      low.fn = system.file("extdata", "hg38.dust50.bed.gz", package="codac"),
      igx.fn = system.file("extdata", "hg38.igx.bed", package="codac"),
      sno.fn = system.file("extdata", "hg38.sno.gff.gz", package="codac")
    )
    if (tolower(genome) == "hg38") {
        par <- par38
    }
    par <- c(par, pres)
    for (name in names(opts)) {
        val <- opts[[name]]
        par[[name]] <- val
    }
    par$version <- packageVersion("codac")
    return(par)
}
