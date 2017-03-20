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
makeParams <- function(gtf.fn, genome="grch38", config.file="built-in", preset="longread.balanced",
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
      gme = "grch38",
      chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
              "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"),
      chr.len = c(248956422L, 242193529L, 198295559L, 190214555L, 181538259L, 170805979L, 159345973L, 145138636L, 138394717L,
                  133797422L, 135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L, 83257441L, 80373285L,
                  58617616L, 64444167L, 46709983L, 50818468L, 156040895L, 57227415L, 16569L),
      gtf.fn = gtf.fn,
      cut.fn = system.file("extdata", sprintf("qc-cutoffs-%s.conf", pres$lib.type), package="codac"),
      goi.fn = system.file("extdata", "grch38.genes-of-interest.txt", package="codac"),
      loi.fn = system.file("extdata", "grch38.loci-of-interest.txt", package="codac"),
      gtb.fn = system.file("extdata", "grch38.genes-to-blacklist.txt", package="codac"),
      art.fn = system.file("extdata", sprintf("grch38.artifacts-%s.csv.gz", pres$read.length), package="codac"),
      rep.fn = system.file("extdata", "grch38.rm.bed.gz", package="codac"),
      ret.fn = system.file("extdata", "grch38.retro.bed.gz", package="codac"),
      cyt.fn = system.file("extdata", "grch38.cyt", package="codac"),
      seg.fn = system.file("extdata", "grch38.segdup.gz", package="codac"),
      alt.fn = system.file("extdata", "grch38.scaffold.placement.txt", package="codac"),
      sv.fn  = system.file("extdata", "grch38.sv.1000g.v2.vcf.gz", package="codac"),
      low.fn = system.file("extdata", "grch38.dust50.bed.gz", package="codac"),
      igx.fn = system.file("extdata", "grch38.igx.bed", package="codac"),
      sno.fn = system.file("extdata", "grch38.sno.gff.gz", package="codac")
    )
    if (tolower(genome) == "grch38") {
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
