#' @export
config <- function() {
    option_list = list(
      optparse::make_option(c("-c", "--config"), type="character",
                            default="built-in",
                            help="preset"),
      optparse::make_option(c("-p", "--preset"), type="character",
                            default="longread.balanced",
                            help="preset"),
      optparse::make_option(c("-u", "--unstranded"), action="store_true",
                            default=FALSE,
                            help="set if library is unstranded"),
      optparse::make_option(c("-l", "--libtype"), type="character",
                            default="poly",
                            help="library type: poly, capt, ribo"),
      optparse::make_option(c("-k", "--keep.enc.bpt"), action="store_true",
                            default=FALSE,
                            help="keep encompassing breakpoints"),
      optparse::make_option(c("-a", "--keep.all.bpt"), action="store_true",
                            default=FALSE,
                            help="keep all breakpoints"),
      optparse::make_option(c("-g", "--genome"), type="character",
                            default="grch38",
                            help="set genome version: grch38"),
      optparse::make_option(c("-o", "--opts"), type="character",
                            default=NA_character_,
                            help="additional options")
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::config()' [options] gtf_file out_file",
      description=c("\n"),
      epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2017\n"),
      option_list=option_list
    )
    
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    ## numeric parameters
    if(!is.na(opt$options$opts)) {
        tmp <- str_split(str_split(opt$options$opts,",")[[1]], ":")
        tmp.n <- lapply(tmp, "[", 1)
        opts <- lapply(lapply(tmp, "[", 2), as.numeric)
        names(opts) <- tmp.n
    } else {
        opts <- list()
    }
    ##
    if (length(opt$args) < 1) {
        optparse::print_help(parser)
        write("gtf_file is missing.\n", stderr())
        quit("no", 1)
    }
    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("out_file is missing.\n", stderr())
        quit("no", 1)
    }
    gtf <- opt$args[1]
    out <- opt$args[2]
    gme <- opt$options$genome
    par <- makeParams(gtf, genome=gme, config.file=opt$options$config,
                      preset=opt$options$preset,
                      lib.type=opt$options$libtype,
                      stranded=!opt$options$unstranded,
                      only.spn.bpt=!opt$options$keep.enc.bpt,
                      only.hx.bpt=!opt$options$keep.all.bpt,
                      opts)
    ann <- makeAnnotations(par)
    saveRDS(ann, out)
}

#' @export
detect <- function() {
    option_list = list(
      optparse::make_option(c("-s", "--nostat"), action="store_true",
                            default=FALSE,
                            help="disable bundle stats"),
      optparse::make_option(c("-p", "--nospl"), action="store_true",
                            default=FALSE,
                            help="disable splicing report"),
      optparse::make_option(c("-v", "--nosv"), action="store_true",
                            default=FALSE,
                            help="disable structural variant calls"),
      optparse::make_option(c("-b", "--nobs"), action="store_true",
                            default=FALSE,
                            help="disable back-splice calls"),
      optparse::make_option(c("-l", "--nosl"), action="store_true",
                            default=FALSE,
                            help="disable stem-loop calls"),
      optparse::make_option(c("-t", "--nots"), action="store_true",
                            default=FALSE,
                            help="disable trans-splice calls"),
      optparse::make_option(c("-f", "--full"), action="store_true",
                            default=FALSE,
                            help="set flag to store temporary (unfiltered) bundle")
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::detect()' [options] cfg_file inp_dir [out_dir]",
      description=c("\n"),
      epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2017\n"),
      option_list=option_list
    )
    
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    ##
    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("required arguments missing missing.\n", stderr())
        quit("no", 1)
    }
    ## paths
    cfg.pth <- opt$args[1]
    inp.pth <- opt$args[2]
    out.dir <- ifelse(is.na(opt$args[3]), getwd(), opt$args[3])
    sid <- basename(inp.pth)
    ## annotation / configuration
    ann <- readRDS(cfg.pth)
    if (packageVersion("codac")[,1:2] != ann$par$version[,1:2]) {
        write("Configuration and program version mismatch.\n", stderr())
        quit("no", 1)
    }
    ## create bundle / stats
    beg.time <- Sys.time()
    bun <- makeBundle(ann, inp.pth)
    end.time <- Sys.time()
    ## FULL
    if (opt$options$full) {
        saveRDS(bun, file.path(out.dir, paste0(sid, "-codac-sv-full.rds")))
    }
    ## STAT
    if (!opt$options$nostat) {
        stat <- statBundle(bun, ann)
        stat$log$beg.time <- beg.time
        stat$log$end.time <- end.time
        stat$log$gc <- gc()
        saveRDS(stat, file.path(out.dir, paste0(sid, "-codac-stat.rds")))
    }
    ## SPL
    if (!opt$options$nospl) {
        saveRDS(bun$spl, file.path(out.dir, paste0(sid, "-codac-spl.rds")))
    }
    ## BS
    if (!opt$options$nobs) {
        bs.bun <- bsBundle(bun, ann)
        saveRDS(bs.bun, file.path(out.dir, paste0(sid, "-codac-bs.rds")))
    }
    ## SL
    if (!opt$options$nosl) {
        sl.bun <- slBundle(bun, ann)
        saveRDS(sl.bun, file.path(out.dir, paste0(sid, "-codac-sl.rds")))
    }
    ## SV
    if (!opt$options$nosv) {
        sv.bun <- svBundle(bun, ann)
        saveRDS(sv.bun, file.path(out.dir, paste0(sid, "-codac-sv.rds")))
    }
    ## TS
    if (!opt$options$nots) {
        ts.bun <- tsBundle(bun, ann)
        saveRDS(ts.bun, file.path(out.dir, paste0(sid, "-codac-ts.rds")))
    }

}

#' @export
assemble <- function() {
    option_list = list(
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=8L,
                              help="set the number of cores")
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::assemble()' [options] cfg_file inp_file [out_dir]",
      description=c("Assemble detected breakpoints\n"),
      epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2017\n"),
      option_list=option_list
    )
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    ## input check
    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("required arguments missing missing.\n", stderr())
        quit("no", 1)
    }
    cfg.pth <- opt$args[1]
    inp.pth <- opt$args[2]
    out.dir <- ifelse(is.na(opt$args[3]), getwd(), opt$args[3])
    ann <- readRDS(cfg.pth)
    bun <- readRDS(inp.pth)
    name <- str_replace(basename(inp.pth), ".rds$", "")
    options(mc.cores=opt$options$j)
    asm <- assembleBundle(bun, ann)
    out.pth <- file.path(out.dir, paste0(name, "-asm.rds"))
    saveRDS(asm, out.pth)
}

#' @export
report <- function() {
    option_list = list(
      optparse::make_option(c("-n", "--noncoding"), action="store_true",
                            default=FALSE,
                            help="include non-coding gene names in output")
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::report()' [options] cfg_file spl_file bun_file [out_dir]",
      description=c("Export \n"),
      epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2017\n"),
      option_list=option_list
    )
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    ## input check
    if (length(opt$args) < 3) {
        optparse::print_help(parser)
        write("required arguments missing missing.\n", stderr())
        quit("no", 1)
    }
    ##
    cfg.pth <- opt$args[1]
    spl.pth <- opt$args[2]
    inp.pth <- opt$args[3]
    out.dir <- ifelse(is.na(opt$args[4]), getwd(), opt$args[4])
    ##
    ann <- readRDS(cfg.pth)
    spl <- readRDS(spl.pth)
    bun <- readRDS(inp.pth)
    name <- str_replace(basename(inp.pth), ".rds$", "")
    if (bun$type == "bs") {
        rep <- bsFormat(bsReport(bun, spl, ann, !opt$options$noncoding))
    } else if (bun$type == "sv") {
        rep <- svFormat(svReport(bun, spl, ann, !opt$options$noncoding))
    } else if (bun$type == "sl") {
        rep <- slFormat(slReport(bun, spl, ann, !opt$options$noncoding))
    } else if (bun$type == "ts") {
        rep <- tsFormat(tsReport(bun, spl, ann, !opt$options$noncoding))
    }
    else {
        write("Unknown bundle type.\n", stderr())
        quit("no", 1)
    }
    write.csv(rep, file.path(out.dir, paste0(name, "-rep.csv")), row.names=FALSE, quote=FALSE)
}

#' @export
qc.report <- function() {
    option_list = list(
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::qc.report()' [options] stat_file [out_dir]",
      description=c("Export \n"),
      epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2017\n"),
      option_list=option_list
      )
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    ## input check
    if (length(opt$args) < 1) {
        optparse::print_help(parser)
        write("required arguments missing missing.\n", stderr())
        quit("no", 1)
    }
    ##
    cfg.pth <- opt$args[1]
    stat.pth <- opt$args[2]
    out.dir <- ifelse(is.na(opt$args[3]), getwd(), opt$args[3])
    ann <- readRDS(cfg.pth)
    stat <- readRDS(stat.pth)
    stat.fmt <- qcFormat(stat, ann)
    name <- str_replace(basename(stat.pth), ".rds$", "")
    out.fn <- file.path(out.dir, paste0(name, "-rep.csv"))
    write.csv(stat.fmt, out.fn, row.names=FALSE, quote=FALSE)
}
