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
      optparse::make_option(c("-k", "--keep.enc.bpt"), action="store_true",
                            default=FALSE,
                            help="keep encompassing breakpoints"),
      optparse::make_option(c("-a", "--keep.all.bpt"), action="store_true",
                            default=FALSE,
                            help="keep all breakpoints"),
      optparse::make_option(c("-g", "--genome"), type="character",
                            default="hg38",
                            help="set genome version: hg38"),
      optparse::make_option(c("-o", "--opts"), type="character",
                            default=NA_character_,
                            help="additional options")
    )
    parser = optparse::OptionParser(
      "Rscript -e 'library(methods);codac::config()' [options] gtf_file gmap_index mm2_index out_file",
      description=c("Create configuration file.\n"),
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
    if (length(opt$args) < 4) {
        optparse::print_help(parser)
        write("Some of gtf_file gmap_genome out_file are missing.\n", stderr())
        quit("no", 1)
    }
    gtf <- opt$args[1]
    gmap.index <- opt$args[2]
    mm2.index <- opt$args[3]
    out <- opt$args[4]
    par <- makeParams(gtf, gmap.index=gmap.index, mm2.index=mm2.index, genome=opt$options$genome,
                      config.file=opt$options$config,
                      preset=opt$options$preset,
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
      optparse::make_option(c("-j", "--cores"), type="integer",
                            default=8L,
                            help="set the number of cores"),
      optparse::make_option(c("-s", "--nostat"), action="store_true",
                            default=FALSE,
                            help="disable bundle stats"),
      optparse::make_option(c("-p", "--nospl"), action="store_true",
                            default=FALSE,
                            help="do not save splicing counts"),
      optparse::make_option(c("-c", "--nocts"), action="store_true",
                            default=FALSE,
                            help="do not save gene counts"),
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
      optparse::make_option(c("-a", "--noasm"), action="store_true",
                            default=FALSE,
                            help="disable assembly and realignment"),
      optparse::make_option(c("-f", "--full"), action="store_true",
                            default=FALSE,
                            help="set flag to store temporary (unfiltered) bundle")
    )
    parser = optparse::OptionParser(
      "Rscript -e 'library(methods);codac::detect()' [options] cfg_file inp_dir [out_dir]",
      description=c("Detect all types of chimeric RNAs.\n"),
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
    ## cores
    options(mc.cores=opt$options$cores)
    ## paths
    cfg.pth <- opt$args[1]
    inp.pth <- opt$args[2]
    out.dir <- ifelse(is.na(opt$args[3]), getwd(), opt$args[3])
    run <- runDetect(inp.pth, cfg.pth)
    #### output
    sid <- basename(inp.pth)
    ## STAT
    if (!opt$options$nostat) {
        saveRDS(run$stat, file.path(out.dir, paste0(sid, "-codac-stat.rds")))
    }
    ## FULL
    if (opt$options$full) {
        saveRDS(run$bun, file.path(out.dir, paste0(sid, "-codac-full.rds")))
    }
    ## SPL
    if (!opt$options$nospl) {
        saveRDS(run$spl, file.path(out.dir, paste0(sid, "-codac-spl.rds")))
    }
    ## CTS
    if (!opt$options$nocts) {
        saveRDS(run$cts, file.path(out.dir, paste0(sid, "-codac-cts.rds")))
    }
    ## SV
    if (!opt$options$nosv) {
        sv.bun <- svBundle(run$bun, run$ann)
        if (!opt$options$noasm) {
            sv.bun <- assembleBreakpoints(sv.bun, run$ann)
            sv.bun <- alignBreakpoints(sv.bun, run$ann)
            sv.bun <- validateBreakpoints(sv.bun, run$ann)
        }
        saveRDS(sv.bun, file.path(out.dir, paste0(sid, "-codac-sv.rds")))
    }
    ## BS
    if (!opt$options$nobs) {
        bs.bun <- bsBundle(run$bun, run$ann)
        saveRDS(bs.bun, file.path(out.dir, paste0(sid, "-codac-bs.rds")))
    }
    ## SL
    if (!opt$options$nosl) {
        sl.bun <- slBundle(run$bun, run$ann)
        saveRDS(sl.bun, file.path(out.dir, paste0(sid, "-codac-sl.rds")))
    }
    ## TS
    if (!opt$options$nots) {
        ts.bun <- tsBundle(run$bun, run$ann)
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
    parser = optparse::OptionParser(
      "Rscript -e 'library(methods);codac::assemble()' [options] cfg_file inp_file [out_dir]",
      description=c("Assemble detected chimeric breakpoints.\n"),
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
    options(mc.cores=opt$options$cores)
    sew <- sewBreakpoints(bun, ann)
    out.pth <- file.path(out.dir, paste0(name, "-asm.rds"))
    saveRDS(sew, out.pth)
}

#' @export
report <- function() {
    option_list = list(
      optparse::make_option(c("-n", "--noncoding"), action="store_true",
                            default=FALSE,
                            help="include non-coding gene names in output")
    )
    parser = optparse::OptionParser(
      "Rscript -e 'library(methods);codac::report()' [options] cfg_file spl_file cts_file inp_file [out_dir]",
      description=c("Produce chimeric RNA reports in CSV format.\n"),
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
    cts.pth <- opt$args[3]
    inp.pth <- opt$args[4]
    out.dir <- ifelse(is.na(opt$args[5]), getwd(), opt$args[5])
    ##
    ann <- readRDS(cfg.pth)
    spl <- readRDS(spl.pth)
    cts <- readRDS(cts.pth)
    bun <- readRDS(inp.pth)
    name <- str_replace(basename(inp.pth), ".rds$", "")
    if (bun$type == "bs") {
        rep <- bsFormat(bsReport(bun, spl, cts, ann, !opt$options$noncoding))
    } else if (bun$type == "sv") {
        rep <- svFormat(svReport(bun, spl, cts, ann, !opt$options$noncoding))
    } else if (bun$type == "sl") {
        rep <- slFormat(slReport(bun, spl, cts, ann, !opt$options$noncoding))
    } else if (bun$type == "ts") {
        rep <- tsFormat(tsReport(bun, spl, cts, ann, !opt$options$noncoding))
    }
    else {
        write("Unknown bundle type.\n", stderr())
        quit("no", 1)
    }
    fwrite(rep, file.path(out.dir, paste0(name, "-rep.csv")), row.names=FALSE, quote=FALSE)
}

#' @export
qc.report <- function() {
    option_list = list(
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::qc.report()' [options] cfg_file stat_file [out_dir]",
      description=c("Produce QC report in CSV format.\n"),
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
    ##
    cfg.pth <- opt$args[1]
    stat.pth <- opt$args[2]
    out.dir <- ifelse(is.na(opt$args[3]), getwd(), opt$args[3])
    stat <- readRDS(stat.pth)
    stat.fmt <- qcFormat(stat)
    name <- str_replace(basename(stat.pth), ".rds$", "")
    out.fn <- file.path(out.dir, paste0(name, "-rep.csv"))
    fwrite(stat.fmt, out.fn, row.names=FALSE, quote=FALSE)
}

#' @export
neo.report <- function() {
    option_list = list(
    )
    parser = optparse::OptionParser("Rscript -e 'library(methods);codac::neo.report()' [options] cfg_file asm_file [out_dir]",
      description=c("Produce neoantigen report in CSV format.\n"),
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
    ##
    cfg.pth <- opt$args[1]
    asm.pth <- opt$args[2]
    out.dir <- ifelse(is.na(opt$args[3]), getwd(), opt$args[3])
    ##
    ann <- readRDS(cfg.pth)
    asm <- readRDS(asm.pth)
    neo.rep <- neoReport(asm, ann)
    neo.fmt <- neoFormat(neo.rep)
    name <- str_replace(basename(asm.pth), ".rds$", "")
    out.fn <- file.path(out.dir, paste0(name, "-neo.csv"))
    fwrite(neo.fmt, out.fn, row.names=FALSE, quote=FALSE)
}
