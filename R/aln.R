BAM.EMPTY.GMAP <- data.table(
    contig_id=character(),flag=integer(),rname=character(),strand=character(),pos=integer(),
    qwidth=integer(),mapq=integer(),cigar=character(),NM=integer(),SM=integer(),XQ=integer(),
    X2=integer(),XO=character(),NH=integer(),XT=character()
)

BAM.EMPTY.MINIMAP2 <- data.table(
    contig_id=character(),flag=integer(),rname=character(),strand=character(),pos=integer(),
    qwidth=integer(),mapq=integer(),cigar=character(),NM=integer(),ms=integer(),AS=integer(),
    nn=integer(),tp=character(),cm=integer(),s1=integer(),s2=integer(),ts=character(),SA=character()
)


.writeFastas <- function(ctg) {
    ctg.seq <- DNAStringSet(ctg$ctg.seq)
    names(ctg.seq) <- ctg$contig_id
    fn <- tempfile()
    writeXStringSet(ctg.seq, fn)
    return(fn)
}

.readBam <- function(bamFile, blanks) {
    bam.tag <- names(blanks)
    bam <- scanBam(bamFile,
                   param=ScanBamParam(what=c("qname", "flag", "rname", "strand", "pos",
                                             "qwidth", "mapq", "cigar"),
                                      tag=bam.tag
                                      ))[[1]]
    tbl.q <- as.data.table(bam[!(names(bam) %in% c("tag"))])
    tbl.q[,":="(rname=as.character(rname), strand=as.character(strand))]
    tbl.t <- as.data.table(bam$tag)
    for (tag in names(blanks)) {
        if (is.null(tbl.t[[tag]])) {
            tbl.t[,eval(tag):=blanks[[tag]]]
        }
    }
    tbl.t <- tbl.t[,names(blanks),with=FALSE]
    res <- cbind(tbl.q, tbl.t)
    setnames(res, "qname", "contig_id")
    setkey(res, contig_id)
    return(res)
}

.runGmap <- function(ctg, ann) {
    if (nrow(ctg)>0) {
        gfn <- .writeFastas(ctg)
        bamFile <- paste0(gfn, ".gmap")
        OPT <- sprintf("-f samse -n0 -D %s -d %s -x %s -K %d",
                       dirname(ann$par$gmap.index), basename(ann$par$gmap.index),
                       ann$par$asm.gmap.minchi, ann$par$asm.gmap.maxint)
        PIPE <- sprintf("%s 2> /dev/null | samtools view -b > %s", gfn, bamFile)
        err <- system2("gmap", c(OPT, PIPE), stdout=TRUE)
        blanks <- list(
            NM=NA_integer_,SM=NA_integer_,XQ=NA_integer_,
            X2=NA_integer_,XO=NA_character_,NH=NA_integer_,
            XT=NA_character_
        )
        bam <- .readBam(bamFile, blanks)
        unlink(bamFile)
        unlink(gfn)
    } else {
        bam <- BAM.EMPTY.GMAP
    }
    return(bam)
}

.runMinimap2 <- function(ctg, ann) {
    if (nrow(ctg)>0) {
        gfn <- .writeFastas(ctg)
        bamFile <- paste0(gfn, ".m2")
        OPT <- sprintf("-t %s -a -x splice %s", getOption("mc.cores"), ann$par$mm2.index)
        PIPE <- sprintf("%s 2> /dev/null | samtools view -F256 -b > %s", gfn, bamFile)
        err <- system2("minimap2", c(OPT, PIPE), stdout=TRUE)
        blanks <- list(
            NM=NA_integer_,ms=NA_integer_,AS=NA_integer_,
            nn=NA_integer_,tp=NA_character_,cm=NA_integer_,
            s1=NA_integer_,s2=NA_integer_,ts=NA_character_,
            SA=NA_character_
        )
        bam <- .readBam(bamFile, blanks)
        unlink(bamFile)
        unlink(gfn)
    } else {
        bam <- BAM.EMPTY.MINIMAP2
    }
    return(bam)
}

.alignWrapperBreakpoint <- function(f, ctg, ann) {
    bam <- f(ctg, ann)
    setkey(ctg, contig_id)
    setkey(bam, contig_id)
    aln <- bam[ctg]
    CTG.KEY <- colnames(CTG.EMPTY)
    BAM.KEY <- setdiff(intersect(colnames(BAM.EMPTY.MINIMAP2), colnames(BAM.EMPTY.GMAP)), CTG.KEY)
    UNQ.KEY <- setdiff(colnames(aln), c(BPT.KEY, CTG.KEY, BAM.KEY))
    aln <- aln[,c(BPT.KEY, CTG.KEY, BAM.KEY, UNQ.KEY),with=FALSE]
    return(aln)
}

.positionAlignments <- function(aln) {
    if (nrow(aln)>0) {
        ## compute cliping
        clip.l <- as.integer(str_match(aln$cigar, "^([0-9]+)[SH]")[,2])
        clip.l[!is.na(aln$cigar) & is.na(clip.l)] <- 0L
        clip.r <- as.integer(str_match(aln$cigar, "([0-9]+)[SH]$")[,2])
        clip.r[!is.na(aln$cigar) & is.na(clip.r)] <- 0L
        ## update alignment
        aln[,n.seg:=.N, by=contig_id]
        aln[,aligned_bases:=cigarWidthAlongQuerySpace(cigar, after.soft.clipping = TRUE)]
        aln[,first_base_r:=pos]
        aln[,first_base_q:=ifelse(aln$strand=="+", clip.l+1L, clip.r+1L)]
        aln[,last_base_r:=first_base_r + cigarWidthAlongReferenceSpace(cigar)]
        aln[,last_base_q:=first_base_q + aligned_bases]
    } else {
        aln[,":="(n.seg=integer(), aligned_bases=integer(),
                  first_base_r=integer(), first_base_q=integer(),
                  last_base_r=integer(), last_base_q=integer())]
    }
    return(aln)
}

#' @export
alignBreakpoints <- function(bun, ann) {
    aln.mm2 <- .alignWrapperBreakpoint(.runMinimap2, bun$ctg, ann)
    aln.gmap <- .alignWrapperBreakpoint(.runGmap, bun$ctg, ann)
    aln.mm2 <- .positionAlignments(aln.mm2)
    aln.gmap <- .positionAlignments(aln.gmap)
    bun$aln.mm2 <- aln.mm2
    bun$aln.gmap <- aln.gmap
    return(bun)    
}
