.getREC <- function(fn) {
    tmp <- tempfile()
    system2("zcat", sprintf("%s > %s", fn, tmp), stdout=TRUE)
    tmp.tbl <- tryCatch(fread(tmp, header=FALSE, showProgress=FALSE), error=function(e)
        data.table(matrix(nrow=0,ncol=length(REC.FIELDS))))
    REC.FIELDS <- c("chr.5", "pos.5", "str.5", "chr.3", "pos.3", "str.3", "type", "rep.5", "rep.3", "qname",
                    "seg.5.pos", "seg.5.cig", "seg.3.pos", "seg.3.cig")
    setnames(tmp.tbl, REC.FIELDS)
    unlink(tmp)
    return(tmp.tbl)
}

.getBam <- function(fn, pipe, nfields) {
    tmp <- tempfile()
    system2("samtools", c("view", paste(sprintf(pipe, fn), ">", tmp)), stdout=TRUE)
    tmp.tbl <- tryCatch(fread(tmp, header=FALSE, showProgress=FALSE), error=function(e)
        data.table(matrix(nrow=0,ncol=nfields)))
    unlink(tmp)
    return(tmp.tbl)
}

.readREC <- function(dir) {
    ## import 
    jnc.se.rec <- .getREC(dir$jnc.se.fn)
    jnc.se.rec[,ftype:="se",]
    jnc.pe.rec <- .getREC(dir$jnc.pe.fn)
    jnc.pe.rec[,ftype:="pe",]
    jnc.rec <- rbind(jnc.se.rec, jnc.pe.rec)
    rm("jnc.se.rec", "jnc.pe.rec")
    ## clipping
    jnc.rec[,
            ":="(
              clp.5.l=as.integer(str_match(seg.5.cig, "^([0-9]+)S")[,2]),
              clp.5.r=as.integer(str_match(seg.5.cig, "([0-9]+)S$")[,2]),
              clp.3.l=as.integer(str_match(seg.3.cig, "^([0-9]+)S")[,2]),
              clp.3.r=as.integer(str_match(seg.3.cig, "([0-9]+)S$")[,2])
            )]
    jnc.rec <- .naTo0(jnc.rec)
    setkeyv(jnc.rec, "qname")
    return(jnc.rec)
}

.readTAG <- function(dir) {
    tmp.se <- .getBam(dir$bam.se.fn, "%s | cut -f 1,13,14,15", 4)
    tmp.se$V5 <- .getBam(dir$bam.se.fn, "%s | cut -f 10 | awk '{print length}'", 1)$V1
    tmp <- tmp.se[,list(V3[1],V4[1],V5[1]),by=list(V1,V2)][order(V1, V2)]
    tmp <- suppressWarnings(tmp[,N:=.N,by=V1])[N==2][,1:(ncol(tmp)-1),with=FALSE]
    bam.se.rec <- cbind(tmp[seq(1,nrow(tmp),2)]$V1,
                        tmp[seq(1,nrow(tmp),2), 3:5, with=FALSE],
                        tmp[seq(2,nrow(tmp),2), 3:5, with=FALSE])
    tmp.pe <- .getBam(dir$bam.pe.fn, "%s | cut -f 1,13,14,15", 4)
    tmp.pe$V5 <- .getBam(dir$bam.pe.fn, "%s | cut -f 10 | awk '{print length}'", 1)$V1
    tmp <- tmp.pe[,list(V3[1],V4[1],V5[1]),by=list(V1,V2)][order(V1, V2)]
    tmp <- tmp[,N:=.N,by=V1][N==2][,1:(ncol(tmp)-1),with=FALSE]
    bam.pe.rec <- cbind(tmp[seq(1,nrow(tmp),2)]$V1,
                        tmp[seq(1,nrow(tmp),2), 3:5, with=FALSE],
                        tmp[seq(2,nrow(tmp),2), 3:5, with=FALSE])
    bam.rec <- rbind(bam.se.rec, bam.pe.rec)
    rm("bam.se.rec", "bam.pe.rec")
    BAM.FIELDS <- c("qname", "seg.5.as", "seg.5.nm", "seg.5.len", "seg.3.as", "seg.3.nm", "seg.3.len")
    setnames(bam.rec, BAM.FIELDS)
    ## fread... is faster than str_match
    bam.rec$seg.5.as <- fread(paste(bam.rec$seg.5.as, collapse="\n"), head=FALSE, sep=":", showProgress=FALSE)$V3
    bam.rec$seg.3.as <- fread(paste(bam.rec$seg.3.as, collapse="\n"), head=FALSE, sep=":", showProgress=FALSE)$V3
    bam.rec$seg.5.nm <- fread(paste(bam.rec$seg.5.nm, collapse="\n"), head=FALSE, sep=":", showProgress=FALSE)$V3
    bam.rec$seg.3.nm <- fread(paste(bam.rec$seg.3.nm, collapse="\n"), head=FALSE, sep=":", showProgress=FALSE)$V3
    ## estimated overlap
    bam.rec[, ":="(
        ovr.5=seg.5.as + (2*seg.5.nm),
        ovr.3=seg.3.as + (2*seg.3.nm)
                                  
    )]
    ## estimated alignment error
    bam.rec[, ":="(
        err.5=seg.5.nm / ovr.5,
        err.3=seg.3.nm / ovr.3
    )]
    setkeyv(bam.rec, "qname")
    return(bam.rec)
}

.readSEQ <- function(dir) {
    seq.pe <- .getBam(dir$bam.pe.fn, "%s | cut -f 1,2,10,13", 4)
    setnames(seq.pe, c("qname", "flag", "seq", "tag"))
    seq.pe[,n:=.N,by="qname"] # replacement for seq.pe <- seq.pe[,.SD[1],by=list(qname,tag)]
    seq.pe.2 <- seq.pe[n==2]
    seq.pe.3 <- seq.pe[n==3]
    setkey(seq.pe.3, "qname", "tag", "flag")
    row1and3 <- setdiff(seq_len(nrow(seq.pe.3)), seq(2, nrow(seq.pe.3), 3))
    seq.pe.3 <- seq.pe.3[row1and3]
    seq.pe <- rbind(seq.pe.2, seq.pe.3)[,list(qname, tag, seq)]
    seq.se <- .getBam(dir$bam.se.fn, "%s | cut -f 1,13,10", 3)
    setnames(seq.se, c("qname", "seq", "tag"))
    setkey(seq.se, "qname", "tag")
    seq.se <- seq.se[,.SD[1],by=list(qname,tag)]
    seq.rec <- rbind(seq.pe, seq.se)
    seq.rec$tag <- factor(seq.rec$tag, levels=c("HI:i:1", "HI:i:2"))
    setkeyv(seq.rec, "qname")
    return(seq.rec)
}

.mergeJNC <- function(dir, ann) {
    rec <- .readREC(dir)
    tag <- .readTAG(dir)
    seq <- .readSEQ(dir)
    jnc <- merge(rec, tag)
    rm(rec)
    rm(tag)
    jnc[seq[tag=="HI:i:1"], seq.5:=seq]
    jnc[seq[tag=="HI:i:2"], seq.3:=seq]
    rm(seq)
    jnc <- jnc[!is.na(seq.5) & !is.na(seq.3)]
    jnc <- jnc[chr.5 %in% ann$par$chr & chr.3 %in% ann$par$chr]
    jnc.5 <- GRanges(jnc$chr.5, IRanges(jnc$pos.5, width=1),jnc$str.5)
    jnc.3 <- GRanges(jnc$chr.3, IRanges(jnc$pos.3, width=1),jnc$str.3)
    jnc.ok <- (jnc.5 %over% ann$loci) & (jnc.3 %over% ann$loci)
    jnc <- jnc[jnc.ok]
    jnc[, ":="(chr.5=factor(chr.5, levels=ann$par$chr, ordered=TRUE),
               chr.3=factor(chr.3, levels=ann$par$chr, ordered=TRUE),
               str.5=factor(str.5, levels=c("+", "-"), ordered=TRUE),
               str.3=factor(str.3, levels=c("+", "-"), ordered=TRUE)
               )]
    return(jnc)
}

.dustJunctions <- function(jnc) {
    ## maximim possible dusty score for a given sequence length
    max.low <- data.table(
      len=4:1000,
      low=dustyScore(DNAStringSet(sapply(4:1000, function(x) paste0(rep("A", x), collapse="")))))
    setkey(max.low, len)
    ## 1. split dustyScore into chunks to save some memory?
    tmp.tbl <- jnc[,.(
      qname,
      len.5=nchar(seq.5)-clp.5.l-clp.5.r,
      len.3=nchar(seq.3)-clp.3.l-clp.3.r,
      raw.low.5=.dustySplit(DNAStringSet(str_sub(seq.5, 1+clp.5.l, -1-clp.5.r)), 100000),
      raw.low.3=.dustySplit(DNAStringSet(str_sub(seq.3, 1+clp.3.l, -1-clp.3.r)), 100000),
      low.5=0,
      low.3=0
    )]
    setkey(tmp.tbl, len.5)
    tmp.tbl[max.low, low.5:=raw.low.5 / low]
    setkey(tmp.tbl, len.3);
    tmp.tbl[max.low, low.3:=raw.low.3 / low]
    low.tbl <- tmp.tbl[,.(qname, low.5, low.3)]
    setkeyv(low.tbl, "qname")
    jnc <- jnc[low.tbl]
    return(jnc)
}

.orientJunctions <- function(jnc, ann, debug=FALSE) {
    ## canonical orientation
    if (!ann$par$stranded) {
        mot.cmp <- jnc[,ifelse(type==1, FALSE,
                        ifelse(type==2, TRUE, NA))]
        
        ts <- c("D", "A", "E", "G", ".")
        score <- data.table(expand.grid(ts, ts))
        setnames(score, c("5", "3"))
                       # D   A   E   G   .
        score$val <- c(  1,  1,  1,  1,  1, # ->D
                        16,  1, 14, 12, 10, # ->A
                        15,  1,  9,  7,  5, # ->E
                        13,  1,  8,  4,  2, # ->G
                        11,  1,  6,  3,  0  # ->.
                       )
        setkeyv(score,  c("5", "3"))
        ## 1. split into 100k chunks to avoid NCList error, 2. do not creat temp. variables to save some memory?
        cls <- ann$features[ann$features$type %in% ts]
        typ <- data.table(
          cls.5.1 = unlist(unname(lapply(split(with(jnc,
            GRanges(chr.5, IRanges(pos.5-2, pos.5+2), str.5)), factor(ceiling(seq_len(nrow(jnc)) / 100000))),
            function(j) {cls[findOverlaps(j, cls, select="first")]$type}))),
          cls.3.1 = unlist(unname(lapply(split(with(jnc,
            GRanges(chr.3, IRanges(pos.3-2, pos.3+2), str.3)), factor(ceiling(seq_len(nrow(jnc)) / 100000))),
            function(j) {cls[findOverlaps(j, cls, select="first")]$type}))),
          cls.5.2 = unlist(unname(lapply(split(with(jnc,
                                                    GRanges(chr.3, IRanges(pos.3-2, pos.3+2), ifelse(str.3=="+", "-", "+"))),
                                               factor(ceiling(seq_len(nrow(jnc)) / 100000))),
            function(j) {cls[findOverlaps(j, cls, select="first")]$type}))),
          cls.3.2 = unlist(unname(lapply(split(with(jnc,
                                                    GRanges(chr.5, IRanges(pos.5-2, pos.5+2), ifelse(str.5=="+", "-", "+"))),
                                               factor(ceiling(seq_len(nrow(jnc)) / 100000))),
            function(j) {cls[findOverlaps(j, cls, select="first")]$type}))),
          idx = seq_len(nrow(jnc))
        )
        rm(cls); gc() ## save some memory
        setkey(typ, cls.5.1, cls.3.1)
        typ[score, val.1:=val]
        setkey(typ, cls.5.2, cls.3.2)
        typ[score, val.2:=val]
        ##
        jnc.cmp <- typ[order(idx), ifelse(val.1 > val.2, FALSE,
                                   ifelse(val.2 > val.1, TRUE , NA))]
        rm(typ); gc() ## save some memory
        ##
        chr.cmp <- jnc[,ifelse(chr.5 < chr.3, FALSE,
                        ifelse(chr.3 < chr.5, TRUE , NA))]
        pos.cmp <- jnc[,ifelse(pos.5 < pos.3, FALSE,
                        ifelse(pos.3 < pos.5, TRUE , NA))]
        str.cmp <- jnc[,ifelse(str.5 < str.3, FALSE,
                        ifelse(str.3 < str.5, TRUE , NA))]
        ## strand orientation logic: pick DA motif, pick highest score, pick S-L, else pick by rule
        fin.cmp <- ifelse(!is.na(mot.cmp), mot.cmp,
                   ifelse(!is.na(jnc.cmp), jnc.cmp,
                   ifelse(!is.na(chr.cmp), chr.cmp,
                   ifelse(!is.na(pos.cmp), pos.cmp,
                   ifelse(!is.na(str.cmp), str.cmp, FALSE)))))
        if (debug) {
            jnc$orient <- ifelse(!is.na(mot.cmp), ifelse(mot.cmp, "mot-rev", "mot-fwd"),
                          ifelse(!is.na(jnc.cmp), ifelse(jnc.cmp, "jnc-rev", "jnc-fwd"),
                          ifelse(!is.na(chr.cmp), ifelse(chr.cmp, "chr-rev", "chr-fwd"),
                          ifelse(!is.na(pos.cmp), ifelse(pos.cmp, "pos-rev", "pos-fwd"),
                          ifelse(!is.na(str.cmp), ifelse(str.cmp, "str-rev", "str-fwd"), "error")))))
        }
        rm(jnc.cmp, chr.cmp, str.cmp, pos.cmp); gc() ## save some memory
        jnc[fin.cmp,
            ":="(
               chr.5=chr.3, chr.3=chr.5, pos.5=pos.3, pos.3=pos.5,
               str.5=ifelse(str.3=="+", "-", "+"), str.3=ifelse(str.5=="+", "-", "+"),
               rep.5=rep.3, rep.3=rep.5, err.5=err.3, err.3=err.5, ovr.5=ovr.3, ovr.3=ovr.5,
               seg.5.pos=seg.3.pos, seg.3.pos=seg.5.pos, seg.5.cig=seg.3.cig, seg.3.cig=seg.5.cig,
               clp.5.l=clp.3.l, clp.3.l=clp.5.l, clp.5.r=clp.3.r, clp.3.r=clp.5.r,
               seg.5.as=seg.3.as, seg.3.as=seg.5.as, seg.5.nm=seg.3.nm, seg.3.nm=seg.5.nm,
               seg.5.len=seg.3.len, seg.3.len=seg.5.len, type=ifelse(type==2L, 1L, type),
               seq.5=seq.3, seq.3=seq.5, low.5=low.3, low.3=low.5
            )]
    }
    ## orient PE sequences
    jnc[ftype=="pe" & str.5=="-", seq.5:=as.character(reverseComplement(DNAStringSet(seq.5)))]
    jnc[ftype=="pe" & str.3=="+", seq.3:=as.character(reverseComplement(DNAStringSet(seq.3)))]
    ## orient SE sequences
    jnc[ftype=="se" & str.5=="-" & str.3=="+",":="(seq.5=seq.3, seq.3=seq.5)]    
    jnc[ftype=="se" & str.5=="-" & str.3=="-", seq.5:=as.character(reverseComplement(DNAStringSet(seq.5)))]
    jnc[ftype=="se" & str.5=="+" & str.3=="+", seq.3:=as.character(reverseComplement(DNAStringSet(seq.3)))]
    return(jnc)
}

.defineHQjunctions <- function(jnc, ann) {
    hq.seg=ann$par$jnc.ovr
    hq.ovr=ann$par$jnc.ovr.pct
    hq.err=ann$par$jnc.err
    hq.low=ann$par$jnc.low
    ## calculate junction uniqueness
    jnc[, unq.jnc:=NA_real_]
    jnc[, unq.jnc:=1/.N, by=c("type","chr.5","pos.5","str.5","chr.3","pos.3","str.3","seg.5.cig","seg.3.cig")]
    ## definie HQ encompassing fragments
    encs <- jnc[type < 0]
    encs <- encs[,
      hq.jnc:=(
          (
            (ovr.5 > (hq.ovr * seg.5.len)) & (ovr.3 > (hq.ovr * seg.3.len))
          ) & 
          (
            (err.5 < hq.err) & (err.3 < hq.err) 
          ) &
          (
            (low.5 < hq.low) & (low.3 < hq.low)
          )
    )]
    ## definie HQ spanning fragments
    spns <- jnc[type > -1]
    spns <- spns[,
      hq.jnc:=(
        (
          ((ftype == "pe") & ((ovr.5 + ovr.3) > (hq.ovr * (seg.5.len + seg.3.len)))) |
          ((ftype == "se") & ((ovr.5 + ovr.3) > (hq.ovr * (seg.5.len            ))))
        ) &
        (
          (ovr.5 > hq.seg) & (ovr.3 > hq.seg)
        ) &
        (  
          (err.5 < hq.err) & (err.3 < hq.err)
        ) &
        (
          (low.5 < hq.low) & (low.3 < hq.low)
        )
    )]
    ## merge back
    jnc.hq <- rbind(encs, spns)
    return(jnc.hq)
}

#' @export
readJunctions <- function(dir, ann, debug=FALSE) {
    jnc <- .mergeJNC(dir, ann)
    jnc <- .dustJunctions(jnc)
    jnc <- .orientJunctions(jnc, ann, debug)
    jnc <- .defineHQjunctions(jnc, ann)
    jnc$sid <- factor(dir$sid)
    setkeyv(jnc, JNC.KEY)
    return(jnc)
}
