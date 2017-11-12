.makeFeatures <- function(par, seqi) {
    txdb <- suppressMessages(makeTxDbFromGFF(par$gtf.fn))
    seqlevels(txdb, pruning.mode="coarse") <- seqlevels(seqi)
    gme <- GRanges(seqnames(seqi), IRanges(start=1, end=seqlengths(seqi)), strand="*")
    names(gme) <- seqnames(seqi)
    gme$type <- "."
    gen <- genes(txdb)
    mcols(gen) <- NULL
    gen$type <- "G"
    txs <- transcripts(txdb)
    names(txs) <- mcols(txs)$tx_name
    mcols(txs) <- NULL
    txs$type <- "T"
    ins <- unlist(intronsByTranscript(txdb, use.names=TRUE))
    mcols(ins) <- NULL
    ins$type <- "I"
    exs <- unlist(exonsBy(txdb, use.names=TRUE))
    mcols(exs) <- NULL
    exs$type <- "E"
    u5s <- unlist(fiveUTRsByTranscript(txdb, use.names=TRUE))
    mcols(u5s) <- NULL
    u5s$type <- "5"
    u3s <- unlist(threeUTRsByTranscript(txdb, use.names=TRUE))
    mcols(u3s) <- NULL
    u3s$type <- "3"
    cds <- unlist(cdsBy(txdb, use.names=TRUE))
    mcols(cds) <- NULL
    cds$type <- "C"
    ssl <- ins
    end(ssl) <- start(ssl)
    ssl$type <- ifelse(strand(ssl)=="+", "D", "A")
    ssr <- ins
    start(ssr) <- end(ssr)
    ssr$type <- ifelse(strand(ssr)=="+", "A", "D")
    sss <- c(ssl, ssr)
    fts <- c(sss, cds, u3s, u5s, exs, ins, txs, gen, gme)
    fts$split_id <- names(fts)
    names(fts) <- NULL
    seqinfo(fts) <- seqi
    return(fts)
}

.makeGenes <- function(par, seqi) {
    tmp <- fread(paste(system2("grep", sprintf("'\tgene\t' %s", par$gtf.fn), stdout=TRUE),
                       collapse="\n"), showProgress=FALSE, header=FALSE)
    gen <- with(tmp, GRanges(V1, IRanges(V4, V5), V7,
                      gene_id = fread(paste(str_match(V9, "gene_id[^;]+"),
                        collapse="\n"), showProgress=FALSE, header=FALSE)$V2,
                      gene_name = fread(paste(str_match(V9, "gene_name[^;]+"),
                        collapse="\n"), showProgress=FALSE, header=FALSE)$V2,
                      gene_type = fread(paste(str_match(V9, "gene_type[^;]+"),
                        collapse="\n"), showProgress=FALSE, header=FALSE)$V2
                      ))
    seqlevels(gen, pruning.mode="coarse") <- seqlevels(seqi)
    seqinfo(gen) <- seqi
    gene.names <- fread(par$goi.fn, head=FALSE)$V1
    mcols(gen)$goi <- FALSE
    mcols(gen[gen$gene_name %in% gene.names])$goi <- TRUE
    mcols(gen)$protein <- (gen$gene_type == "protein_coding")
    names(gen) <- gen$gene_id
    return(gen)
}

.makeTranscripts <- function(par, fts, seqi) {
    tmp <- fread(paste(system2("grep", sprintf("'\ttranscript\t' %s", par$gtf.fn), stdout=TRUE),
                       collapse="\n"), showProgress=FALSE, header=FALSE)
    tx <- with(tmp, GRanges(V1, IRanges(V4, V5), V7,
                      gene_id = fread(paste(str_match(V9, "gene_id[^;]+"),
                        collapse="\n"), showProgress=FALSE, header=FALSE)$V2,
                      transcript_id = fread(paste(str_match(V9, "transcript_id[^;]+"),
                        collapse="\n"), showProgress=FALSE, header=FALSE)$V2,
                      tags = fread(paste(str_match(tmp$V9, "tags[^;]+"),
                        collapse="\n"), showProgress=FALSE, header=FALSE)$V2
                      ))
    seqlevels(tx, pruning.mode="coarse") <- seqlevels(seqi)
    seqinfo(tx) <- seqi
    names(tx) <- tx$transcript_id
    ebt <- split(fts[fts$type=="E"], fts[fts$type=="E"]$split_id)
    ebs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, ebt)
    tx$tx.seq <- unname(ebs[names(tx)])
    fcex <- fts[fts$type=="C"]
    rcex <- rev(fcex)
    cex1 <- fcex[!duplicated(fcex$split_id)]
    cexn <- rev(rcex[!duplicated(rcex$split_id)])
    strand(cexn) <- ifelse(strand(cexn)=="+", "-", "+")
    base1 <- promoters(cex1, 0, 1)
    basen <- promoters(cexn, 0, 1)
    strand(basen) <- ifelse(strand(basen)=="+", "-", "+")
    mcols(tx)$cds1 <- 1
    mcols(tx)$cdsN <- 0
    tx[base1$split_id]$cds1 <- start(pmapToTranscripts(base1, ebt[base1$split_id]))
    tx[basen$split_id]$cdsN <- start(pmapToTranscripts(basen, ebt[basen$split_id]))
    return(tx)
}

.makeFrames <- function(fts) {
    cds <- fts[fts$type=="C"]
    cds <- split(cds, cds$split_id)
    end.ovr <- cumsum(width(cds)) %% 3
    end.ovr <- unlist(end.ovr)
    beg.ovr <- (cumsum(width(cds)) - width(cds) + 1) %% 3
    beg.ovr <- unlist(beg.ovr)
    beg.ovr[beg.ovr==0] <- 3
    beg.ovr <- beg.ovr - 1
    cds.exons <- unlist(cds, use.names=FALSE)
    mcols(cds.exons) <- NULL
    start(cds.exons) <- start(cds.exons) - 1
    end(cds.exons) <- end(cds.exons) + 1
    cds.exons$a.frame <- unname(beg.ovr)
    cds.exons$d.frame <- unname(end.ovr)
    names(cds.exons) <- NULL
    frame <- cds.exons
    tmp.p <- frame[strand(frame)=="+"]
    start(tmp.p) <- end(tmp.p)
    tmp.m <- frame[strand(frame)=="-"]
    end(tmp.m) <- start(tmp.m)
    frame.d <- c(tmp.p, tmp.m)
    frame.d$frame <- frame.d$d.frame
    frame.d$type <- "D"
    frame.d$a.frame <- NULL
    frame.d$d.frame <- NULL
    tmp.p <- frame[strand(frame)=="+"]
    end(tmp.p) <- start(tmp.p)
    tmp.m <- frame[strand(frame)=="-"]
    start(tmp.m) <- end(tmp.m)
    frame.a <- c(tmp.p, tmp.m)
    frame.a$frame <- frame.a$a.frame
    frame.a$type <- "A"
    frame.a$a.frame <- NULL
    frame.a$d.frame <- NULL
    frame.da <- c(frame.d, frame.a)
    frame.da <- frame.da[!duplicated(frame.da)]
    return(list(frame=frame, frame.da=frame.da))
}

.fixRanges <- function(rng, seqi, stranded=FALSE) {
    if (!stranded) {
        strand(rng) <- "*"
    }
    seqlevels(rng, pruning.mode="coarse") <- seqlevels(seqi)
    seqinfo(rng) <- seqi
    fix <- sort(rng)
    return(fix)
}

.parseSVs <- function(par, core) {
    tmp.sv <- fread(paste0(system2("zcat", par$sv.fn, stdout=TRUE), collapse="\n"),
                    showProgress=FALSE, skip="#CHROM", colClasses=list(character=1))
    tmp.sv <- data.table(
      chr = paste0("chr", tmp.sv$`#CHROM`),
      beg = tmp.sv$POS,
      typ = str_match(tmp.sv$INFO, "SVTYPE=([^;,]*)")[,2],
      end = as.integer(str_match(tmp.sv$INFO, ";END=([^;,]*)")[,2]),
      af = as.numeric(str_match(tmp.sv$INFO, ";AF=([^;,]*)")[,2]),
      cs = str_match(tmp.sv$INFO, ";CS=([^;,]*)")[,2],
      info=tmp.sv$INFO
    )
    tmp.sv[is.na(end), end:=beg]
    svs <- with(tmp.sv, GRanges(chr, IRanges(beg, end), "*", type=typ, meth=cs, af=af))
    svs <- .fixRanges(svs, core$seqi, stranded=FALSE)
    return(svs)
}

.makeCoreAnnotations <- function(par) {
    ##### seqi,features,genes,transcripts,loci,frames,cytobands
    ## SEQI
    seqi <- keepStandardChromosomes(Seqinfo(genome=par$gme))
    ## FEATURES
    fts <- .makeFeatures(par, seqi)
    ## GENES
    genes <- .makeGenes(par, seqi)
    transcripts <- .makeTranscripts(par, fts, seqi)
    ## FRAMES
    frames <- .makeFrames(fts)
    ## LOCI
    tmp <- reduce(genes)
    loci <- c(tmp, gaps(tmp))
    loci <- loci[strand(loci) != "*"]
    mcols(loci)$locus_id <- sprintf("L%06d", 1:length(loci))
    mcols(loci)$genic <- FALSE
    mcols(loci[unique(sort(subjectHits(findOverlaps(genes, loci))))])$genic <- TRUE
    ## CYTOBANDS
    tmp = fread(par$cyt.fn)[V1 %in% seqlevels(loci)[1:24]]
    cyt <- with(tmp, GRanges(V1, IRanges(V2+1, V3), "*", cytoband_id=paste(V1, V4, sep="_"),
                             arm=str_sub(V4, 1, 1),
                             band=str_sub(V4, 2), stain=V5))
    seqlevels(cyt, pruning.mode="coarse") <- seqlevels(seqi)
    seqinfo(cyt) <- seqi
    core <- list(seqi=seqi, features=fts, genes=genes, transcripts=transcripts, loci=loci, frames=frames, cytobands=cyt)
    return(core)
}

.makeOverlapAnnotations <- function(core) {
    ## gene-locus overlaps
    tmp <- findOverlaps(core$genes, core$loci)
    grp.map <- data.table(gene_id=core$genes[queryHits(tmp)]$gene_id, gene_locus=subjectHits(tmp))
    go <- data.table(locus_id=mcols(core$loci)$locus_id[grp.map$gene_locus],
                     gene_id=grp.map$gene_id)
    ## cytoband-locus overlaps
    tmp <- data.table(as.data.frame(findOverlaps(core$cytobands, core$loci)))
    tmp$overlap <- width(pintersect(core$cytobands[tmp$queryHits], core$loci[tmp$subjectHits]))
    tmp <- tmp[order(-overlap)]
    tmp <- tmp[!duplicated(tmp$subjectHits)]
    co <- data.table(locus_id=factor(core$loci[tmp$subjectHits]$locus_id, levels=core$loci$locus_id),
                     cytoband_id=core$cyt[tmp$queryHits]$cytoband_id
                     )
    ## locus-locus overlaps
    tmp <- findOverlaps(core$loci, core$loci, ignore.strand=TRUE)
    tmp <- tmp[queryHits(tmp) != subjectHits(tmp)]
    lo <- data.table(locus_id.5=factor(mcols(core$loci)$locus_id[queryHits(tmp)], levels=core$loci$locus_id),
                     locus_id.3=factor(mcols(core$loci)$locus_id[subjectHits(tmp)], levels=core$loci$locus_id))
    ## locus-locus adjacency
    tmp <- findOverlaps(core$loci, drop.self=TRUE, maxgap=0)
    la <- data.table(locus_id.5=factor(mcols(core$loci)$locus_id[queryHits(tmp)], levels=core$loci$locus_id),
                     locus_id.3=factor(mcols(core$loci)$locus_id[subjectHits(tmp)], levels=core$loci$locus_id))
    ## locus-to-index
    l2i <- data.table(locus_id=core$loci$locus_id, idx=1:length(core$loci))
    ## gene-index
    g2i <- data.table(gene_id=core$genes$gene_id, idx=1:length(core$genes))
    ## cytoband-index
    c2i <- data.table(cytoband_id=core$cytobands$cytoband_id, idx=1:length(core$cytobands))
    ## ann object
    setkeyv(lo, c("locus_id.5", "locus_id.3"))
    setkeyv(la, c("locus_id.5", "locus_id.3"))
    setkey(go, "locus_id")
    setkey(co, "locus_id")
    setkey(l2i, "locus_id")
    setkey(g2i, "gene_id")
    setkey(c2i, "cytoband_id")
    over <- list(gene2idx = g2i,
                 locus2idx = l2i,
                 cytoband2idx = c2i,
                 gene.ovr = go,
                 cytoband.ovr = co,
                 locus.ovr = lo,
                 locus.adj = la)
    return(over)
}

.makeMaskAnnotations <- function(par, core) {
    ## retrotransposed
    rets <- .fixRanges(import(par$ret.fn), core$seqi)
    ## SNORNAs
    snos <- .fixRanges(import(par$sno.fn), core$seqi)
    ## repeats
    reps <- .fixRanges(import(par$rep.fn), core$seqi)
    ## SVs
    svs <- .parseSVs(par, core)
    ## altlocs
    alts <- with(fread(par$alt.fn)[,.(chr=paste0("chr", parent_name), beg=parent_start, end=parent_stop)],
                 GRanges(chr, IRanges(beg, end), strand="*"))
    alts <- .fixRanges(alts, core$seqi)
    ## Low-comlexity
    tmp.low <- fread(paste0(system2("zcat", par$low.fn, stdout=TRUE), collapse="\n"), showProgress=FALSE, header=FALSE)
    tmp.low <- with(tmp.low, GRanges(V1, IRanges(V2+1, V3), "*"))
    low <- .fixRanges(tmp.low, core$seqi)
    ## IGX
    igx <- .fixRanges(with(fread(par$igx.fn, header=FALSE), GRanges(V1, IRanges(V2, V3))), core$seqi)
    mask <- list(repeats=reps, retros=rets, snornas=snos, altlocs=alts, svs=svs, igx=igx, low=low)
    return(mask)
}

.makeMiscAnnotations <- function(par, core) {
    ## known artifacts
    tmp <- fread(paste0(system2("zcat", par$art.fn, stdout=TRUE), collapse="\n"),
                 showProgress=FALSE, header=TRUE)
    arts <- tmp[(chr.5 %in% seqlevels(core$seqi)) & (chr.3 %in% seqlevels(core$seqi))]
    ## segments
    tmp <- fread(paste0(system2("zcat", par$seg.fn, stdout=TRUE), collapse="\n"),
                 showProgress=FALSE, header=TRUE)[,c(2:4,8:10),with=FALSE]
    segs <- tmp[(chrom %in% seqlevels(core$seqi)) & (otherChrom %in% seqlevels(core$seqi))]
    setnames(segs, c("chr.5", "start.5", "end.5", "chr.3", "start.3", "end.3"))
    ## genes to blacklist
    gtb <- core$genes[core$genes$gene_name %in% readLines(par$gtb.fn)]
    ## loci of interest
    tmp <- with(fread(par$loi.fn, showProgress=FALSE, header=FALSE), GRanges(V1, IRanges(V2, V3), name=V4))
    loi <- .fixRanges(tmp, core$seqi)
    ## misc
    misc <- list(
        artifacts = arts,
        segdups = segs,
        gtb = gtb,
        loi = loi
    )
    return(misc)
}

#' @export
makeAnnotations <- function(par) {
    core <- .makeCoreAnnotations(par)
    over <- .makeOverlapAnnotations(core)
    mask <- .makeMaskAnnotations(par, core)
    misc <- .makeMiscAnnotations(par, core)
    ann <- c(core, over, mask, misc, list(par=par))
    return(ann)
}
