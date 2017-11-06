#' @import Biostrings
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom parallel mclapply
#' @importFrom stringr str_sub str_split str_match str_replace str_replace_all
#' @importFrom igraph graph.data.frame V clusters
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels keepStandardChromosomes Seqinfo seqlengths
#' @importFrom IRanges IRanges %over%
#' @importFrom ShortRead dustyScore
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom rtracklayer import import.bw
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BSgenome getSeq
#' @importFrom GenomicAlignments cigarWidthAlongQuerySpace cigarWidthAlongReferenceSpace
#' @importFrom data.table data.table fread setnames setkey setkeyv copy as.data.table uniqueN rbindlist key melt dcast.data.table
NULL

MAXINT <- 999999999L
JNC.KEY <- c("sid", "qname")
BPT.KEY <- c("sid", "type", "chr.5", "pos.5", "str.5", "chr.3", "pos.3", "str.3")
CHM.KEY <- c("sid", "locus_id.5.1", "locus_id.3.1", "locus_id.5.2", "locus_id.3.2")

COLS <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", 
    "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#FDC086", "#E6AB02", "#A6761D", "#7FC97F", "#BEAED4", "#386CB0",
    "#F0027F", "#BF5B17", "#666666", "#999999", "#FFFF33", "#FFFF99"
)
