#' Phasing SNPs on segments across samples
#'
#' This function phases SNPs located on the same segment across samples
#'
#' @param snpTable A table of SNPs with their positions.
#' @param segTable A table of segments to consider including chromosome start and stop positions as well as a BAF value per segment.
#' @param samples A vector of sample names to include for phasing. 
#'     Defaults to NULL in which case all samples are considered.
#' @return A table similar to the SNPtable provided as input, but also including 
#'     the segment the SNPs can be found on as well as classification whether they
#'     are likely part of allele A or B.
#' @export
phaseSegmentsWithSegs <- function(snpTable, segTable, samples = NULL) {
    ### if no samples specified, use all samples for phasing
    if (is.null(samples)) samples <- unique(segTable$region)
    segTable <- segTable %>%
        dplyr::filter(region %in% samples)
    ### create minimum consistent regions data frame
    ### remove segments that don't have a BAF_segmented value to minimize small segments
    ### TODO check if this makes sense or whether we need to use all segments
    segTable      <- segTable %>% dplyr::filter(!is.na(BAF_segmented))
    segTable.gr   <- GenomicRanges::makeGRangesFromDataFrame(segTable, start.field = "startpos", end.field = "endpos")
    minConsReg.gr <- disjoin(segTable.gr)
    minConsReg    <- as.data.frame(minConsReg.gr)
    minConsReg$segment <- seq(1, nrow(minConsReg), by = 1)

    ### allocate which minimum consistent segment snps belong to
    snpTable.gr <- GenomicRanges::makeGRangesFromDataFrame(snpTable, start.field = "pos", end.field = "pos")
    overlap <- GenomicRanges::findOverlaps(minConsReg.gr, snpTable.gr)
    snpTable$segment <- NA
    snpTable$segment[S4Vectors::subjectHits(overlap)] <- minConsReg$segment[S4Vectors::queryHits(overlap)]

    newSegTable <- lapply(samples, function(samp) {
        newSegTableRegion <- cbind(region = rep(samp, nrow(minConsReg)), minConsReg)
        ### allocate which minimum consistent segment snps belong to
        newSegTableRegion.gr <- GenomicRanges::makeGRangesFromDataFrame(newSegTableRegion, start.field = "start", end.field = "end")
        segTableRegion <- segTable %>% dplyr::filter(region == samp)
        if (nrow(segTableRegion) == 0) {
            newSegTableRegion$BAF_segmented <- 0.5
            newSegTableRegion$nMajor <- 1
            newSegTableRegion$nMinor <- 1
            return(newSegTableRegion)
        } else {
            segTableRegion.gr    <- GenomicRanges::makeGRangesFromDataFrame(segTableRegion, start.field = "startpos", end.field = "endpos")
            overlap <- GenomicRanges::findOverlaps(segTableRegion.gr, newSegTableRegion.gr)
            newSegTableRegion$BAF_segmented <- NA
            newSegTableRegion$BAF_segmented[S4Vectors::subjectHits(overlap)] <- segTableRegion$BAF_segmented[S4Vectors::queryHits(overlap)]
            newSegTableRegion$nMajor <- NA
            newSegTableRegion$nMajor[S4Vectors::subjectHits(overlap)] <- segTableRegion$nMajor[S4Vectors::queryHits(overlap)]
            newSegTableRegion$nMinor <- NA
            newSegTableRegion$nMinor[S4Vectors::subjectHits(overlap)] <- segTableRegion$nMinor[S4Vectors::queryHits(overlap)]
            return(newSegTableRegion)
        }
    })
    newSegTable <- Reduce(rbind, newSegTable)

    phasingTable <- newSegTable %>% 
        dplyr::mutate(BAFdiff = abs(0.5 - BAF_segmented)) %>% 
        dplyr::group_by(segment) %>% 
        dplyr::filter(BAFdiff == max(BAFdiff, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::select(region, segment, BAFdiff)

    tmp <- phasingTable %>% 
        dplyr::left_join(snpTable, by = c("segment", "region")) %>%
        dplyr::filter(!is.na(BAF)) %>%
        dplyr::mutate(class = ifelse(BAFdiff == 0, "mixed", ifelse(BAF < 0.5, "bottom", "top"))) %>%
        dplyr::select(chrom, pos, ref, alt, class) %>%
        unique() %>%
        dplyr::left_join(snpTable, by = c("chrom", "pos", "ref", "alt"))
    return(list(phasingTable = tmp, newSegTable = newSegTable))
}
