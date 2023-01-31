#' Testing presence of allelic imbalance using BAF estimates of SNPs
#'
#' This function tests whether for a given segment and region a BAF seperation of SNPs can be detected
#' based on previous characterization of SNPs as either allele A or B.
#'
#' @param phasingTable A table of SNPs with their positions along with their characterization.
#'     Additionally, their BAF and the segment they are assigned to are also included.
#' @param newSegTable A table of segments for all regions with chromosome, start, end, width.
#'     Additionally, the segment ID as well as the BAF of the segment 
#'     and major and minor copy number are included.
#' @param seg.totalSNPs An integer used to filter minimum size of segments to test. 
#'     Defaults to 8 meaning segments with less than 8 SNPs are not tested.
#' @param seg.nSNPs An integer used to filter minimum size of SNPs per allele.
#'     Defaults to 4 meaning at least 4 SNPs need to be characterized per allele to test segment. 
#' @param testMSAI A flag whether called or rescued segments should be tested for MSAI.
#'     Defaults to TRUE meaning MSAI is tested in all called or rescued segments.
#' @return A table similar to the phasingTable provided as input, but also including 
#'     the characterization of whether allelic imbalance is detected/rescued in the segment.
#' @export
testingSCNA <- function(phasingTable, newSegTable, seg.totalSNPs = 8, seg.nSNPs = 4, testMSAI = TRUE) {
    testingTable <- phasingTable %>% dplyr::filter(class != "mixed")

    foo <- testingTable %>% 
        dplyr::group_by(region, segment, class) %>%
        dplyr::mutate(nSNPs = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(region, segment) %>%
        dplyr::mutate(totalSNPs = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::select(region, segment, class, nSNPs, totalSNPs) %>%
        unique()
    segmentsToRemove <- foo %>% 
        dplyr::filter(nSNPs < seg.nSNPs | totalSNPs < seg.totalSNPs | nSNPs == totalSNPs) %>% 
        dplyr::pull(segment) %>% 
        unique()

    testingTable <- testingTable %>% dplyr::filter(!segment %in% segmentsToRemove)
    if (nrow(testingTable) == 0) return(NULL)

    segmentsOfInterest <- unique(testingTable$segment)
    samplesOfInterest <- unique(testingTable$region)

    combinationsToTest <- newSegTable %>% 
        dplyr::filter(BAF_segmented == 0.5) %>% 
        dplyr::select(region, segment) %>% 
        dplyr::filter(segment %in% segmentsOfInterest, region %in% samplesOfInterest)

    combinationsToTest$pval <- sapply(1:nrow(combinationsToTest), function(x) {
        tmpTestingTable <- testingTable %>% dplyr::filter(segment %in% combinationsToTest$segment[x], region %in% combinationsToTest$region[x])
        tmp <- wilcox.test(tmpTestingTable %>% dplyr::filter(class == "top") %>% dplyr::pull(BAF), tmpTestingTable %>% dplyr::filter(class == "bottom") %>% dplyr::pull(BAF))
        return(tmp$p.value)
    })

    newSegTable <- newSegTable %>% 
        dplyr::left_join(combinationsToTest, by = c("region", "segment")) %>% 
        dplyr::mutate(type = ifelse(BAF_segmented != 0.5, "called", ifelse(BAF_segmented == 0.5, ifelse(pval < 0.005, "SCNA rescued", "not detected"), NA)))
    newPhasingTable <- phasingTable %>% 
        dplyr::left_join(newSegTable %>% dplyr::select(region, segment, type), by = c("region", "segment")) %>%
        dplyr::mutate(type = ifelse(class == "mixed", "not called", type))

    if (testMSAI) {
        possibleMSAItoTest <- newPhasingTable %>% 
            dplyr::filter(type %in% c("called", "SCNA rescued")) %>% 
            dplyr::group_by(region, segment, class) %>% 
            dplyr::summarize(n = dplyr::n()) %>% 
            dplyr::filter(n > 8) %>%
            dplyr::group_by(region, segment) %>%
            dplyr::tally() %>%
            dplyr::filter(n == 2) %>% 
            dplyr::ungroup() %>%
            dplyr::select(region, segment) %>% 
            unique()

        possibleMSAItoTest$pval <- sapply(1:nrow(possibleMSAItoTest), function(x) {
            tmpTestingTable <- testingTable %>% dplyr::filter(segment %in% possibleMSAItoTest$segment[x], region %in% possibleMSAItoTest$region[x])
            tmp <- wilcox.test(tmpTestingTable %>% dplyr::filter(class == "top") %>% dplyr::pull(BAF), tmpTestingTable %>% dplyr::filter(class == "bottom") %>% dplyr::pull(BAF), alternative = "less")
            return(tmp$p.value)
        })
        possibleMSAItoTest <- possibleMSAItoTest %>%
            dplyr::mutate(MSAI = ifelse(pval < 0.01, "yes", "no")) %>%
            dplyr::select(-pval)
    
        newPhasingTable <- newPhasingTable %>% 
            dplyr::left_join(possibleMSAItoTest, by = c("region", "segment")) %>%
            dplyr::mutate(type2 = type) %>%
            dplyr::mutate(type = ifelse(MSAI %in% "yes", "MSAI", type2)) 
    }
    return(newPhasingTable)
}