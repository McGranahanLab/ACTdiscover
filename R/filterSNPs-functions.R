    filterProblematicSNPs <- function(snpTable, problematicLoci = NULL) {
        if (is.null(problematicLoci)) {
            problematicLoci <- ACTdiscover::problematicLoci
        }
        problematicLoci$chr <- paste0("chr", problematicLoci$chr)
        names(problematicLoci)[1] <- "chrom"
        problematicLoci$problematic <- "yes"

        snpTable <- snpTable %>% 
            dplyr::left_join(problematicLoci, by = c("chrom", "pos")) %>%
            dplyr::filter(is.na(problematic)) %>%
            dplyr::select(-problematic)
        return(snpTable)
    }

    filterBlacklistRegions <- function(snpTable, blacklistRegions = NULL) {
        if (is.null(blacklistRegions)) {
            blacklistRegions <- ACTdiscover::blacklistRegions
        }
        snpTable.gr <- GenomicRanges::makeGRangesFromDataFrame(snpTable, start.field = "pos", end.field = "pos")
        
        overlap <- GenomicRanges::findOverlaps(blacklistRegions, snpTable.gr)

        snpTable <- snpTable[-S4Vectors::subjectHits(overlap),]
        return(snpTable)
    }