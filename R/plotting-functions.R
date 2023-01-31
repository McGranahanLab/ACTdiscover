#' Plotting phasing of SNPs across the genome for multiple samples
#'
#' This function plots the phased SNPs for all samples of interest.
#'
#' @param snpTable A table of SNPs with their positions along with their characterization.
#' @param genome A genome version to use, for example BSgenome.Hsapiens.UCSC.hg19.
#' @param samples A vector of sample names to include for plotting in the order they appear in. 
#'     Defaults to NULL in which case all samples are considered.    
#' @param col.palette A named vector of colors for plotting the SNPs.
#'     Defaults to orange and purple for allele A and B, respectively.
#' @return A ggplot object that can be plotted directly or saved as a pdf.
#' @export
plotPhasedSNPs <- function(snpTable, genome, samples = NULL, col.palette = setNames(c("#f1a340", "#998ec3", "#bdbdbd"), c("top", "bottom", "mixed"))) { 
    if (is.null(samples)) samples <- unique(unique(snpTable$region))

    snpTable.plot <- snpTable %>% 
        dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
        dplyr::filter(chrom %in% 1:22) %>%
        dplyr::filter(region %in% samples) %>%
        dplyr::mutate(region = factor(region, levels = samples))

    ### getting positions of tick marks for chromosome labels
    chr.sizes       <- seqlengths(genome)[paste0('chr', c(1:22, 'X', 'Y'))]
    cumsum.chr.size <- cumsum(as.numeric(chr.sizes))
    chr_tick_marks <- as.numeric(chr.sizes) / 2
    chr_tick_marks[2:length(chr_tick_marks)] <- chr_tick_marks[2:length(chr_tick_marks)] + cumsum.chr.size[1:(length(cumsum.chr.size) - 1)]

    ### adjust position to continuous variable instead of resetting at 0 every chromosome
    snpTable.plot$adjPos <- snpTable.plot$pos
    for (i in 2:24) {
        snpTable.plot$adjPos[snpTable.plot$chrom == i] <- snpTable.plot$adjPos[snpTable.plot$chrom == i] + cumsum.chr.size[i - 1]
    }
    snpTable.plot$adjPos <- as.numeric(as.character(snpTable.plot$adjPos))

    p <- ggplot2::ggplot(snpTable.plot, ggplot2::aes(adjPos, BAF, color = class)) + 
            ggplot2::geom_point() + 
            ggplot2::geom_vline(xintercept = cumsum.chr.size[1:22], linetype = 2) +
            ggplot2::geom_hline(yintercept = 0.5, linetype = 2, color = "gray") +
            ggplot2::xlab("") + 
            ggplot2::scale_x_continuous(breaks = chr_tick_marks[1:22], labels = paste0('chr', c(1:22)), expand = c(0, 0)) +
            ggplot2::scale_color_manual(values = col.palette, guide = "none") +
            ggplot2::facet_wrap(~ region, ncol = 1) +
            cowplot::theme_cowplot(font_size = 12) +
            ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    return(p)
}

#' Plotting phasing of SNPs across the genome for multiple samples with track per sample
#'
#' This function plots the phased SNPs for all samples of interest.
#' Additionally, for each sample a track is added at the bottom of the plot
#' highlighting the classification of the segment within that sample.
#' Standard classification includes not called, not detected, called, SCNA rescued, MSAI.
#'
#' @param snpTable A table of SNPs with their positions along with their characterization.
#' @param genome A genome version to use, for example BSgenome.Hsapiens.UCSC.hg19.
#' @param samples A vector of sample names to include for plotting in the order they appear in. 
#'     Defaults to NULL in which case all samples are considered.    
#' @param col.palette A named vector of colors for plotting the SNPs.
#'     Defaults to orange and purple for allele A and B, respectively.
#' @param fill.palette A named vector of colors for plotting the characterization of each segment.
#'     Defaults to orange and purple for allele A and B, respectively.
#' @return A ggplot object that can be plotted directly or saved as a pdf.
#' @export
plotPhasedSNPsWithTrack <- function(snpTable, genome, samples = NULL, col.palette = setNames(c("#f1a340", "#998ec3", "#bdbdbd"), c("top", "bottom", "mixed")), fill.palette = setNames(c("#D95980", "#335120", "#9DCC5F", "#63AAC0", "#bdbdbd"), c("MSAI", "SCNA rescued", "called", "not detected", "not called"))) { 
    if (is.null(samples)) samples <- unique(unique(snpTable$region))

    snpTable.plot <- snpTable %>% 
        dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
        dplyr::filter(chrom %in% 1:22) %>%
        dplyr::filter(region %in% samples) %>%
        dplyr::mutate(region = factor(region, levels = samples)) %>%
        dplyr::mutate(type = factor(type, levels = rev(c("MSAI", "SCNA rescued", "called", "not detected", "not called", "NA"))))

    ### getting positions of tick marks for chromosome labels
    chr.sizes       <- seqlengths(genome)[paste0('chr', c(1:22, 'X', 'Y'))]
    cumsum.chr.size <- cumsum(as.numeric(chr.sizes))
    chr_tick_marks <- as.numeric(chr.sizes) / 2
    chr_tick_marks[2:length(chr_tick_marks)] <- chr_tick_marks[2:length(chr_tick_marks)] + cumsum.chr.size[1:(length(cumsum.chr.size) - 1)]

    ### adjust position to continuous variable instead of resetting at 0 every chromosome
    snpTable.plot$adjPos <- snpTable.plot$pos
    for (i in 2:24) {
        snpTable.plot$adjPos[snpTable.plot$chrom == i] <- snpTable.plot$adjPos[snpTable.plot$chrom == i] + cumsum.chr.size[i - 1]
    }
    snpTable.plot$adjPos <- as.numeric(as.character(snpTable.plot$adjPos))
    snpTable.plot <- snpTable.plot %>% 
        dplyr::group_by(region, segment) %>% 
        dplyr::mutate(startpos = min(adjPos), endpos = max(adjPos))

    p <- ggplot2::ggplot(snpTable.plot) + 
            ggplot2::geom_point(ggplot2::aes(adjPos, BAF, color = class)) + 
            ggplot2::geom_vline(xintercept = cumsum.chr.size[1:22], linetype = 2) +
            ggplot2::geom_hline(yintercept = 0.5, linetype = 2, color = "gray") +
            ggplot2::geom_rect(data = snpTable.plot %>% dplyr::select(region, segment, startpos, endpos, type) %>% unique(), ggplot2::aes(xmin = startpos, xmax = endpos, ymin = -0.125, ymax = -0.025, fill = type)) +
            ggplot2::xlab("") + 
            ggplot2::scale_x_continuous(breaks = chr_tick_marks[1:22], labels = paste0('chr', c(1:22)), expand = c(0, 0)) +
            ggplot2::scale_color_manual(values = col.palette, guide = "none") +
            ggplot2::scale_fill_manual(name = "Sample level classification",
                              values = rev(fill.palette), 
                              labels = rev(c("Mirrored subclonal allelic imbalance", "Allelic imbalance rescued", "Allelic imbalance in single region", "Allelic imbalance not detected", "Haplotype information unavailable"))) +
            ggplot2::facet_wrap(~ region, ncol = 1) +
            cowplot::theme_cowplot(font_size = 12) +
            ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    return(p)
}

#' Plotting phasing of SNPs across the genome for multiple samples with track per sample and summary track
#'
#' This function plots the phased SNPs for all samples of interest.
#' Additionally, for each sample a track is added at the bottom of the plot
#' highlighting the classification of the segment within that sample.
#' Standard classification includes not called, not detected, called, SCNA rescued, MSAI.
#' At the bottom of the plot a summary track is included highlighting whether each segment
#' has an SCNA detected in any region, whether this is homogeneous and detected or rescued in all regions,
#' whether it is heterogeneous and only detected in a subset of regions, or whether there is MSAI.
#'
#' @param newPhasingTable A table of SNPs with their positions along with their characterization.
#' @param genome A genome version to use, for example BSgenome.Hsapiens.UCSC.hg19.
#' @param samples A vector of sample names to include for plotting in the order they appear in. 
#'     Defaults to NULL in which case all samples are considered.    
#' @param col.palette A named vector of colors for plotting the SNPs.
#'     Defaults to orange and purple for allele A and B, respectively.
#' @param fill.palette A named vector of colors for plotting the characterization of each segment.
#'     Defaults to orange and purple for allele A and B, respectively.
#' @return A ggplot object that can be plotted directly or saved as a pdf.
#' @export
plotPhasedSNPsWithTrackAndSummary <- function(newPhasingTable, genome, samples = NULL, col.palette = setNames(c("#f1a340", "#998ec3", "#bdbdbd"), c("top", "bottom", "mixed")), fill.palette = setNames(c("#D95980", "#335120", "#9DCC5F", "#63AAC0", "#bdbdbd"), c("MSAI", "SCNA rescued", "called", "not detected", "not called"))) { 
    if (is.null(samples)) samples <- unique(unique(newPhasingTable$region))

    if (length(grep("simulations", names(newPhasingTable))) == 0) {
        newPhasingTable$simulations <- "not tested"
    }

    snpTable.plot <- newPhasingTable %>% 
        dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
        dplyr::filter(chrom %in% 1:22) %>%
        dplyr::filter(region %in% samples) %>%
        dplyr::mutate(region = factor(region, levels = samples)) %>%
        dplyr::mutate(newType = ifelse(type == "not detected", ifelse(simulations == "not detected", "unsure", "not detected"), type)) %>%
        dplyr::mutate(type = factor(type, levels = rev(c("MSAI", "SCNA rescued", "called", "not detected", "not called", "NA"))))

    ### getting positions of tick marks for chromosome labels
    chr.sizes       <- seqlengths(genome)[paste0('chr', c(1:22, 'X', 'Y'))]
    cumsum.chr.size <- cumsum(as.numeric(chr.sizes))
    chr_tick_marks <- as.numeric(chr.sizes) / 2
    chr_tick_marks[2:length(chr_tick_marks)] <- chr_tick_marks[2:length(chr_tick_marks)] + cumsum.chr.size[1:(length(cumsum.chr.size) - 1)]

    ### adjust position to continuous variable instead of resetting at 0 every chromosome
    snpTable.plot$adjPos <- snpTable.plot$pos
    for (i in 2:24) {
        snpTable.plot$adjPos[snpTable.plot$chrom == i] <- snpTable.plot$adjPos[snpTable.plot$chrom == i] + cumsum.chr.size[i - 1]
    }
    snpTable.plot$adjPos <- as.numeric(as.character(snpTable.plot$adjPos))
    snpTable.plot <- snpTable.plot %>% 
        dplyr::group_by(region, segment) %>% 
        dplyr::mutate(startpos = min(adjPos), endpos = max(adjPos))

    segSummary <- newPhasingTable %>% 
        dplyr::filter(region %in% samples) %>% 
        dplyr::select(region, segment, type, simulations) %>% 
        unique() %>% 
        dplyr::mutate(newType = ifelse(type == "not detected", ifelse(simulations == "not detected", "unsure", "not detected"), type)) %>%
        dplyr::filter(newType != "unsure") %>%
        dplyr::group_by(segment) %>% 
        dplyr::summarize(characterization = ifelse(any(type == "MSAI"), 
                                                   "MSAI", 
                                                   ifelse(all(type %in% "not called"), 
                                                          "no SCNA", 
                                                          ifelse(all(type %in% c("called", "SCNA rescued")), 
                                                                 "homogeneous", 
                                                                 "heterogeneous"))))

    p1 <- ggplot2::ggplot(snpTable.plot) + 
            ggplot2::geom_point(ggplot2::aes(adjPos, BAF, color = class)) + 
            ggplot2::geom_vline(xintercept = cumsum.chr.size[1:22], linetype = 2) +
            ggplot2::geom_hline(yintercept = 0.5, linetype = 2, color = "gray") +
            ggplot2::geom_rect(data = snpTable.plot %>% dplyr::select(region, segment, startpos, endpos, type) %>% unique(), ggplot2::aes(xmin = startpos, xmax = endpos, ymin = -0.125, ymax = -0.025, fill = type)) +
            ggplot2::xlab("") + 
            ggplot2::scale_x_continuous(breaks = chr_tick_marks[1:22], labels = paste0('chr', c(1:22)), expand = c(0, 0)) +
            ggplot2::scale_color_manual(values = col.palette, guide = "none") +
            ggplot2::scale_fill_manual(name = "Sample level classification",
                              values = rev(fill.palette), 
                              labels = rev(c("Mirrored subclonal allelic imbalance", "Allelic imbalance rescued", "Allelic imbalance in single region", "Allelic imbalance not detected", "Haplotype information unavailable"))) +
            ggplot2::facet_wrap(~ region, ncol = 1) +
            cowplot::theme_cowplot(font_size = 12) +
            ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    segSummary.plot <- snpTable.plot %>% 
        dplyr::ungroup() %>% 
        dplyr::select(segment, startpos, endpos) %>% 
        unique() %>% 
        dplyr::left_join(segSummary, by = "segment") %>%
        dplyr::mutate(characterization = factor(characterization, levels = rev(c("MSAI", "heterogeneous", "homogeneous", "no SCNA"))))

    p2 <- ggplot2::ggplot(segSummary.plot) +
            ggplot2::geom_rect(ggplot2::aes(xmin = startpos, xmax = endpos, ymin = 0.5, ymax = 1.5, fill = characterization)) +
            ggplot2::geom_vline(xintercept = cumsum.chr.size[1:22], linetype = 2) +
            ggplot2::xlab("") + 
            ggplot2::scale_x_continuous(breaks = chr_tick_marks[1:22], labels = paste0('chr', c(1:22)), expand = c(0, 0)) +
            ggplot2::scale_y_continuous(breaks = 1, labels = "") +
            ggplot2::scale_fill_manual(name = "Overall classification",
                              values = rev(c("MSAI" = "#D95980", "heterogeneous" = "#92271B", "homogeneous" = "#4C8EA5", "no SCNA" = "#bdbdbd")),
                              labels = rev(c("MSAI" = "MSAI", "heterogeneous" = "Heterogeneous SCNA", "homogeneous" = "Homogeneous SCNA", "no SCNA" = "No SCNA detected"))) +
            cowplot::theme_cowplot(font_size = 12) +
            ggplot2::theme(legend.position = "bottom", legend.justification = "center")
    
    ### set margins
    g.p1 <- ggplot2::ggplotGrob(p1 + ggplot2::theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm")) + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()))
    g.p2 <- ggplot2::ggplotGrob(p2 + ggplot2::theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))

    ### align widths
    plot_list  <- list(g.p1, 
                       g.p2)
    all_widths <- lapply(plot_list, function(x) x$widths)
    plot_list_alignedWidths <- lapply(plot_list, function(x) {
        x$widths <- do.call(grid::unit.pmax, all_widths)
        return(x)
    })

    aa <- length(samples) * 2 + 2
    full_plot <- cowplot::ggdraw() +
        cowplot::draw_plot(plot_list_alignedWidths[[2]], x = 0, y = 0, width = 1, height = 1/aa) +
        cowplot::draw_plot(plot_list_alignedWidths[[1]], x = 0, y = 1.25/aa, width = 1, height = 1-1.25/aa)

    return(full_plot)
}