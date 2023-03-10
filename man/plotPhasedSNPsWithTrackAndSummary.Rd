% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting-functions.R
\name{plotPhasedSNPsWithTrackAndSummary}
\alias{plotPhasedSNPsWithTrackAndSummary}
\title{Plotting phasing of SNPs across the genome for multiple samples with track per sample and summary track}
\usage{
plotPhasedSNPsWithTrackAndSummary(
  newPhasingTable,
  genome,
  samples = NULL,
  col.palette = setNames(c("#f1a340", "#998ec3", "#bdbdbd"), c("top", "bottom", "mixed")),
  fill.palette = setNames(c("#D95980", "#335120", "#9DCC5F", "#63AAC0", "#bdbdbd"),
    c("MSAI", "SCNA rescued", "called", "not detected", "not called"))
)
}
\arguments{
\item{newPhasingTable}{A table of SNPs with their positions along with their characterization.}

\item{genome}{A genome version to use, for example BSgenome.Hsapiens.UCSC.hg19.}

\item{samples}{A vector of sample names to include for plotting in the order they appear in.
Defaults to NULL in which case all samples are considered.}

\item{col.palette}{A named vector of colors for plotting the SNPs.
Defaults to orange and purple for allele A and B, respectively.}

\item{fill.palette}{A named vector of colors for plotting the characterization of each segment.
Defaults to orange and purple for allele A and B, respectively.}
}
\value{
A ggplot object that can be plotted directly or saved as a pdf.
}
\description{
This function plots the phased SNPs for all samples of interest.
Additionally, for each sample a track is added at the bottom of the plot
highlighting the classification of the segment within that sample.
Standard classification includes not called, not detected, called, SCNA rescued, MSAI.
At the bottom of the plot a summary track is included highlighting whether each segment
has an SCNA detected in any region, whether this is homogeneous and detected or rescued in all regions,
whether it is heterogeneous and only detected in a subset of regions, or whether there is MSAI.
}
