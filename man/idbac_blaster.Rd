% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blaster.R
\name{idbac_blaster}
\alias{idbac_blaster}
\title{"BLAST" MALDI-TOF MS peak lists from IDBac databases}
\usage{
idbac_blaster(
  queryPool,
  subjectPool,
  peakPercentPresence = 0.7,
  lowerMassCutoff = 3000,
  upperMassCutoff = 15000,
  minSNR = 4,
  tolerance = 0.002,
  protein = TRUE,
  chunksize = 10,
  similarity_cutoff = NA
)
}
\arguments{
\item{queryPool}{IDBac connecttion to query database (peaks to search)}

\item{subjectPool}{IDBac connecttion to subject database (peaks to search against)}

\item{peakPercentPresence}{percent of replicates a peak must be present in (on scale of 0 to 1)}

\item{lowerMassCutoff}{masses below this number will not be considered}

\item{upperMassCutoff}{masses above this number will not be considered}

\item{minSNR}{masses with SNR below this number will not be considered}

\item{tolerance}{
    \code{double}, maximal relative deviation of a peak position (mass) to be
    considered as identical. Must be multiplied by \code{10^-6} for ppm,
    e.g. use \code{tolerance=5e-6} for 5 ppm.
  }

\item{protein}{TRUE/FALSE protein or small molecule peaks are being compared}

\item{chunksize}{numeric, bin spectrum into X-sized bins}

\item{similarity_cutoff}{numeric to remove similarities above threshold, NA to return everythingg (default)}
}
\value{
data.table of pairwise similarities between query and subject peak lists
}
\description{
"BLAST" MALDI-TOF MS peak lists from IDBac databases
}
