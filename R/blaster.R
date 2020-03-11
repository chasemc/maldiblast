#' "BLAST" MALDI-TOF MS peak lists from IDBac databases
#'
#' @param queryPool IDBac connecttion to query database (peaks to search)
#' @param subjectPool IDBac connecttion to subject database (peaks to search against)
#' @param peakPercentPresence percent of replicates a peak must be present in (on scale of 0 to 1)
#' @param lowerMassCutoff masses below this number will not be considered
#' @param upperMassCutoff masses above this number will not be considered
#' @param minSNR masses with SNR below this number will not be considered
#' @param protein  TRUE/FALSE protein or small molecule peaks are being compared
#' @param chunksize numeric, bin spectrum into X-sized bins
#' @param similarity_cutoff numeric to remove similarities above threshold, NA to return everythingg (default)
#'
#' @inheritParams MALDIquant::binPeaks
#'
#' @return data.table of pairwise similarities between query and subject peak lists
#' @export
#'

idbac_blaster <- function(queryPool,
                          subjectPool,
                          peakPercentPresence = 0.7,
                          lowerMassCutoff = 3000,
                          upperMassCutoff = 15000,
                          minSNR = 4,
                          tolerance = .002,
                          protein = TRUE,
                          chunksize = 10,
                          similarity_cutoff = NA){


  result <- idbac_get_peaks(queryPool = queryPool,
                            subjectPool = subjectPool,
                            peakPercentPresence = peakPercentPresence,
                            lowerMassCutoff = lowerMassCutoff,
                            upperMassCutoff = upperMassCutoff,
                            minSNR = minSNR,
                            tolerance = tolerance,
                            protein = protein)

  peak_blaster(query_peak_list = lapply(result$query_vectors, function(x) x@mass),
               subject_peak_list = lapply(result$subject_vectors, function(x) x@mass),
               lowerMassCutoff = lowerMassCutoff,
               upperMassCutoff = upperMassCutoff,
               chunksize = chunksize)


}


#' Get peaks from IDBac query and subject databases
#'
#' @inheritParams .create_large_sparse
#' @inheritParams idbac_blaster
#' @return query and subject peak lists
#' @export
#'

idbac_get_peaks <- function(queryPool,
                            subjectPool,
                            peakPercentPresence = 0.7,
                            lowerMassCutoff = 3000,
                            upperMassCutoff = 15000,
                            minSNR = 4,
                            tolerance = .002,
                            protein = TRUE){

  query_vectors <- IDBacApp::idbac_get_peaks(pool = queryPool,
                                             sampleIDs = IDBacApp::idbac_available_samples(queryPool, allSamples = T),
                                             peakPercentPresence = peakPercentPresence,
                                             lowerMassCutoff = lowerMassCutoff,
                                             upperMassCutoff = upperMassCutoff,
                                             minSNR = minSNR,
                                             tolerance = tolerance,
                                             protein = protein)

  subject_vectors <- IDBacApp::idbac_get_peaks(pool = subjectPool,
                                               sampleIDs = IDBacApp::idbac_available_samples(subjectPool, allSamples = T),
                                               peakPercentPresence = peakPercentPresence,
                                               lowerMassCutoff = lowerMassCutoff,
                                               upperMassCutoff = upperMassCutoff,
                                               minSNR = minSNR,
                                               tolerance = tolerance,
                                               protein = protein)
  return(list(query_vectors = query_vectors,
              subject_vectors = subject_vectors))

}





#' MALDI-TOF MS peak "BLAST"er
#'
#' @inheritParams .create_large_sparse
#' @inheritParams idbac_blaster
#'
#' @return data.table of pairwise similarities between query and subject peak lists
#' @export

peak_blaster <- function(query_peak_list,
                         subject_peak_list,
                         lowerMassCutoff = 3000,
                         upperMassCutoff = 15000,
                         chunksize = 10,
                         similarity_cutoff = NA){



  create_large_sparse_output <- .create_large_sparse(query_peak_list = query_peak_list,
                                                     subject_peak_list = subject_peak_list,
                                                     lowerMassCutoff = lowerMassCutoff,
                                                     upperMassCutoff = upperMassCutoff,
                                                     chunksize = chunksize)

  big_dist <- .create_dist_from_sparse(create_large_sparse_output = create_large_sparse_output)

  .trim_and_melt_dist(big_dist = big_dist,
                      create_large_sparse_output = create_large_sparse_output,
                      similarity_cutoff = similarity_cutoff)


}




#' Create sparse matrix from two sets of peak lists
#'
#' @param query_peak_list list of numeric vectors representing query peaks
#' @param subject_peak_list list of numeric vectors representing subject peaks
#'
#' @inheritParams idbac_blaster
#'
#' @return list of length 3, containing length of queries, length of subjects, and the sparse matrix
#'
.create_large_sparse <- function(query_peak_list,
                                 subject_peak_list,
                                 lowerMassCutoff,
                                 upperMassCutoff,
                                 chunksize){
  query_vectors <- IDBacApp:::createFuzzyVectorUnit(massStart = lowerMassCutoff,
                                                    massEnd = upperMassCutoff,
                                                    chunksize = chunksize,
                                                    massList = query_peak_list)

  subject_vectors <- IDBacApp:::createFuzzyVectorUnit(massStart = lowerMassCutoff,
                                                      massEnd = upperMassCutoff,
                                                      chunksize = chunksize,
                                                      massList = subject_peak_list)

  return(list(query_length = attributes(query_vectors)$Dim[[2]], # number of cols in sparse matrix
              subject_length = attributes(subject_vectors)$Dim[[2]], # number of cols in sparse matrix
              sparse_matrix = cbind(query_vectors, subject_vectors)))

}


#' Calculate cosine distance from the sparse matrix
#'
#' @param create_large_sparse_output output from create_large_sparse()
#' @importFrom coop cosine
#' @return square matrix
.create_dist_from_sparse <- function(create_large_sparse_output){

  1 - cosine(create_large_sparse_output$sparse_matrix)

}


#' Melt distance matrix into pairwise data.table
#'
#' @param big_dist output from create_dist_from_sparse()
#' @inheritParams .create_dist_from_sparse
#' @inheritParams idbac_blaster
#' @return pairwise distances as data.table
.trim_and_melt_dist <- function(big_dist,
                                create_large_sparse_output,
                                similarity_cutoff){

  big_dist <- big_dist[1:create_large_sparse_output$query_length, seq(create_large_sparse_output$query_length + 1,
                                                                      create_large_sparse_output$query_length + create_large_sparse_output$subject_length)]

  #big_dist <-as.matrix(big_dist)
  rows_big_dist <- attributes(big_dist)$dim[[1]] #faster nrow

  melted <- data.table::data.table(query = rep(attributes(big_dist)$dimnames[[1]], # faster row names
                                               attributes(big_dist)$dim[[2]]), # faster ncol
                                   subject = unlist(lapply(attributes(big_dist)$dimnames[[2]],# faster col names
                                                           function(x){
                                                             rep(x,
                                                                 rows_big_dist)
                                                           })),
                                   value = as.vector(big_dist))

  if (is.numeric(similarity_cutoff)) {
    melted <- melted[melted$value < similarity_cutoff]
  }

  return(melted)
}
