blaster <- function(queryPool,
                    subjectPool,
                    peakPercentPresence = 0.7,
                    lowerMassCutoff = 3000,
                    upperMassCutoff = 15000,
                    minSNR = 4,
                    tolerance = .002,
                    protein = TRUE,
                    chunksize = 10,
                    similarity_cutoff = NA){

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

  query_vectors <- IDBacApp:::createFuzzyVectorUnit(massStart = lowerMassCutoff,
                                                    massEnd = upperMassCutoff,
                                                    chunksize = chunksize,
                                                    massList = lapply(query_vectors, function(x) x@mass))

  subject_vectors <- IDBacApp:::createFuzzyVectorUnit(massStart = lowerMassCutoff,
                                                      massEnd = upperMassCutoff,
                                                      chunksize = chunksize,
                                                      massList = lapply(subject_vectors, function(x) x@mass))

  query_length <- attributes(query_vectors)$Dim[[2]] # number of cols in sparse matrix
  subject_length <- attributes(subject_vectors)$Dim[[2]] # number of cols in sparse matrix

  big_freakin_dist <- cbind(query_vectors, subject_vectors)

  big_freakin_dist <- 1 - coop::cosine(big_freakin_dist)

  big_freakin_dist <- big_freakin_dist[1:query_length, seq(query_length + 1, query_length + subject_length)]


  rows_big_freakin_dist <- attributes(big_freakin_dist)$dim[[1]] #faster nrow

  melted <- data.table::data.table(query = rep(attributes(big_freakin_dist)$dimnames[[1]], # faster row names
                                               attributes(big_freakin_dist)$dim[[2]]), # faster ncol
                                   subject = unlist(lapply(attributes(big_freakin_dist)$dimnames[[2]],# faster col names
                                                           function(x){
                                                             rep(x,
                                                                 rows_big_freakin_dist)
                                                           })),
                                   value = as.vector(big_freakin_dist),
                                   key = 'query')

  if (is.numeric(similarity_cutoff)) {
    melted <- melted[value < similarity_cutoff]
  }

  return(melted)

}
