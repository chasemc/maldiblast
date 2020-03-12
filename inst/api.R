library(plumber)

#* @apiTitle Simple API
#* Echo provided text
#* @param query The text to be echoed in the response
#* @get /ech3o
function(query = "") {

  if(query == ""){
    stop("No input provided")
  }
  if (isFALSE(jsonlite::validate(query)[[1]])){
    stop("Expected json as input")
  }

  query <- jsonlite::fromJSON(query)

  if(!is.list(query)){
    stop("Input should be list of masses")
  }
  if(length(query) > 10) {
    stop("Number of queries must be less than 10")
  }
  if(any(lengths(query) > 300)) {
    stop("Each query is limited to a maximum of 300 m/z values")
  }

  if(any(sapply(unlist(query, recursive = FALSE), is.list))){
    stop("Input should be list of masses, with no nested lists")
  }
  if(!all(sapply(query, is.numeric))){
    temp <- sapply(query, class)
    stop("Input mass values should be numeric....\n ", paste0(names(temp), ": ",  temp, collapse = "\n "))
  }


  peak_blaster(query_peak_list = query,
               subject_peak_list = maldiblast:::subject_masses,
               lowerMassCutoff = 3000,
               upperMassCutoff = 15000,
               chunksize = 10,
               similarity_cutoff = NA)

}

