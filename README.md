# Work in Progress
 This repo may disappear  at any time


API can be started with `plumb(file='R/api.R')$run()`, requires R {plumber}



```{r}
example_pool <- IDBacApp::idbac_connect("apostle_islands","C:/Users/chase/Documents/GitHub/sponge_meta/data/sqlite")[[1]]
example_pool2 <- IDBacApp::idbac_connect("two_sponge_project","C:/Users/chase/Documents/GitHub/sponge_meta/data/sqlite")[[1]]
```



```{r}


z <- maldiblast::blaster(queryPool = example_pool,
                         subjectPool = example_pool2,
                         peakPercentPresence = 0.6,
                         lowerMassCutoff = 3000,
                         upperMassCutoff = 15000,
                         minSNR = 4,
                         tolerance = .002,
                         protein = TRUE,
                         chunksize = 10,
                         similarity_cutoff = NA)
```


 temp <- .retrieve_peaks_from_pool(pool = example_pool,
                                    sampleIDs =  IDBacApp::idbac_available_samples(example_pool, allSamples = T)[1:5],
                                    protein = T,
                                    minSNR = 5)


```{r}
z <- z[order(z$value, decreasing = F)][, head(.SD, 3), by = "query"]

z <- split(z, z$query)
```


```{r}
my_plot <- IDBacApp:::assembleMirrorPlots(sampleID1 = "SP1155",
                                          sampleID2 = "tsp_927",
                                          peakPercentPresence = 0.6,
                                          lowerMassCutoff = 4000,
                                          upperMassCutoff = 15000,
                                          minSNR = 10, 
                                          tolerance = 0.002,
                                          pool1 = example_pool,
                                          pool2 = example_pool2,
                                          normalizeSpectra = F)
IDBacApp:::mirrorPlot(my_plot)
```

