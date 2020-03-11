# Work in Progress
 This repo may disappear  at any time

```{r}
source("R/blaster.R")
```


```{r}
#example_pool <- IDBacApp::idbac_connect("VN2018_update","C:\\Users\\chase\\Documents\\GitHub\\vietnam\\data")[[1]]



example_pool <- IDBacApp::idbac_connect("apostle_islands","C:/Users/chase/Documents/GitHub/sponge_meta/data/sqlite")[[1]]
example_pool2 <- IDBacApp::idbac_connect("two_sponge_project","C:/Users/chase/Documents/GitHub/sponge_meta/data/sqlite")[[1]]
```



```{r}


z <- blaster(queryPool = example_pool,
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

