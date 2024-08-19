# Data Availability
The dataset generated and analyzed during the current study is not publicly available due to [reason, e.g., privacy concerns, proprietary restrictions, etc.], but it will be made available upon reasonable request. Researchers who are interested in accessing the dataset for replication or further study can contact the corresponding author at westshine28@snu.ac.kr.

# Simulation Studies
## Requirements
* simulation study requires [wdm, locfdr, VineCopula, EnvStats, BiocManager, devtools, reticulate, mixtools, truncnorm, splines] packages from cran R, [IHW, swfdr] package from Bioconductor and [RadaFDR, DESeq2] package from github.
* To download packages in cran R, use the following command : install.packages("packageName")
* To download package in Bioconductor, call Bioconductor first and then use the following commmand : BiocManager::install("packageName")
* To download package in Github, use the following command : devtools::install_github("fxia22/RadaFDR", force = TRUE)
* Installing all packages would take ~1 hour.


```R
library(wdm)
library(locfdr)
library(VineCopula)
library(EnvStats)
library(BiocManager)
library(devtools)
library(reticulate)
library(IHW)
library(swfdr)
library(RadaFDR)
library(truncnorm)
library(mixtools)
library(splines)
library(DESeq2)
```

    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
        tapply, union, unique, unsplit, which.max, which.min
    
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following object is masked from ‘package:utils’:
    
        findMatches
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following object is masked from ‘package:EnvStats’:
    
        iqr
    
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    



```R
sessionInfo()
```


    R version 4.4.0 (2024-04-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 22.04.4 LTS
    
    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    time zone: Etc/UTC
    tzcode source: system (glibc)
    
    attached base packages:
    [1] stats4    splines   stats     graphics  grDevices utils     datasets 
    [8] methods   base     
    
    other attached packages:
     [1] DESeq2_1.44.0               SummarizedExperiment_1.34.0
     [3] Biobase_2.64.0              MatrixGenerics_1.16.0      
     [5] matrixStats_1.3.0           GenomicRanges_1.56.1       
     [7] GenomeInfoDb_1.40.1         IRanges_2.38.1             
     [9] S4Vectors_0.42.1            BiocGenerics_0.50.0        
    [11] mixtools_2.0.0              truncnorm_1.0-9            
    [13] reticulate_1.38.0           RadaFDR_0.1.5              
    [15] IHW_1.32.0                  locfdr_1.1-8               
    [17] wdm_0.2.4                   devtools_2.4.5             
    [19] usethis_2.2.2               swfdr_1.30.0               
    [21] BiocManager_1.30.23         EnvStats_2.8.1             
    [23] VineCopula_2.5.0           
    
    loaded via a namespace (and not attached):
     [1] remotes_2.4.2.1         fdrtool_1.2.17          rlang_1.1.2            
     [4] magrittr_2.0.3          compiler_4.4.0          png_0.1-8              
     [7] callr_3.7.3             vctrs_0.6.5             stringr_1.5.1          
    [10] profvis_0.3.8           pkgconfig_2.0.3         crayon_1.5.2           
    [13] fastmap_1.1.1           XVector_0.44.0          ellipsis_0.3.2         
    [16] utf8_1.2.4              promises_1.2.1          sessioninfo_1.2.2      
    [19] UCSC.utils_1.0.0        ps_1.7.5                purrr_1.0.2            
    [22] zlibbioc_1.50.0         cachem_1.0.8            jsonlite_1.8.8         
    [25] later_1.3.2             DelayedArray_0.30.1     BiocParallel_1.38.0    
    [28] uuid_1.1-1              parallel_4.4.0          R6_2.5.1               
    [31] stringi_1.8.3           pkgload_1.3.3           Rcpp_1.0.12            
    [34] IRkernel_1.3.2          base64enc_0.1-3         httpuv_1.6.13          
    [37] Matrix_1.4-0            tidyselect_1.2.1        abind_1.4-5            
    [40] codetools_0.2-18        miniUI_0.1.1.1          curl_5.2.0             
    [43] processx_3.8.3          pkgbuild_1.4.3          lattice_0.20-45        
    [46] tibble_3.2.1            shiny_1.8.0             evaluate_0.23          
    [49] desc_1.4.3              survival_3.2-13         urlchecker_1.0.1       
    [52] kernlab_0.9-32          pillar_1.9.0            lpsymphony_1.32.0      
    [55] plotly_4.10.4           generics_0.1.3          IRdisplay_1.1          
    [58] ggplot2_3.5.0           munsell_0.5.0           scales_1.3.0           
    [61] xtable_1.8-4            glue_1.6.2              slam_0.1-52            
    [64] lazyeval_0.2.2          ADGofTest_0.3           tools_4.4.0            
    [67] data.table_1.15.4       locfit_1.5-9.10         pbdZMQ_0.3-10          
    [70] fs_1.6.3                mvtnorm_1.2-5           grid_4.4.0             
    [73] tidyr_1.3.1             colorspace_2.1-0        nlme_3.1-155           
    [76] GenomeInfoDbData_1.2.12 repr_1.1.6              cli_3.6.2              
    [79] fansi_1.0.6             S4Arrays_1.4.1          segmented_2.1-1        
    [82] viridisLite_0.4.2       dplyr_1.1.4             gtable_0.3.4           
    [85] digest_0.6.33           SparseArray_1.4.8       htmlwidgets_1.6.4      
    [88] memoise_2.0.1           htmltools_0.5.7         lifecycle_1.0.4        
    [91] httr_1.4.7              mime_0.12               MASS_7.3-61            

