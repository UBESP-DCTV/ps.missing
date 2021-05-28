
# Propensity Score and Missing data

This repository contains the R scripts to reproduce the analysis
workflow of the study “**Propensity Score Analysis with partially
observed baseline covariates: a practical comparison of methods for
handling missing data**” on a simulated dataset.

Before running the scripts, please install the following R packages:

``` r
install.packages(
  c("tidyverse", "BaBooN", "mice", "miceadds", "optmatch",
    "MatchIt", "WeightIt", "rms", "misaem", "twang", "CBPS", "cobalt",
    "assertive", "usethis", "here", "naniar", "tictoc", "MASS",
    "rpart", "rpart.plot", "knitr", "rmarkdown", "kableExtra",
    "gtsummary", "corrplot", "readr"), 
    dependencies = TRUE
)
```

The project contains the following folders:

-   *data*, containing the R script
    \_simultate\_synthetic\_dataset.R\_\_ that can be used to simulate
    the synthetic dataset for the analysis.

-   *analysis.R*, containing the R scripts *misc\_functions.R* and
    *main\_analysis*. The former contains the functions that are used to
    run the analysis, the latter can be used to perform the main
    analysis

-   *results.R*, containing the R markdown file *main\_results.Rmd*.
    that can be used to visualize the results of the main analysis.

The analysis workflow can be reproduced as follows:

1.  Run the script *data/simultate\_synthetic\_dataset.R*, which
    simulates the dataset and saves it in **data/synthetic\_data.rds**.

2.  Run the script *analysis/main\_analysis.R*, which runs the main
    analysis and save the results.

3.  Knitr the *results/main\_results.Rmd* file to produce a report of
    the results with tables and plots.

The workflow of the analysis is guaranteed to be reproducible with R
version *4.1.0* and the following packages versions:

``` r
library(tidyverse)
library(BaBooN)
library(mice)
library(miceadds)
library(optmatch)
library(MatchIt)
library(WeightIt)
library(rms)
library(misaem)
library(twang)
library(CBPS)
library(rmarkdown)
library(knitr)
library(kableExtra)
library(gtsummary)
library(cobalt)
library(assertive)
library(usethis)
library(here)
library(naniar)
library(tictoc)
library(MASS)
library(rpart)
library(rpart.plot)
library(corrplot)
library(readr)

si <- sessionInfo()
print(si, locale = FALSE)
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19043)
## 
## Matrix products: default
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] corrplot_0.88       rpart.plot_3.0.9    rpart_4.1-15       
##  [4] tictoc_1.0.1        naniar_0.6.1        here_1.0.1         
##  [7] usethis_2.0.1       assertive_0.3-6     cobalt_4.3.1       
## [10] gtsummary_1.4.1     kableExtra_1.3.4    knitr_1.33         
## [13] rmarkdown_2.8       CBPS_0.22           glmnet_4.1-1       
## [16] Matrix_1.3-3        numDeriv_2016.8-1.1 nnet_7.3-16        
## [19] MASS_7.3-54         twang_2.0           misaem_1.0.1       
## [22] rms_6.2-0           SparseM_1.81        Hmisc_4.5-0        
## [25] Formula_1.2-4       lattice_0.20-44     WeightIt_0.12.0    
## [28] MatchIt_4.1.0       optmatch_0.9-13     survival_3.2-11    
## [31] miceadds_3.11-6     mice_3.13.0         BaBooN_0.2-0       
## [34] Rcpp_1.0.6          forcats_0.5.1       stringr_1.4.0      
## [37] dplyr_1.0.6         purrr_0.3.4         readr_1.4.0        
## [40] tidyr_1.1.3         tibble_3.1.2        ggplot2_3.3.3      
## [43] tidyverse_1.3.1    
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1               backports_1.2.1           
##   [3] systemfonts_1.0.2          assertive.files_0.0-2     
##   [5] splines_4.1.0              TH.data_1.0-10            
##   [7] digest_0.6.27              foreach_1.5.1             
##   [9] htmltools_0.5.1.1          svd_0.5                   
##  [11] fansi_0.4.2                magrittr_2.0.1            
##  [13] checkmate_2.0.0            assertive.datetimes_0.0-3 
##  [15] assertive.numbers_0.0-2    cluster_2.1.2             
##  [17] modelr_0.1.8               matrixStats_0.58.0        
##  [19] sandwich_3.0-1             svglite_2.0.0             
##  [21] jpeg_0.1-8.1               colorspace_2.0-1          
##  [23] rvest_1.0.0                mitools_2.4               
##  [25] assertive.strings_0.0-3    haven_2.4.1               
##  [27] xfun_0.23                  crayon_1.4.1              
##  [29] jsonlite_1.7.2             zoo_1.8-9                 
##  [31] iterators_1.0.13           glue_1.4.2                
##  [33] gtable_0.3.0               webshot_0.5.2             
##  [35] MatrixModels_0.5-0         shape_1.4.6               
##  [37] abind_1.4-5                scales_1.1.1              
##  [39] mvtnorm_1.1-1              DBI_1.1.1                 
##  [41] assertive.data.uk_0.0-2    assertive.models_0.0-2    
##  [43] assertive.code_0.0-3       viridisLite_0.4.0         
##  [45] xtable_1.8-4               htmlTable_2.2.1           
##  [47] foreign_0.8-81             assertive.data.us_0.0-2   
##  [49] survey_4.0                 htmlwidgets_1.5.3         
##  [51] httr_1.4.2                 RColorBrewer_1.1-2        
##  [53] ellipsis_0.3.2             pkgconfig_2.0.3           
##  [55] dbplyr_2.1.1               utf8_1.2.1                
##  [57] tidyselect_1.1.1           rlang_0.4.11              
##  [59] munsell_0.5.0              cellranger_1.1.0          
##  [61] tools_4.1.0                xgboost_1.4.1.1           
##  [63] cli_2.5.0                  generics_0.1.0            
##  [65] assertive.reflection_0.0-5 broom_0.7.6               
##  [67] evaluate_0.14              yaml_2.2.1                
##  [69] fs_1.5.0                   assertive.matrices_0.0-2  
##  [71] visdat_0.5.3               assertive.sets_0.0-3      
##  [73] nlme_3.1-152               quantreg_5.85             
##  [75] RItools_0.1-17             xml2_1.3.2                
##  [77] compiler_4.1.0             rstudioapi_0.13           
##  [79] png_0.1-7                  gt_0.3.0                  
##  [81] reprex_2.0.0               broom.helpers_1.3.0       
##  [83] stringi_1.6.2              assertive.base_0.0-9      
##  [85] gbm_2.1.8                  assertive.data_0.0-3      
##  [87] vctrs_0.3.8                pillar_1.6.1              
##  [89] norm_1.0-9.5               lifecycle_1.0.0           
##  [91] data.table_1.14.0          conquer_1.0.2             
##  [93] assertive.types_0.0-3      R6_2.5.0                  
##  [95] latticeExtra_0.6-29        assertive.properties_0.0-4
##  [97] gridExtra_2.3              codetools_0.2-18          
##  [99] polspline_1.1.19           assertthat_0.2.1          
## [101] rprojroot_2.0.2            withr_2.4.2               
## [103] multcomp_1.4-17            hms_1.1.0                 
## [105] grid_4.1.0                 coda_0.19-4               
## [107] lubridate_1.7.10           base64enc_0.1-3
```
