
source("figure1.R") ## Figure 1

library(splines)

source("presence_absence_only/prevalenceExampleDataSetup.R")

source("presence_absence_only/prevalence_smallNoSpline.R") ## Figure 2a, 2e
source("presence_absence_only/prevalence_small.R") ## Figure 2b, 2f
source("presence_absence_only/prevalence_medium.R") ## Figure 2c, 2g, 3a, 3b
source("presence_absence_only/prevalence_large.R") ## Figure 2d, 2h
## note starts to get slow here ^

library(unmarked) ## double visit
library(detect) ## single visit

source("single_double_visit/occurDetectDataSetup.R")

source("single_double_visit/occurDetect_small_noSpline.R") ## Figure 4a, 4e, 5i, 5m
source("single_double_visit/occurDetect_small.R") ## Figure 4b, 4f, 4j, 4n
source("single_double_visit/occurDetect_medium.R") ## Figure 4c, 4g, 4k, 4o
source("single_double_visit/occurDetect_large.R") ## Figure 4d, 4h, 4l, 4p
## this will be slow ^

## Supporting Materials

source("single_double_visit_abundance/abunDetectDataSetup.R")
source("single_double_visit_abundance/abunDetect_small_noSpline.R") ## Figure S1e, S2e, S1a, S2a
source("single_double_visit_abundance/abunDetect_small.R") ## Figure S1f, S2f, S1b, S2b
source("single_double_visit_abundance/abunDetect_medium.R") ## Figure S1g, S2g, S1c, S2c
## starts to slow down ^
source("single_double_visit_abundance/abunDetect_large.R") ## Figure S1h, S2h, S1d, S2d
## this is REALLY slow
