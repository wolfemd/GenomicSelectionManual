# Preliminary field trial analysis



-   **Context and Purpose:**

-   **Upstream:** Section \@ref(prepare-genotypic-data) - quality control steps for genotypic data

-   **Downstream:** Genomic prediction and related analyses.

-   **Inputs:**

-   **Expected outputs:**

## One-stage or multi-stage?

We often have very large, multi-year, multi-location, multi-trial-type (**MET**) datasets we use to train genomic prediction models. The number of plots can be in the range of 10- or even 100,000 plots with many thousands of unique genotypes observed in unbalanced fashion across heterogenous experimental designs. All that to say, the computational burden and level of expertise required to execute genomic prediction analyses directly on these data is very great.

For the sake of semi-automation and computational efficiency, the standard genomic prediction pipeline I have implemented for NextGen Cassava, which I demonstrate below, has two stages:

**Stage 1. Get BLUPs**

-   Conduct preliminary analysis of the trial data *without* genomic relatedeness / marker-information.

    Identify the best-fitting model for the data and potentially curate the raw data, removing outliers.

-   Fit mixed-model to **MET** data. Extract **`BLUPs`**, **`PEVs`** and variance components (**`VarComps`**). Compute de-regressed BLUPs (**`drgBLUPs`**) and weights (**`WTS`**) for "Stage 2".

**Stage 2. Cross-validation and genomic prediction**

Conduct weighted cross-validation and genomic prediction analyses using **`de-regressed BLUPs`** (alternative **`BLUEs`**) as the response variable and the weights from Stage 1. This effectively reduces the number of "observations" in the analysis to the number of unique genotypes (i.e. clones, inbred lines, etc.), leading to lower computational burden.

**Note of advice for single-stage analyses:** I strongly recommend conducting preliminary analysis of the trial data *without* genomic relatedeness / marker-information. Ensure the model will converge, residuals look acceptable and otherwise assess the best fitting model *before* commiting to lengthy computationally intensive analyses. :)

## Set-up training datasets


```r
phenos<-readRDS(here::here("output","phenotypes_cleaned.rds"))
```

The pipeline version of this analysis will use the `TRUE/FALSE` values of **`CompleteBlocks`** and **`IncompleteBlocks`** ([Preliminary analysis of trial data][#detect_designs]).


```r
phenos %>% 
     count(CompleteBlocks,IncompleteBlocks,locationName) %>% 
     spread(locationName,n)
#> # A tibble: 3 × 4
#>   CompleteBlocks IncompleteBlocks Ibadan Ubiaja
#>   <lgl>          <lgl>             <int>  <int>
#> 1 FALSE          TRUE                777     NA
#> 2 TRUE           FALSE                NA    125
#> 3 TRUE           TRUE                407    414
```

Convert the data to "long format" . Remove missing values. "Nest" the data by Trait. 


```r
traits<-c("DM","MCMDS","logFYLD","logDYLD")
phenos<-phenos %>% 
     # Convert the data to "long format" 
     pivot_longer(cols = all_of(traits), 
                  names_to = "Trait", values_to = "Value") %>%
     # Remove missing values
     filter(!is.na(Value)) %>%
     # Nest the MultiEnvironmentTrial data by trait
     nest(METdata=c(-Trait))
phenos %>% 
     mutate(N_plots=map_dbl(METdata,nrow))
#> # A tibble: 4 × 3
#>   Trait   METdata               N_plots
#>   <chr>   <list>                  <dbl>
#> 1 DM      <tibble [1,448 × 35]>    1448
#> 2 MCMDS   <tibble [1,436 × 35]>    1436
#> 3 logFYLD <tibble [1,590 × 35]>    1590
#> 4 logDYLD <tibble [1,443 × 35]>    1443
```

Where previously there was one row per plot and a large number of columns, now things are simple and tidy. 

One row per trait. The actual plot-basis data are now contained within **METdata**, a list-type column, each element containing a `tibble`.

To demonstrate, check the contents of row 1 of the **METdata** column:


```r
phenos$METdata[[1]] %>% head
#> # A tibble: 6 × 35
#>   studyYear programName locationName studyName   studyDesign
#>       <int> <chr>       <chr>        <chr>       <chr>      
#> 1      2019 IITA        Ubiaja       19.GS.C1.C… Alpha      
#> 2      2019 IITA        Ubiaja       19.GS.C1.C… Alpha      
#> 3      2019 IITA        Ubiaja       19.GS.C1.C… Alpha      
#> 4      2019 IITA        Ubiaja       19.GS.C1.C… Alpha      
#> 5      2019 IITA        Ubiaja       19.GS.C1.C… Alpha      
#> 6      2019 IITA        Ubiaja       19.GS.C1.C… Alpha      
#> # … with 30 more variables: plotWidth <int>,
#> #   plotLength <dbl>, fieldSize <dbl>, plantingDate <chr>,
#> #   harvestDate <chr>, germplasmName <chr>,
#> #   observationUnitDbId <int>, replicate <int>,
#> #   blockNumber <int>, plotNumber <int>, rowNumber <int>,
#> #   colNumber <int>, entryType <chr>, trialType <chr>,
#> #   numberBlocks <int>, numberReps <int>, …
```

## Compare models

There is a "standard" model I have deployed reliably on Cassava MET data for the GS pipeline. 

For the sake of demonstration, below I pick out one trait, and fit 3 models to it.

In the past, I have used both `lmer()` (`lme4` package) and `asreml()` (asreml R package) for this sort of analysis. `lmer()` is limited to homogenous error variance structures and only certain design models / variance structures. `asreml()` 


```r
library(sommer)
#> Loading required package: Matrix
#> 
#> Attaching package: 'Matrix'
#> The following objects are masked from 'package:tidyr':
#> 
#>     expand, pack, unpack
#> Loading required package: MASS
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
#> Loading required package: lattice
#> Loading required package: crayon
#> 
#> Attaching package: 'crayon'
#> The following object is masked from 'package:ggplot2':
#> 
#>     %+%
METdata<-phenos$METdata[[1]] %>% 
     mutate(CompleteBlocks=as.character(CompleteBlocks),
            IncompleteBlocks=as.character(IncompleteBlocks))
m1<-mmer(Value~yearInLoc,
         random=~germplasmName + 
              vs(at(CompleteBlocks,"TRUE"),repInTrial) + 
              vs(at(IncompleteBlocks,"TRUE"),blockInRep),
         rcov=~vs(ds(studyName),units),
         data=METdata, 
         verbose = FALSE)
```


```r
summary(m1)
#> ============================================================================
#>                  Multivariate Linear Mixed Model fit by REML                 
#> ******************************  sommer 4.1  ****************************** 
#> ============================================================================
#>         logLik     AIC      BIC Method Converge
#> Value -324.346 652.692 663.2478     NR     TRUE
#> ============================================================================
#> Variance-Covariance components:
#>                                            VarComp
#> germplasmName.Value-Value                    5.220
#> TRUE:repInTrial.Value-Value                 11.186
#> TRUE:blockInRep.Value-Value                  1.330
#> 19.GS.C1.C2.C3.AYT.42.UB:units.Value-Value   8.353
#> 19.GS.C2.UYT.36.setA.IB:units.Value-Value    8.322
#> 19.GS.C2.UYT.36.setA.UB:units.Value-Value    4.133
#> 19.GS.C2.UYT.36.setB.IB:units.Value-Value    4.376
#> 19.GS.C2.UYT.36.setB.UB:units.Value-Value    5.441
#> 19.GS.C4B.PYT.135.UB:units.Value-Value      11.149
#> 19.GS.C4B.PYT.140.IB:units.Value-Value       3.918
#> 19.GS.C4C.CE.822.IB:units.Value-Value        6.478
#>                                            VarCompSE Zratio
#> germplasmName.Value-Value                     0.5807  8.989
#> TRUE:repInTrial.Value-Value                   4.5829  2.441
#> TRUE:blockInRep.Value-Value                   0.2954  4.503
#> 19.GS.C1.C2.C3.AYT.42.UB:units.Value-Value    1.6652  5.016
#> 19.GS.C2.UYT.36.setA.IB:units.Value-Value     1.6767  4.963
#> 19.GS.C2.UYT.36.setA.UB:units.Value-Value     0.9170  4.507
#> 19.GS.C2.UYT.36.setB.IB:units.Value-Value     0.9845  4.445
#> 19.GS.C2.UYT.36.setB.UB:units.Value-Value     1.1190  4.862
#> 19.GS.C4B.PYT.135.UB:units.Value-Value        1.3080  8.524
#> 19.GS.C4B.PYT.140.IB:units.Value-Value        0.4906  7.987
#> 19.GS.C4C.CE.822.IB:units.Value-Value         0.7080  9.149
#>                                            Constraint
#> germplasmName.Value-Value                    Positive
#> TRUE:repInTrial.Value-Value                  Positive
#> TRUE:blockInRep.Value-Value                  Positive
#> 19.GS.C1.C2.C3.AYT.42.UB:units.Value-Value   Positive
#> 19.GS.C2.UYT.36.setA.IB:units.Value-Value    Positive
#> 19.GS.C2.UYT.36.setA.UB:units.Value-Value    Positive
#> 19.GS.C2.UYT.36.setB.IB:units.Value-Value    Positive
#> 19.GS.C2.UYT.36.setB.UB:units.Value-Value    Positive
#> 19.GS.C4B.PYT.135.UB:units.Value-Value       Positive
#> 19.GS.C4B.PYT.140.IB:units.Value-Value       Positive
#> 19.GS.C4C.CE.822.IB:units.Value-Value        Positive
#> ============================================================================
#> Fixed effects:
#>   Trait                    Effect Estimate Std.Error
#> 1 Value               (Intercept)   38.147    0.2314
#> 2 Value yearInLocIITA_Ubiaja_2019   -5.293    1.2302
#>   t.value
#> 1 164.880
#> 2  -4.302
#> ============================================================================
#> Groups and observations:
#>                 Value
#> germplasmName     811
#> TRUE:repInTrial    15
#> TRUE:blockInRep   126
#> ============================================================================
#> Use the '$' sign to access results and parameters
```


```r
m2<-mmer(Value~yearInLoc,
         random=~germplasmName + 
              vs(at(CompleteBlocks,"TRUE"),repInTrial) + 
              vs(at(IncompleteBlocks,"TRUE"),blockInRep),
         data=METdata, 
         verbose = FALSE)

m3<-mmer(Value~yearInLoc,
         random=~germplasmName + repInTrial + blockInRep,
         data=METdata, 
         verbose = FALSE)
```


```r
anova(m1,m2)
#> Likelihood ratio test for mixed models
#> ==============================================================
#>      Df      AIC      BIC     loLik    Chisq ChiDf
#> mod1 17 652.6920 663.2478 -324.3460               
#> mod2 10 701.9414 712.4973 -348.9707 49.24947     7
#>                      PrChisq
#> mod1                        
#> mod2 2.0273987089018e-08 ***
#> ==============================================================
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>      Df      AIC      BIC     loLik    Chisq ChiDf
#> mod1 17 652.6920 663.2478 -324.3460               
#> mod2 10 701.9414 712.4973 -348.9707 49.24947     7
#>                      PrChisq
#> mod1                        
#> mod2 2.0273987089018e-08 ***
anova(m1,m3)
#> Likelihood ratio test for mixed models
#> ==============================================================
#>      Df      AIC      BIC     loLik    Chisq ChiDf
#> mod1 17 652.6920 663.2478 -324.3460               
#> mod2 10 695.2303 705.7862 -345.6152 42.53834     7
#>                       PrChisq
#> mod1                         
#> mod2 4.09508831534527e-07 ***
#> ==============================================================
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>      Df      AIC      BIC     loLik    Chisq ChiDf
#> mod1 17 652.6920 663.2478 -324.3460               
#> mod2 10 695.2303 705.7862 -345.6152 42.53834     7
#>                       PrChisq
#> mod1                         
#> mod2 4.09508831534527e-07 ***
anova(m2,m3)
#> Likelihood ratio test for mixed models
#> ==============================================================
#>      Df      AIC      BIC     loLik   Chisq ChiDf PrChisq
#> mod1 10 701.9414 712.4973 -348.9707                      
#> mod2 10 695.2303 705.7862 -345.6152 6.71113     0   0 ***
#> ==============================================================
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>      Df      AIC      BIC     loLik   Chisq ChiDf PrChisq
#> mod1 10 701.9414 712.4973 -348.9707                      
#> mod2 10 695.2303 705.7862 -345.6152 6.71113     0   0 ***
```

