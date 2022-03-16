# Standard K-fold Cross-validation



-   **Context and Purpose:**

-   **Upstream:** Section \@ref() -

-   **Downstream:**

-   **Inputs:**

-   **Expected outputs:**

In this section we will run K-fold cross-validation to evaluate the accuracy of predicting the performance of candidate parents (**GEBV**) who have not been phenotyped.

This is always recommended, and there are alternative kinds of predictions that could be set-up to measure this.

Important is the distinction from analyses we will do downstream to assess the accuracy of predicting the performance of *crosses* (i.e. mates).

We will use the `runCrossVal()` function.

We will demonstrate a few of the additional features that it provides in the process:

1.  Support for multiple traits
2.  Computing selection index accuracy

Finally, we'll make a simple plot of the results.

## Process Map

![](images/kfold_crossval_process_map.png){width=100%}

## Set-up for the cross-validation


```r
blups<-readRDS(here::here("output","blups.rds"))
A<-readRDS(file=here::here("output","kinship_add.rds"))
```


```r
blups %<>% 
     # need to rename the "blups" list to comply with the runCrossVal function
     rename(TrainingData=blups) %>% 
     dplyr::select(Trait,TrainingData) %>% 
     # need also to remove phenotyped-but-not-genotyped lines
     # couldn't hurt to also subset the kinship to only phenotyped lines... would save RAM
     mutate(TrainingData=map(TrainingData,
                             ~filter(.,germplasmName %in% rownames(A)) %>% 
                                  # rename the germplasmName column to GID
                                  rename(GID=germplasmName)))

blups
#> # A tibble: 4 × 2
#>   Trait   TrainingData      
#>   <chr>   <list>            
#> 1 DM      <tibble [346 × 6]>
#> 2 MCMDS   <tibble [292 × 6]>
#> 3 logFYLD <tibble [350 × 6]>
#> 4 logDYLD <tibble [348 × 6]>
```

The steps above set-us up almost all the way.


```r
# For fastest, lightest compute of accuracy, remove non-phenotyped from kinship

gids<-blups %>% 
     unnest(TrainingData) %$% unique(GID)
# dim(A) [1] 963 963

A<-A[gids,gids]
```

## Selection indices

Last thing: Let's include selection index weights. You can find an excellent, detailed, open-source chapter from Walsh & Lynch on Selection Index Theory by [**clicking here**](http://nitro.biosci.arizona.edu/zdownload/Volume2/Chapter23.pdf).

$$SI = WT_1 \times Trait_1 + \dots + WT_t \times Trait_t$$ Or in vector form:

$$SI = \boldsymbol{\hat{g}b}$$

$SI$ is the selection index, dimension $[n \times 1]$. $b$ is a $[t \times 1]$ vector of selection index "economic weights" designed to value each trait relative to its impact on the economic potential of changing the corresponding trait by one unit. Finally, $\boldsymbol{\hat{g}$ is a matrix $[n \times t]$ with the (in this case) **GEBV** for each trait on the columns.

`runCrossVal()` will accept a named vector of selection index weights where names must match the "Trait" variable in `blups` using the `SIwts=` argument and setting `selInd=TRUE`.

Here are example weights, I'll use. These are *not* to be taken as canonical. *Weights should be determined for each target population of environments and product profile!*


```r
# I chose to remove MCMDS 
## our preliminary analysis showed it to have ~0 heritability in this dataset
## initial test of cross-val. showed the models do not fit
SIwts<-c(DM=15,
         #MCMDS=-10,
         logFYLD=20,
         logDYLD=20)
SIwts
#>      DM logFYLD logDYLD 
#>      15      20      20
```

I'll run a meager 2 repetitions of 5-fold cross-validation, which means 10 predictions per trait overall. I've got a 16-core laptop so I can use `ncores=10` to do all 10 predictions per trait at the same time. `runCrossVal()` will process all four traits *and* compute the selection index accuracy at the end.

## Execute cross-validation


```r
starttime<-proc.time()[3]
standardCV<-runCrossVal(blups=blups %>% filter(Trait != "MCMDS"),
                        modelType="A",
                        selInd=TRUE,SIwts=SIwts,
                        grms=list(A=A),
                        nrepeats=2,nfolds=5,
                        gid="GID",seed=424242,
                        ncores=10)
#> Loading required package: rsample
#> Loading required package: furrr
#> Loading required package: future
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -150.932   13:5:19      0           0
#>     2      -150.587   13:5:19      0           0
#>     3      -150.456   13:5:19      0           0
#>     4      -150.431   13:5:19      0           0
#>     5      -150.429   13:5:19      0           0
#>     6      -150.429   13:5:19      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -109.584   13:5:19      0           0
#>     2      -109.57   13:5:20      1           0
#>     3      -109.562   13:5:20      1           0
#>     4      -109.56   13:5:20      1           0
#>     5      -109.559   13:5:20      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -115.829   13:5:20      0           0
#>     2      -115.829   13:5:20      0           0
#>     3      -115.828   13:5:20      0           0
#>     4      -115.828   13:5:20      0           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -153.247   13:5:19      0           0
#>     2      -153.244   13:5:19      0           0
#>     3      -153.243   13:5:20      1           0
#>     4      -153.243   13:5:20      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -108.226   13:5:20      0           0
#>     2      -108.147   13:5:20      0           0
#>     3      -108.101   13:5:20      0           0
#>     4      -108.087   13:5:20      0           0
#>     5      -108.085   13:5:20      0           0
#>     6      -108.085   13:5:20      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -117.592   13:5:21      0           0
#>     2      -117.537   13:5:21      0           0
#>     3      -117.513   13:5:21      0           0
#>     4      -117.509   13:5:21      0           0
#>     5      -117.508   13:5:21      0           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -150.198   13:5:20      0           0
#>     2      -149.363   13:5:20      0           0
#>     3      -148.987   13:5:20      0           0
#>     4      -148.881   13:5:20      0           0
#>     5      -148.865   13:5:20      0           0
#>     6      -148.863   13:5:20      0           0
#>     7      -148.862   13:5:21      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -106.107   13:5:21      0           0
#>     2      -105.581   13:5:21      0           0
#>     3      -105.152   13:5:21      0           0
#>     4      -104.92   13:5:21      0           0
#>     5      -104.852   13:5:21      0           0
#>     6      -104.832   13:5:21      0           0
#>     7      -104.827   13:5:21      0           0
#>     8      -104.825   13:5:21      0           0
#>     9      -104.825   13:5:21      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -118.481   13:5:22      1           0
#>     2      -118.255   13:5:22      1           0
#>     3      -118.106   13:5:22      1           0
#>     4      -118.047   13:5:22      1           0
#>     5      -118.035   13:5:22      1           0
#>     6      -118.032   13:5:22      1           0
#>     7      -118.032   13:5:22      1           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -144.958   13:5:21      0           0
#>     2      -144.946   13:5:21      0           0
#>     3      -144.94   13:5:21      0           0
#>     4      -144.939   13:5:21      0           0
#>     5      -144.939   13:5:21      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -107.241   13:5:21      0           0
#>     2      -107.24   13:5:21      0           0
#>     3      -107.24   13:5:21      0           0
#>     4      -107.24   13:5:21      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -114.776   13:5:22      1           0
#>     2      -114.775   13:5:22      1           0
#>     3      -114.775   13:5:22      1           0
#>     4      -114.775   13:5:22      1           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -150.502   13:5:21      0           0
#>     2      -150.404   13:5:21      0           0
#>     3      -150.354   13:5:21      0           0
#>     4      -150.339   13:5:21      0           0
#>     5      -150.336   13:5:21      0           0
#>     6      -150.336   13:5:21      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -112.48   13:5:22      0           0
#>     2      -112.42   13:5:22      0           0
#>     3      -112.38   13:5:22      0           0
#>     4      -112.364   13:5:22      0           0
#>     5      -112.36   13:5:22      0           0
#>     6      -112.358   13:5:22      0           0
#>     7      -112.358   13:5:22      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -118.347   13:5:22      0           0
#>     2      -118.041   13:5:23      1           0
#>     3      -117.869   13:5:23      1           0
#>     4      -117.803   13:5:23      1           0
#>     5      -117.787   13:5:23      1           0
#>     6      -117.784   13:5:23      1           0
#>     7      -117.783   13:5:23      1           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -150.226   13:5:22      0           0
#>     2      -149.466   13:5:22      0           0
#>     3      -149.138   13:5:22      0           0
#>     4      -149.063   13:5:22      0           0
#>     5      -149.056   13:5:22      0           0
#>     6      -149.055   13:5:22      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -111.205   13:5:22      0           0
#>     2      -111.2   13:5:22      0           0
#>     3      -111.196   13:5:22      0           0
#>     4      -111.193   13:5:22      0           0
#>     5      -111.193   13:5:23      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -115.15   13:5:23      0           0
#>     2      -115.132   13:5:23      0           0
#>     3      -115.119   13:5:23      0           0
#>     4      -115.114   13:5:23      0           0
#>     5      -115.113   13:5:23      0           0
#>     6      -115.112   13:5:23      0           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -150.983   13:5:22      0           0
#>     2      -150.511   13:5:22      0           0
#>     3      -150.265   13:5:22      0           0
#>     4      -150.179   13:5:22      0           0
#>     5      -150.162   13:5:22      0           0
#>     6      -150.158   13:5:23      1           0
#>     7      -150.157   13:5:23      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -109.264   13:5:23      0           0
#>     2      -109.264   13:5:23      0           0
#>     3      -109.264   13:5:23      0           0
#>     4      -109.263   13:5:23      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -116.271   13:5:23      0           0
#>     2      -116.238   13:5:23      0           0
#>     3      -116.225   13:5:23      0           0
#>     4      -116.223   13:5:24      1           0
#>     5      -116.223   13:5:24      1           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -146.729   13:5:23      0           0
#>     2      -146.707   13:5:23      0           0
#>     3      -146.695   13:5:23      0           0
#>     4      -146.691   13:5:23      0           0
#>     5      -146.691   13:5:23      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -105.14   13:5:23      0           0
#>     2      -105.116   13:5:23      0           0
#>     3      -105.101   13:5:23      0           0
#>     4      -105.095   13:5:23      0           0
#>     5      -105.095   13:5:24      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -116.469   13:5:24      0           0
#>     2      -116.439   13:5:24      0           0
#>     3      -116.428   13:5:24      0           0
#>     4      -116.426   13:5:24      0           0
#>     5      -116.426   13:5:24      0           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -146.167   13:5:23      0           0
#>     2      -145.784   13:5:23      0           0
#>     3      -145.645   13:5:23      0           0
#>     4      -145.618   13:5:23      0           0
#>     5      -145.616   13:5:24      1           0
#>     6      -145.616   13:5:24      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -108.335   13:5:24      0           0
#>     2      -108.255   13:5:24      0           0
#>     3      -108.205   13:5:24      0           0
#>     4      -108.187   13:5:24      0           0
#>     5      -108.184   13:5:24      0           0
#>     6      -108.184   13:5:24      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -115.606   13:5:24      0           0
#>     2      -115.563   13:5:24      0           0
#>     3      -115.541   13:5:24      0           0
#>     4      -115.535   13:5:25      1           0
#>     5      -115.534   13:5:25      1           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -152   13:5:24      0           0
#>     2      -151.698   13:5:24      0           0
#>     3      -151.579   13:5:24      0           0
#>     4      -151.555   13:5:24      0           0
#>     5      -151.554   13:5:24      0           0
#>     6      -151.553   13:5:24      0           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -107.98   13:5:24      0           0
#>     2      -107.972   13:5:24      0           0
#>     3      -107.968   13:5:24      0           0
#>     4      -107.967   13:5:25      1           0
#> [1] "GBLUP model complete - one trait"
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -119.501   13:5:25      0           0
#>     2      -119.452   13:5:25      0           0
#>     3      -119.431   13:5:25      0           0
#>     4      -119.426   13:5:25      0           0
#>     5      -119.425   13:5:25      0           0
#> [1] "GBLUP model complete - one trait"
#> [1] "Genomic predictions done for all traits in one repeat-fold"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
#> Joining, by = "GID"
timeelapsed<-proc.time()[3]-starttime; 
timeelapsed/60
#>   elapsed 
#> 0.1487667
```

Save the results


```r
saveRDS(standardCV,file = here::here("output","standardCV.rds"))
```

## Plot results


```r
standardCV %>% 
     unnest(accuracyEstOut) %>% 
     dplyr::select(repeats,id,predOf,Trait,Accuracy) %>% 
     ggplot(.,aes(x=Trait,y=Accuracy,fill=Trait)) + 
     geom_boxplot() + theme_bw()
```

<img src="07-kFoldCrossVal_files/figure-html/unnamed-chunk-7-1.png" width="672" />

This result is not what I would expect. SELIND should be similar to individual trait accuracies.

Best guess: SELIND requires the BLUPs for each trait to be observed, so only the clones with complete data will be included.

Your results will be different if you choose a different dataset, hopefully better.
