# Prepare training data

-   **Context and Purpose:** In this step, we do quality control, clean and format training data for further analysis.
-   **Upstream:** Section \@ref(download-training-data) - training data download
-   **Downstream:** pretty much everything
-   **Inputs:** "Raw" field trial data
-   **Expected outputs:** "Cleaned" field trial data


```r
suppressMessages(library(tidyverse)); 
suppressMessages(library(genomicMateSelectR))
```

## Read DB data

Load the phenotype and metadata downloads into R.

I built a function `readDBdata` that simply wraps around `read.csv`, reads and merges the metadata to the plot-basis data. The `metadataFile=` argument can be left NULL.


```r
dbdata<-readDBdata(phenotypeFile = here::here("data","phenotype.csv"),
                   metadataFile = here::here("data","metadata.csv"))
#> Joining, by = c("studyYear", "programDbId", "programName", "programDescription", "studyDbId", "studyName", "studyDescription", "studyDesign", "plotWidth", "plotLength", "fieldSize", "fieldTrialIsPlannedToBeGenotyped", "fieldTrialIsPlannedToCross", "plantingDate", "harvestDate", "locationDbId", "locationName")
```

**HINT:** At any point in the manual, if I reference or use a custom function in the `genomicMateSelectR`, I encourage you to check out the reference page for that function, e.g. [`readDBdata()`](https://wolfemd.github.io/genomicMateSelectR/reference/readDBdata.html). Or look at the code yourself by typing e.g. `readDBdata` at the R console or heading to the GitHub repo.

## Group and select trials to analyze

This step is present in my standard pipeline because I often bulk-download data for a breeding program and then sort through the trials-of-interest after the fact.
 
I think the better future way is to use trial-lists on BreedBase to communicate exactly which trials, years, locations are desired at download, maybe also which traits.

Unecessary for the example dataset chosen.

## Traits and Trait Abbreviations

Cassavabase downloads use very long column-names corresponding to the [full trait-ontology name](https://cropontology.org/ontology/CO_334). For convenience, I replace these names with abbreviations, documented here. For eventual upload of analysis results, names will need to be restored to ontology terms.

I also use this opportunity to subselect traits.

```r
traitabbrevs<-tribble(~TraitAbbrev,~TraitName,
        "CMD1S","cassava.mosaic.disease.severity.1.month.evaluation.CO_334.0000191",
        "CMD3S","cassava.mosaic.disease.severity.3.month.evaluation.CO_334.0000192",
        "CMD6S","cassava.mosaic.disease.severity.6.month.evaluation.CO_334.0000194",
        "DM","dry.matter.content.percentage.CO_334.0000092",
        "RTWT","fresh.storage.root.weight.per.plot.CO_334.0000012",
        "NOHAV","plant.stands.harvested.counting.CO_334.0000010")
traitabbrevs %>% gt::gt()#rmarkdown::paged_table()
```

```{=html}
<div id="csyvgbtwtr" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#csyvgbtwtr .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#csyvgbtwtr .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#csyvgbtwtr .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#csyvgbtwtr .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#csyvgbtwtr .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#csyvgbtwtr .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#csyvgbtwtr .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#csyvgbtwtr .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#csyvgbtwtr .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#csyvgbtwtr .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#csyvgbtwtr .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#csyvgbtwtr .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#csyvgbtwtr .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#csyvgbtwtr .gt_from_md > :first-child {
  margin-top: 0;
}

#csyvgbtwtr .gt_from_md > :last-child {
  margin-bottom: 0;
}

#csyvgbtwtr .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#csyvgbtwtr .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#csyvgbtwtr .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#csyvgbtwtr .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#csyvgbtwtr .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#csyvgbtwtr .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#csyvgbtwtr .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#csyvgbtwtr .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#csyvgbtwtr .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#csyvgbtwtr .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#csyvgbtwtr .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#csyvgbtwtr .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#csyvgbtwtr .gt_left {
  text-align: left;
}

#csyvgbtwtr .gt_center {
  text-align: center;
}

#csyvgbtwtr .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#csyvgbtwtr .gt_font_normal {
  font-weight: normal;
}

#csyvgbtwtr .gt_font_bold {
  font-weight: bold;
}

#csyvgbtwtr .gt_font_italic {
  font-style: italic;
}

#csyvgbtwtr .gt_super {
  font-size: 65%;
}

#csyvgbtwtr .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">TraitAbbrev</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">TraitName</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">CMD1S</td>
<td class="gt_row gt_left">cassava.mosaic.disease.severity.1.month.evaluation.CO_334.0000191</td></tr>
    <tr><td class="gt_row gt_left">CMD3S</td>
<td class="gt_row gt_left">cassava.mosaic.disease.severity.3.month.evaluation.CO_334.0000192</td></tr>
    <tr><td class="gt_row gt_left">CMD6S</td>
<td class="gt_row gt_left">cassava.mosaic.disease.severity.6.month.evaluation.CO_334.0000194</td></tr>
    <tr><td class="gt_row gt_left">DM</td>
<td class="gt_row gt_left">dry.matter.content.percentage.CO_334.0000092</td></tr>
    <tr><td class="gt_row gt_left">RTWT</td>
<td class="gt_row gt_left">fresh.storage.root.weight.per.plot.CO_334.0000012</td></tr>
    <tr><td class="gt_row gt_left">NOHAV</td>
<td class="gt_row gt_left">plant.stands.harvested.counting.CO_334.0000010</td></tr>
  </tbody>
  
  
</table>
</div>
```
Run function `renameAndSelectCols()` to rename columns and remove unselected traits.

```r
dbdata<-renameAndSelectCols(traitabbrevs,
                            indata=dbdata,
                            customColsToKeep = c("observationUnitName"))
#> Joining, by = "TraitName"
```

## QC Trait Values

At this point in the pipeline, we should check the all trait values are in allowable ranges. 
Different ways to approach this. Feel free to make some plots of your data!

The database also has mechanisms to ensure trait values are only within allowable ranges. 
Nevertheless, as a habit, I have an simple _ad hoc_ approach to this:


```r
# comment out the traits not present in this dataset
dbdata<-dbdata %>% 
  dplyr::mutate(CMD1S=ifelse(CMD1S<1 | CMD1S>5,NA,CMD1S),
         CMD3S=ifelse(CMD3S<1 | CMD3S>5,NA,CMD3S),
         #CMD6S=ifelse(CMD6S<1 | CMD6S>5,NA,CMD6S), 
         #CMD9S=ifelse(CMD9S<1 | CMD9S>5,NA,CMD9S),
         # CGM=ifelse(CGM<1 | CGM>5,NA,CGM),
         # CGMS1=ifelse(CGMS1<1 | CGMS1>5,NA,CGMS1),
         # CGMS2=ifelse(CGMS2<1 | CGMS2>5,NA,CGMS2),
         DM=ifelse(DM>100 | DM<=0,NA,DM),
         RTWT=ifelse(RTWT==0 | NOHAV==0 | is.na(NOHAV),NA,RTWT),
         # SHTWT=ifelse(SHTWT==0 | NOHAV==0 | is.na(NOHAV),NA,SHTWT),
         # RTNO=ifelse(RTNO==0 | NOHAV==0 | is.na(NOHAV),NA,RTNO),
         NOHAV=ifelse(NOHAV==0,NA,NOHAV),
         NOHAV=ifelse(NOHAV>42,NA,NOHAV)
         #RTNO=ifelse(!RTNO %in% 1:10000,NA,RTNO)
         )
```

## Post-QC: composite traits

Now that component traits are QC'd, it's time to compute any composite traits.

By composite traits, I mean traits computed from combinations of other traits.

Examples for cassava: season-wide mean disease severity, harvest index, and fresh root yield.

## Season-wide mean disease severity


```r
# [NEW AS OF APRIL 2021]
## VERSION with vs. without CBSD
## Impervious to particular timepoints between 1, 3, 6 and 9 scores

# Without CBSD (West Africa)
dbdata<-dbdata %>% 
  mutate(MCMDS=rowMeans(.[,colnames(.) %in% c("CMD1S","CMD3S","CMD6S","CMD9S")], na.rm = T)) %>% 
  select(-any_of(c("CMD1S","CMD3S","CMD6S","CMD9S")))

# With CBSD (East Africa)
# dbdata<-dbdata %>% 
#   mutate(MCMDS=rowMeans(.[,colnames(.) %in% c("CMD1S","CMD3S","CMD6S","CMD9S")], na.rm = T),
#          MCBSDS=rowMeans(.[,colnames(.) %in% c("CBSD1S","CBSD3S","CBSD6S","CBSD9S")], na.rm = T)) %>% 
#   select(-any_of(c("CMD1S","CMD3S","CMD6S","CMD9S","CBSD1S","CBSD3S","CBSD6S","CBSD9S")))
```

## Fresh root yield (FYLD)

**RTWT** (fresh root weight per plot in kg) --> **FYLD** (fresh root yield in tons per hectare)

$$FYLD = \frac{RTWT_{kg / plot}}{MaxHarvestedPlantsPerPlot \times PlantSpacing}\times10 $$

**NOTE:** *MaxHarvestedPlantsPerPlot* in formula above is to distinguish from the *plantsPerPlot* meta-data field, in case that a net-plot harvest is used. In other words, the value should be the total number of plants intended for harvest in a plot.

*PlantSpacing* is the area in $m^2$ per plant. 

In the example trial data, we the `plantsPerPlot` meta-data field is empty. Luckily, since there are only two trials, we make a quick summary of the NOHAV data, to determine the correct values. 

**RECOMMEND INPUTING plantsPerPlot meta-data to cassavabase**


```r
dbdata %>% count(studyYear,studyName,studyDesign,plotWidth,plotLength)
#> # A tibble: 9 × 6
#>   studyYear studyName studyDesign plotWidth plotLength     n
#>       <int> <chr>     <chr>           <int>      <dbl> <int>
#> 1      2019 19.GS.C1… Alpha               4        4     125
#> 2      2019 19.GS.C2… Alpha               4        4      68
#> 3      2019 19.GS.C2… Alpha               4        4      72
#> 4      2019 19.GS.C2… RCBD                4        4      66
#> 5      2019 19.GS.C2… Alpha               4        4      72
#> 6      2019 19.GS.C4… Alpha               2        4     270
#> 7      2019 19.GS.C4… Alpha               3        2.5   273
#> 8      2019 19.GS.C4… RCBD                1        2.5   777
#> 9      2019 19geneti… Augmented           1        2.5   810
```

So the GS.C1 trial has 2.5 $m^2$ plots, the GeneticGain trial has 8 $m^2$. 

A quick density plot reveals that the GeneticGain trial was likely planted with 10 plants/plot, and the GS.C1.CE 5 plants/plot. (*Disclosure:* I know this is true from experience.)

```r
dbdata %>% 
     ggplot(.,aes(x=NOHAV, fill=studyName)) + geom_density(alpha=0.75)
#> Warning: Removed 921 rows containing non-finite values
#> (stat_density).
```

<img src="02-training_data_QC_files/figure-html/unnamed-chunk-8-1.png" width="672" />


```r
dbdata %<>% 
     # plot area in meters squared
     mutate(plotArea=plotWidth*plotLength,
            # Number of plants per plot
            plantsPerPlot=ifelse(studyName=="13geneticgainIB",10,5))
```


```r
dbdata %<>% 
     mutate(PlantSpacing=plotArea/plantsPerPlot,
            FYLD=RTWT/(plantsPerPlot*PlantSpacing)*10)
dbdata %>% ggplot(.,aes(x=FYLD,fill=studyName)) + geom_density(alpha=0.75)
#> Warning: Removed 943 rows containing non-finite values
#> (stat_density).
```

<img src="02-training_data_QC_files/figure-html/unnamed-chunk-10-1.png" width="672" />

Additional things to compute:

1. log-transform yield traits: this is a habit based on experience. Linear mixed-models should have normally distributed homoskedastic residuals, if they don't log-transform the response variable often helps. For FYLD and related traits, I always log-transform. 


```r
# I log transform yield traits 
# to satisfy homoskedastic residuals assumption 
# of linear mixed models
dbdata %<>% 
     mutate(DYLD=FYLD*(DM/100),
            logFYLD=log(FYLD),
            logDYLD=log(DYLD),
            PropNOHAV=NOHAV/plantsPerPlot) 
# remove non transformed / per-plot (instead of per area) traits
dbdata %<>% select(-RTWT,-FYLD,-DYLD)
dbdata %>% ggplot(.,aes(x=logFYLD,fill=studyName)) + geom_density(alpha=0.75)
#> Warning: Removed 943 rows containing non-finite values
#> (stat_density).
```

<img src="02-training_data_QC_files/figure-html/unnamed-chunk-11-1.png" width="672" />

Debatable whether this is better. Let's not dwell on it. Onward!

## Check genotype-to-phenotype matches

In this step, in the past (see [here](https://wolfemd.github.io/IITA_2021GS/01-cleanTPdata.html#[User_input]_Assign_genos_to_phenos) for an example), I have had to use multiple external flat files and custom code, built over time, to figure out which DNA samples, if any, match eatch field plot. 

With both genotypes and phenotypes derived from the database, this should not be an issue for this tutorial. 


## Check / "detect" experimental designs? 

In this step, in the past, I have not been certain of the experimental designs of the trials I had downloaded. I was also not certain how the designs were represented in the column-names. For this reason, I developed another _ad hoc_ custom code to "detect" the designs. I built the `genomicMateSelectR` function `detectExptDesigns()`. See an example [here](https://wolfemd.github.io/IITA_2021GS/01-cleanTPdata.html#Detect_experimental_designs).

Should not be necessary for the example dataset. However, let's check:


```r
dbdata %>% count(studyName,studyDesign, numberBlocks,replicate, blockNumber, entryType)
#> # A tibble: 224 × 7
#>    studyName  studyDesign numberBlocks replicate blockNumber
#>    <chr>      <chr>              <int>     <int>       <int>
#>  1 19.GS.C1.… Alpha                 NA         1           1
#>  2 19.GS.C1.… Alpha                 NA         1           1
#>  3 19.GS.C1.… Alpha                 NA         2           2
#>  4 19.GS.C1.… Alpha                 NA         2           2
#>  5 19.GS.C1.… Alpha                 NA         3           3
#>  6 19.GS.C1.… Alpha                 NA         3           3
#>  7 19.GS.C2.… Alpha                  6         1           1
#>  8 19.GS.C2.… Alpha                  6         1           1
#>  9 19.GS.C2.… Alpha                  6         1           2
#> 10 19.GS.C2.… Alpha                  6         1           2
#> # … with 214 more rows, and 2 more variables:
#> #   entryType <chr>, n <int>
```

## Checklist and KPI

- [ ]
