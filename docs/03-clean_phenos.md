# Prepare phenotype data


-   **Context and Purpose:** In this step, we do quality control, clean and format training data for further analysis.
-   **Upstream:** Section \@ref(download-training-data) - training data download
-   **Downstream:** pretty much everything
-   **Inputs:** "Raw" field trial data
-   **Expected outputs:** "Cleaned" field trial data

## Read DB data

Load the phenotype and metadata downloads into R.

I built a function `readDBdata` that simply wraps around `read.csv`, reads and merges the metadata to the plot-basis data. The `metadataFile=` argument can be left NULL.


```r
dbdata<-readDBdata(phenotypeFile = here::here("data","phenotype.csv"),
                   metadataFile = here::here("data","metadata.csv"))
#> Joining, by = c("studyYear", "programDbId", "programName", "programDescription", "studyDbId", "studyName", "studyDescription", "studyDesign", "plotWidth", "plotLength", "fieldSize", "fieldTrialIsPlannedToBeGenotyped", "fieldTrialIsPlannedToCross", "plantingDate", "harvestDate", "locationDbId", "locationName")
```

**HINT:** At any point in the manual, if I reference or use a custom function in the `genomicMateSelectR`, I encourage you to check out the reference page for that function, e.g. [`readDBdata()`](https://wolfemd.github.io/genomicMateSelectR/reference/readDBdata.html). Or look at the code yourself by typing e.g. `readDBdata` at the R console or heading to the GitHub repo.

## Check experimental designs {#detect_designs}

**Checklist:** Are the data plot-basis, plant-basis or a mixture? If plant-basis data are present, should they be converted to plot-basis for further analysis?


```r
dbdata %>% count(observationLevel)
#>   observationLevel    n
#> 1             plot 2533
# table(dbdata$observationLevel)
```
Only plot-basis in this case.

**Checklist:** What experimental designs are present? How are they represented by the variables in the dataset? Are all designs consistent with your expectations, for example relative to the reported “`trialType`”, “`studyName`” and/or “`studyDesign`”?

In this step, in the past, I have not been certain of the experimental designs of the trials I had downloaded. I was also not certain how the designs were represented in the column-names. For this reason, I developed an *ad hoc* custom code to "detect" the designs. I built the `genomicMateSelectR` function `detectExptDesigns()`. See an example [here](https://wolfemd.github.io/IITA_2021GS/01-cleanTPdata.html#Detect_experimental_designs).

**RECOMMENDATION:** Each analyst needs to use exploratory data anlaysis, making summary statistics and plots as necessary to determine how the data should be modelled downstream. If there are missing or incorrectly represented trial design variables, ***get it corrected on the database (contact breeding program data manager, as necessary)***.

Because I have a small example dataset, it is possible to look at 9 trials and evaluate. 

Often, _many_ more trials are part of a genomic prediction. This is why it is essential that trial designs be consistent, and clear to the analyst. You may need to derive a strategy similar to the `detectExptDesigns()` function to semi-automate the process.


```r
library(gt)
dbdata %>% 
     count(studyName,trialType, studyDesign, numberBlocks,numberReps,entryType) %>% 
     spread(entryType,n) %>% 
     gt()  %>% 
     tab_options(table.font.size = pct(75))
```

```{=html}
<div id="bzcneabwdu" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#bzcneabwdu .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 75%;
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

#bzcneabwdu .gt_heading {
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

#bzcneabwdu .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#bzcneabwdu .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#bzcneabwdu .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bzcneabwdu .gt_col_headings {
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

#bzcneabwdu .gt_col_heading {
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

#bzcneabwdu .gt_column_spanner_outer {
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

#bzcneabwdu .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bzcneabwdu .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bzcneabwdu .gt_column_spanner {
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

#bzcneabwdu .gt_group_heading {
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

#bzcneabwdu .gt_empty_group_heading {
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

#bzcneabwdu .gt_from_md > :first-child {
  margin-top: 0;
}

#bzcneabwdu .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bzcneabwdu .gt_row {
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

#bzcneabwdu .gt_stub {
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

#bzcneabwdu .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bzcneabwdu .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#bzcneabwdu .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bzcneabwdu .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bzcneabwdu .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bzcneabwdu .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bzcneabwdu .gt_footnotes {
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

#bzcneabwdu .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#bzcneabwdu .gt_sourcenotes {
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

#bzcneabwdu .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#bzcneabwdu .gt_left {
  text-align: left;
}

#bzcneabwdu .gt_center {
  text-align: center;
}

#bzcneabwdu .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bzcneabwdu .gt_font_normal {
  font-weight: normal;
}

#bzcneabwdu .gt_font_bold {
  font-weight: bold;
}

#bzcneabwdu .gt_font_italic {
  font-style: italic;
}

#bzcneabwdu .gt_super {
  font-size: 65%;
}

#bzcneabwdu .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">trialType</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyDesign</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">numberBlocks</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">numberReps</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">check</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">test</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_left">NA</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">NA</td>
<td class="gt_row gt_right">NA</td>
<td class="gt_row gt_right">15</td>
<td class="gt_row gt_right">110</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_left">Uniform Yield Trial</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">58</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_left">Uniform Yield Trial</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">62</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_left">Uniform Yield Trial</td>
<td class="gt_row gt_left">RCBD</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">56</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_left">Uniform Yield Trial</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">62</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_left">Preliminary Yield Trial</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">30</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">12</td>
<td class="gt_row gt_right">258</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_left">Preliminary Yield Trial</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">28</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">31</td>
<td class="gt_row gt_right">242</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_left">Clonal Evaluation</td>
<td class="gt_row gt_left">RCBD</td>
<td class="gt_row gt_right">42</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">132</td>
<td class="gt_row gt_right">645</td></tr>
    <tr><td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_left">genetic_gain_trial</td>
<td class="gt_row gt_left">Augmented</td>
<td class="gt_row gt_right">11</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">NA</td>
<td class="gt_row gt_right">810</td></tr>
  </tbody>
  
  
</table>
</div>
```

Summary table above shows:

1) **trialType** and **studyDesign** cannot be 100% relied upon, at least not here. 
2) The only trial actually listed as having `studyDesign=="Augmented"` does not have "check" vs. "test" distinguished in the "entryType".
3) A `trialType=="Clonal Evaluation"` with `studyDesign=="RCBD"` but actually only 1 replication.

Next, I'll check if the **replicate** and **blockNumber** columns reliably distinguish complete and incomplete blocks in the data.


```r
dbdata %>% 
     group_by(studyName) %>% 
     summarize(N_replicate=length(unique(replicate)),
               N_blockNumber=length(unique(blockNumber))) %>% 
     gt() %>% tab_options(table.font.size = pct(75))
```

```{=html}
<div id="thxssybctx" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#thxssybctx .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 75%;
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

#thxssybctx .gt_heading {
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

#thxssybctx .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#thxssybctx .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#thxssybctx .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#thxssybctx .gt_col_headings {
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

#thxssybctx .gt_col_heading {
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

#thxssybctx .gt_column_spanner_outer {
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

#thxssybctx .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#thxssybctx .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#thxssybctx .gt_column_spanner {
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

#thxssybctx .gt_group_heading {
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

#thxssybctx .gt_empty_group_heading {
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

#thxssybctx .gt_from_md > :first-child {
  margin-top: 0;
}

#thxssybctx .gt_from_md > :last-child {
  margin-bottom: 0;
}

#thxssybctx .gt_row {
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

#thxssybctx .gt_stub {
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

#thxssybctx .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#thxssybctx .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#thxssybctx .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#thxssybctx .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#thxssybctx .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#thxssybctx .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#thxssybctx .gt_footnotes {
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

#thxssybctx .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#thxssybctx .gt_sourcenotes {
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

#thxssybctx .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#thxssybctx .gt_left {
  text-align: left;
}

#thxssybctx .gt_center {
  text-align: center;
}

#thxssybctx .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#thxssybctx .gt_font_normal {
  font-weight: normal;
}

#thxssybctx .gt_font_bold {
  font-weight: bold;
}

#thxssybctx .gt_font_italic {
  font-style: italic;
}

#thxssybctx .gt_super {
  font-size: 65%;
}

#thxssybctx .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">N_replicate</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">N_blockNumber</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_right">3</td>
<td class="gt_row gt_right">3</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">30</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">28</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">42</td></tr>
    <tr><td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">11</td></tr>
  </tbody>
  
  
</table>
</div>
```

Here, I notice that except 1 trial (**19.GS.C1.C2.C3.AYT.42.UB**) has the same number of reps and blocks. 

The question is, are complete replications of the experiment indicated by **replicate** and incomplete sub-blocks represented by **blockNumber**


```r
dbdata %>% 
     group_by(studyName) %>% 
     summarize(N_replicate=length(unique(replicate)),
               N_blockNumber=length(unique(blockNumber)),
               doRepsEqualBlocks=all(replicate==blockNumber)) %>% 
     gt() %>% tab_options(table.font.size = pct(75))
```

```{=html}
<div id="lamhrhupdh" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#lamhrhupdh .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 75%;
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

#lamhrhupdh .gt_heading {
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

#lamhrhupdh .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#lamhrhupdh .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#lamhrhupdh .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#lamhrhupdh .gt_col_headings {
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

#lamhrhupdh .gt_col_heading {
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

#lamhrhupdh .gt_column_spanner_outer {
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

#lamhrhupdh .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#lamhrhupdh .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#lamhrhupdh .gt_column_spanner {
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

#lamhrhupdh .gt_group_heading {
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

#lamhrhupdh .gt_empty_group_heading {
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

#lamhrhupdh .gt_from_md > :first-child {
  margin-top: 0;
}

#lamhrhupdh .gt_from_md > :last-child {
  margin-bottom: 0;
}

#lamhrhupdh .gt_row {
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

#lamhrhupdh .gt_stub {
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

#lamhrhupdh .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#lamhrhupdh .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#lamhrhupdh .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#lamhrhupdh .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#lamhrhupdh .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#lamhrhupdh .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#lamhrhupdh .gt_footnotes {
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

#lamhrhupdh .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#lamhrhupdh .gt_sourcenotes {
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

#lamhrhupdh .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#lamhrhupdh .gt_left {
  text-align: left;
}

#lamhrhupdh .gt_center {
  text-align: center;
}

#lamhrhupdh .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#lamhrhupdh .gt_font_normal {
  font-weight: normal;
}

#lamhrhupdh .gt_font_bold {
  font-weight: bold;
}

#lamhrhupdh .gt_font_italic {
  font-style: italic;
}

#lamhrhupdh .gt_super {
  font-size: 65%;
}

#lamhrhupdh .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">N_replicate</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">N_blockNumber</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">doRepsEqualBlocks</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_right">3</td>
<td class="gt_row gt_right">3</td>
<td class="gt_row gt_center">TRUE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">30</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">28</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">42</td>
<td class="gt_row gt_center">FALSE</td></tr>
    <tr><td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">11</td>
<td class="gt_row gt_center">FALSE</td></tr>
  </tbody>
  
  
</table>
</div>
```

So for 1 trial, there are 3 complete blocks, no sub-blocks. For 6 trials, there are 2 complete replications and nested sub-blocks represented by the **blockNumber** variable. For 2 trials, there are only incomplete blocks. 


Next, I decided to check that the **replicate** column definitely means complete blocks. The below might look a bit complicated, but I basically merge two summaries: (1) he overall number of accessions per trial, and (2) the average number of accessions per replicate per trial. 

```r
# the overall number of accessions per trial
dbdata %>% 
     group_by(studyName) %>% 
     summarize(N_accession=length(unique(germplasmName))) %>% 
     # the average number of accessions per replicate per trial
     left_join(dbdata %>% 
                    group_by(studyName,replicate) %>% 
                    summarize(N_accession=length(unique(germplasmName))) %>% 
                    group_by(studyName) %>% 
                    summarize(avgAccessionsPerReplicate=ceiling(mean(N_accession)))) %>% 
     gt() %>% tab_options(table.font.size = pct(75))
#> `summarise()` has grouped output by 'studyName'. You can override using the `.groups` argument.
#> Joining, by = "studyName"
```

```{=html}
<div id="ozlethzisw" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#ozlethzisw .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 75%;
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

#ozlethzisw .gt_heading {
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

#ozlethzisw .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ozlethzisw .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ozlethzisw .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ozlethzisw .gt_col_headings {
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

#ozlethzisw .gt_col_heading {
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

#ozlethzisw .gt_column_spanner_outer {
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

#ozlethzisw .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ozlethzisw .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ozlethzisw .gt_column_spanner {
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

#ozlethzisw .gt_group_heading {
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

#ozlethzisw .gt_empty_group_heading {
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

#ozlethzisw .gt_from_md > :first-child {
  margin-top: 0;
}

#ozlethzisw .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ozlethzisw .gt_row {
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

#ozlethzisw .gt_stub {
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

#ozlethzisw .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ozlethzisw .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#ozlethzisw .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ozlethzisw .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ozlethzisw .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ozlethzisw .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ozlethzisw .gt_footnotes {
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

#ozlethzisw .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#ozlethzisw .gt_sourcenotes {
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

#ozlethzisw .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#ozlethzisw .gt_left {
  text-align: left;
}

#ozlethzisw .gt_center {
  text-align: center;
}

#ozlethzisw .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ozlethzisw .gt_font_normal {
  font-weight: normal;
}

#ozlethzisw .gt_font_bold {
  font-weight: bold;
}

#ozlethzisw .gt_font_italic {
  font-style: italic;
}

#ozlethzisw .gt_super {
  font-size: 65%;
}

#ozlethzisw .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">N_accession</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">avgAccessionsPerReplicate</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_right">42</td>
<td class="gt_row gt_right">42</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_right">35</td>
<td class="gt_row gt_right">34</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_right">36</td>
<td class="gt_row gt_right">36</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_right">36</td>
<td class="gt_row gt_right">33</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_right">36</td>
<td class="gt_row gt_right">36</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_right">135</td>
<td class="gt_row gt_right">135</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_right">129</td>
<td class="gt_row gt_right">127</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_right">657</td>
<td class="gt_row gt_right">657</td></tr>
    <tr><td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_right">782</td>
<td class="gt_row gt_right">782</td></tr>
  </tbody>
  
  
</table>
</div>
```
The numbers are very similar for all trials, indicating complete blocks. 

One more: look at the min, mean and max number of accessions per **blockNumber**.

```r
# the overall number of accessions per trial
dbdata %>% 
     group_by(studyName) %>% 
     summarize(N_accession=length(unique(germplasmName))) %>% 
     left_join(dbdata %>% 
     group_by(studyName,replicate,blockNumber) %>% 
     summarize(N_accession=length(unique(germplasmName))) %>% ungroup() %>% 
     group_by(studyName) %>% 
     summarize(minAccessionsPerBlock=ceiling(min(N_accession)),
               avgAccessionsPerBlock=ceiling(mean(N_accession)),
               maxAccessionsPerBlock=ceiling(max(N_accession)))) %>% 
     gt() %>% tab_options(table.font.size = pct(60))
#> `summarise()` has grouped output by 'studyName', 'replicate'. You can override using the `.groups` argument.
#> Joining, by = "studyName"
```

```{=html}
<div id="mriozpdmst" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#mriozpdmst .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 60%;
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

#mriozpdmst .gt_heading {
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

#mriozpdmst .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#mriozpdmst .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#mriozpdmst .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mriozpdmst .gt_col_headings {
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

#mriozpdmst .gt_col_heading {
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

#mriozpdmst .gt_column_spanner_outer {
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

#mriozpdmst .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#mriozpdmst .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#mriozpdmst .gt_column_spanner {
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

#mriozpdmst .gt_group_heading {
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

#mriozpdmst .gt_empty_group_heading {
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

#mriozpdmst .gt_from_md > :first-child {
  margin-top: 0;
}

#mriozpdmst .gt_from_md > :last-child {
  margin-bottom: 0;
}

#mriozpdmst .gt_row {
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

#mriozpdmst .gt_stub {
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

#mriozpdmst .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mriozpdmst .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#mriozpdmst .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mriozpdmst .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#mriozpdmst .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#mriozpdmst .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mriozpdmst .gt_footnotes {
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

#mriozpdmst .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#mriozpdmst .gt_sourcenotes {
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

#mriozpdmst .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#mriozpdmst .gt_left {
  text-align: left;
}

#mriozpdmst .gt_center {
  text-align: center;
}

#mriozpdmst .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#mriozpdmst .gt_font_normal {
  font-weight: normal;
}

#mriozpdmst .gt_font_bold {
  font-weight: bold;
}

#mriozpdmst .gt_font_italic {
  font-style: italic;
}

#mriozpdmst .gt_super {
  font-size: 65%;
}

#mriozpdmst .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">N_accession</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">minAccessionsPerBlock</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">avgAccessionsPerBlock</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">maxAccessionsPerBlock</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_right">42</td>
<td class="gt_row gt_right">41</td>
<td class="gt_row gt_right">42</td>
<td class="gt_row gt_right">42</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_right">35</td>
<td class="gt_row gt_right">11</td>
<td class="gt_row gt_right">12</td>
<td class="gt_row gt_right">12</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_right">36</td>
<td class="gt_row gt_right">12</td>
<td class="gt_row gt_right">12</td>
<td class="gt_row gt_right">12</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_right">36</td>
<td class="gt_row gt_right">9</td>
<td class="gt_row gt_right">11</td>
<td class="gt_row gt_right">12</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_right">36</td>
<td class="gt_row gt_right">12</td>
<td class="gt_row gt_right">12</td>
<td class="gt_row gt_right">12</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_right">135</td>
<td class="gt_row gt_right">9</td>
<td class="gt_row gt_right">9</td>
<td class="gt_row gt_right">9</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_right">129</td>
<td class="gt_row gt_right">8</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">10</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_right">657</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">19</td>
<td class="gt_row gt_right">20</td></tr>
    <tr><td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_right">782</td>
<td class="gt_row gt_right">39</td>
<td class="gt_row gt_right">73</td>
<td class="gt_row gt_right">80</td></tr>
  </tbody>
  
  
</table>
</div>
```
From this, you can see that except for `studyName=="19.GS.C1.C2.C3.AYT.42.UB"` the sub-blocks represented by **blockNumber** have only subsets of the total number of accessions in the trial, as expected. 

Further, except for `studyName=="19geneticgainUB"` all trials have pretty consistently sized sub-blocks.

Now I will _ad hoc_ create two variables (**CompleteBlocks** and **IncompleteBlocks**), indicating (`TRUE/FALSE`) whether to model using the **replicate** and/or **blockNumber** variable. 

I _also_ like to create _explicitly nested_ design variables (**yearInLoc**, **trialInLocYr**, **repInTrial**, **blockInRep**). 


```r
dbdata %<>% 
     group_by(studyName) %>% 
     summarize(N_replicate=length(unique(replicate)),
               N_blockNumber=length(unique(blockNumber)),
               doRepsEqualBlocks=all(replicate==blockNumber)) %>% 
     ungroup() %>% 
     mutate(CompleteBlocks=ifelse(N_replicate>1,TRUE,FALSE),
            IncompleteBlocks=ifelse(N_blockNumber>1 & !doRepsEqualBlocks,TRUE,FALSE)) %>% 
     left_join(dbdata) %>% 
     mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
            trialInLocYr=paste0(yearInLoc,"_",studyName),
            repInTrial=paste0(trialInLocYr,"_",replicate),
            blockInRep=paste0(repInTrial,"_",blockNumber))
#> Joining, by = "studyName"
```

Just to check:

```r
dbdata %>% 
     count(studyName,CompleteBlocks,IncompleteBlocks) %>% 
     left_join(dbdata %>% 
                    group_by(studyName) %>% 
                    summarize(nRepInTrial=length(unique(repInTrial)),
                              nBlockInRep=length(unique(blockInRep)))) %>% 
     gt() %>% tab_options(table.font.size = pct(67))
#> Joining, by = "studyName"
```

```{=html}
<div id="nrhicrtevs" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#nrhicrtevs .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 67%;
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

#nrhicrtevs .gt_heading {
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

#nrhicrtevs .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#nrhicrtevs .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#nrhicrtevs .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nrhicrtevs .gt_col_headings {
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

#nrhicrtevs .gt_col_heading {
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

#nrhicrtevs .gt_column_spanner_outer {
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

#nrhicrtevs .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#nrhicrtevs .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#nrhicrtevs .gt_column_spanner {
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

#nrhicrtevs .gt_group_heading {
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

#nrhicrtevs .gt_empty_group_heading {
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

#nrhicrtevs .gt_from_md > :first-child {
  margin-top: 0;
}

#nrhicrtevs .gt_from_md > :last-child {
  margin-bottom: 0;
}

#nrhicrtevs .gt_row {
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

#nrhicrtevs .gt_stub {
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

#nrhicrtevs .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#nrhicrtevs .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#nrhicrtevs .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#nrhicrtevs .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#nrhicrtevs .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#nrhicrtevs .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nrhicrtevs .gt_footnotes {
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

#nrhicrtevs .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#nrhicrtevs .gt_sourcenotes {
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

#nrhicrtevs .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#nrhicrtevs .gt_left {
  text-align: left;
}

#nrhicrtevs .gt_center {
  text-align: center;
}

#nrhicrtevs .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#nrhicrtevs .gt_font_normal {
  font-weight: normal;
}

#nrhicrtevs .gt_font_bold {
  font-weight: bold;
}

#nrhicrtevs .gt_font_italic {
  font-style: italic;
}

#nrhicrtevs .gt_super {
  font-size: 65%;
}

#nrhicrtevs .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">CompleteBlocks</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">IncompleteBlocks</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">n</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">nRepInTrial</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">nBlockInRep</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_right">125</td>
<td class="gt_row gt_right">3</td>
<td class="gt_row gt_right">3</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">68</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">72</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">66</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">72</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">6</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">270</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">30</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">273</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">28</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">777</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">42</td></tr>
    <tr><td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_center">TRUE</td>
<td class="gt_row gt_right">810</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">11</td></tr>
  </tbody>
  
  
</table>
</div>
```
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
traitabbrevs %>% gt()#rmarkdown::paged_table()
```

```{=html}
<div id="hvhscscudt" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#hvhscscudt .gt_table {
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

#hvhscscudt .gt_heading {
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

#hvhscscudt .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#hvhscscudt .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#hvhscscudt .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hvhscscudt .gt_col_headings {
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

#hvhscscudt .gt_col_heading {
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

#hvhscscudt .gt_column_spanner_outer {
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

#hvhscscudt .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#hvhscscudt .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#hvhscscudt .gt_column_spanner {
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

#hvhscscudt .gt_group_heading {
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

#hvhscscudt .gt_empty_group_heading {
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

#hvhscscudt .gt_from_md > :first-child {
  margin-top: 0;
}

#hvhscscudt .gt_from_md > :last-child {
  margin-bottom: 0;
}

#hvhscscudt .gt_row {
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

#hvhscscudt .gt_stub {
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

#hvhscscudt .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hvhscscudt .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#hvhscscudt .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hvhscscudt .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#hvhscscudt .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#hvhscscudt .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hvhscscudt .gt_footnotes {
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

#hvhscscudt .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#hvhscscudt .gt_sourcenotes {
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

#hvhscscudt .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#hvhscscudt .gt_left {
  text-align: left;
}

#hvhscscudt .gt_center {
  text-align: center;
}

#hvhscscudt .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#hvhscscudt .gt_font_normal {
  font-weight: normal;
}

#hvhscscudt .gt_font_bold {
  font-weight: bold;
}

#hvhscscudt .gt_font_italic {
  font-style: italic;
}

#hvhscscudt .gt_super {
  font-size: 65%;
}

#hvhscscudt .gt_footnote_marks {
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
                            customColsToKeep = c("observationUnitName",
                                                 "CompleteBlocks",
                                                 "IncompleteBlocks",
                                                 "yearInLoc",
                                                 "trialInLocYr",
                                                 "repInTrial","blockInRep"))
#> Joining, by = "TraitName"
```

## QC Trait Values

At this point in the pipeline, we should check the all trait values are in allowable ranges. Different ways to approach this. Feel free to make some plots of your data!

The database also has mechanisms to ensure trait values are only within allowable ranges. 

Nevertheless, as a habit, I have an simple *ad hoc* approach to this:


```r
# comment out the traits not present in this dataset
dbdata<-dbdata %>% 
     dplyr::mutate(CMD1S=ifelse(CMD1S<1 | CMD1S>5,NA,CMD1S),
                   CMD3S=ifelse(CMD3S<1 | CMD3S>5,NA,CMD3S),
                   # CMD6S=ifelse(CMD6S<1 | CMD6S>5,NA,CMD6S), 
                   # CMD9S=ifelse(CMD9S<1 | CMD9S>5,NA,CMD9S),
                   # CGM=ifelse(CGM<1 | CGM>5,NA,CGM),
                   # CGMS1=ifelse(CGMS1<1 | CGMS1>5,NA,CGMS1),
                   # CGMS2=ifelse(CGMS2<1 | CGMS2>5,NA,CGMS2),
                   DM=ifelse(DM>100 | DM<=0,NA,DM),
                   RTWT=ifelse(RTWT==0 | NOHAV==0 | is.na(NOHAV),NA,RTWT),
                   # SHTWT=ifelse(SHTWT==0 | NOHAV==0 | is.na(NOHAV),NA,SHTWT),
                   # RTNO=ifelse(RTNO==0 | NOHAV==0 | is.na(NOHAV),NA,RTNO),
                   NOHAV=ifelse(NOHAV==0,NA,NOHAV),
                   NOHAV=ifelse(NOHAV>42,NA,NOHAV)
                   # RTNO=ifelse(!RTNO %in% 1:10000,NA,RTNO)
     )
```

## Post-QC: composite traits

Now that component traits are QC'd, it's time to compute any composite traits.

By composite traits, I mean traits computed from combinations of other traits.

Examples for cassava: season-wide mean disease severity, harvest index, and fresh root yield.

### Season-wide mean disease severity


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

### Fresh root yield (FYLD)

**RTWT** (fresh root weight per plot in kg) --> **FYLD** (fresh root yield in tons per hectare)

$$FYLD = \frac{RTWT_{kg / plot}}{MaxHarvestedPlantsPerPlot \times PlantSpacing}\times10$$
**NOTE:** *MaxHarvestedPlantsPerPlot* in formula above is to distinguish from the *plantsPerPlot* meta-data field, in case that a net-plot harvest is used. In other words, the value should be the total number of plants intended for harvest in a plot, assuming there were no missing plants in the plot.

*PlantSpacing* is the area in $m^2$ per plant.


```r
dbdata %>% 
     count(studyYear,studyName,studyDesign,plotWidth,plotLength,plantsPerPlot) %>% 
     mutate(plotArea=plotWidth*plotLength) %>% 
     gt() %>% tab_options(table.font.size = pct(67))
```

```{=html}
<div id="kgrrrvxhgd" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#kgrrrvxhgd .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 67%;
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

#kgrrrvxhgd .gt_heading {
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

#kgrrrvxhgd .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#kgrrrvxhgd .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#kgrrrvxhgd .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kgrrrvxhgd .gt_col_headings {
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

#kgrrrvxhgd .gt_col_heading {
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

#kgrrrvxhgd .gt_column_spanner_outer {
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

#kgrrrvxhgd .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#kgrrrvxhgd .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#kgrrrvxhgd .gt_column_spanner {
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

#kgrrrvxhgd .gt_group_heading {
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

#kgrrrvxhgd .gt_empty_group_heading {
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

#kgrrrvxhgd .gt_from_md > :first-child {
  margin-top: 0;
}

#kgrrrvxhgd .gt_from_md > :last-child {
  margin-bottom: 0;
}

#kgrrrvxhgd .gt_row {
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

#kgrrrvxhgd .gt_stub {
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

#kgrrrvxhgd .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#kgrrrvxhgd .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#kgrrrvxhgd .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#kgrrrvxhgd .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#kgrrrvxhgd .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#kgrrrvxhgd .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kgrrrvxhgd .gt_footnotes {
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

#kgrrrvxhgd .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#kgrrrvxhgd .gt_sourcenotes {
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

#kgrrrvxhgd .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#kgrrrvxhgd .gt_left {
  text-align: left;
}

#kgrrrvxhgd .gt_center {
  text-align: center;
}

#kgrrrvxhgd .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#kgrrrvxhgd .gt_font_normal {
  font-weight: normal;
}

#kgrrrvxhgd .gt_font_bold {
  font-weight: bold;
}

#kgrrrvxhgd .gt_font_italic {
  font-style: italic;
}

#kgrrrvxhgd .gt_super {
  font-size: 65%;
}

#kgrrrvxhgd .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">studyYear</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyDesign</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">plotWidth</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">plotLength</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">plantsPerPlot</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">n</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">plotArea</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">4</td>
<td class="gt_row gt_right">4.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">125</td>
<td class="gt_row gt_right">16.0</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">4</td>
<td class="gt_row gt_right">4.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">68</td>
<td class="gt_row gt_right">16.0</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">4</td>
<td class="gt_row gt_right">4.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">72</td>
<td class="gt_row gt_right">16.0</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_left">RCBD</td>
<td class="gt_row gt_right">4</td>
<td class="gt_row gt_right">4.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">66</td>
<td class="gt_row gt_right">16.0</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">4</td>
<td class="gt_row gt_right">4.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">72</td>
<td class="gt_row gt_right">16.0</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">2</td>
<td class="gt_row gt_right">4.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">270</td>
<td class="gt_row gt_right">8.0</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_left">Alpha</td>
<td class="gt_row gt_right">3</td>
<td class="gt_row gt_right">2.5</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">273</td>
<td class="gt_row gt_right">7.5</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_left">RCBD</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">2.5</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">777</td>
<td class="gt_row gt_right">2.5</td></tr>
    <tr><td class="gt_row gt_right">2019</td>
<td class="gt_row gt_left">19geneticgainUB</td>
<td class="gt_row gt_left">Augmented</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">2.5</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">810</td>
<td class="gt_row gt_right">2.5</td></tr>
  </tbody>
  
  
</table>
</div>
```
In the example trial data, the **`plantsPerPlot`** meta-data field is empty. To my knowledge, no meta-data field is available in BreedBase to represent a net-plot harvest. 

***RECOMMEND INPUTING plantsPerPlot meta-data to cassavabase for _your_ breeding program!***

Luckily, since there are only 9 trials and this is a tutorial, we will decisions manually.

Firstly noting that the trial **19geneticgainUB** actually _does not_ have phenotypes (for any trait). It will be excluded downstream. (*~~I might find a substitute genetic gain trial, from an earlier year, for the sake of this example~~*)

To decide what the real **MaxHarvestedPlantsPerPlot** and **plantsPerPlot** were likely to have been, I make two plots below and also compute the maximum NOHAV for each trial. 


```r
dbdata %>% 
     ggplot(.,aes(x=NOHAV, fill=studyName)) + geom_density(alpha=0.75)
#> Warning: Removed 921 rows containing non-finite values
#> (stat_density).
```

<img src="03-clean_phenos_files/figure-html/unnamed-chunk-15-1.png" width="672" />
Maybe clearer to make a boxplot?

```r
dbdata %>% 
     # plot area in meters squared
     mutate(plotArea=plotWidth*plotLength) %>% 
     ggplot(.,aes(x=plotArea,y=NOHAV, fill=studyName)) + 
     geom_boxplot() + theme(axis.text.x = element_blank())
#> Warning: Removed 921 rows containing non-finite values
#> (stat_boxplot).
```

<img src="03-clean_phenos_files/figure-html/unnamed-chunk-16-1.png" width="672" />


```r
plantsPerPlot_choices<-dbdata %>% 
     distinct(studyYear,studyName,plotWidth,plotLength,plantsPerPlot) %>% 
     left_join(dbdata %>% 
                    group_by(studyName) %>% 
                    summarize(MaxNOHAV=max(NOHAV, na.rm=T))) %>% 
          # plot area in meters squared
     mutate(plotArea=plotWidth*plotLength,
            # Number of plants per plot
            plantsPerPlot=MaxNOHAV,
            plantsPerPlot=ifelse(studyName=="19.GS.C2.UYT.36.setA.UB",20,plantsPerPlot)) %>% 
     # exclude the empty genetic gain trial
     filter(studyName!="19geneticgainUB") %>% 
     select(studyName,plotArea,MaxNOHAV,plantsPerPlot)
#> Warning in max(NOHAV, na.rm = T): no non-missing arguments
#> to max; returning -Inf
#> Joining, by = "studyName"
plantsPerPlot_choices %>% gt() #%>% tab_options(table.font.size = pct(67))
```

```{=html}
<div id="sqozbjxjbk" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#sqozbjxjbk .gt_table {
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

#sqozbjxjbk .gt_heading {
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

#sqozbjxjbk .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#sqozbjxjbk .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#sqozbjxjbk .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#sqozbjxjbk .gt_col_headings {
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

#sqozbjxjbk .gt_col_heading {
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

#sqozbjxjbk .gt_column_spanner_outer {
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

#sqozbjxjbk .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#sqozbjxjbk .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#sqozbjxjbk .gt_column_spanner {
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

#sqozbjxjbk .gt_group_heading {
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

#sqozbjxjbk .gt_empty_group_heading {
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

#sqozbjxjbk .gt_from_md > :first-child {
  margin-top: 0;
}

#sqozbjxjbk .gt_from_md > :last-child {
  margin-bottom: 0;
}

#sqozbjxjbk .gt_row {
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

#sqozbjxjbk .gt_stub {
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

#sqozbjxjbk .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#sqozbjxjbk .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#sqozbjxjbk .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#sqozbjxjbk .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#sqozbjxjbk .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#sqozbjxjbk .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#sqozbjxjbk .gt_footnotes {
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

#sqozbjxjbk .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#sqozbjxjbk .gt_sourcenotes {
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

#sqozbjxjbk .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#sqozbjxjbk .gt_left {
  text-align: left;
}

#sqozbjxjbk .gt_center {
  text-align: center;
}

#sqozbjxjbk .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#sqozbjxjbk .gt_font_normal {
  font-weight: normal;
}

#sqozbjxjbk .gt_font_bold {
  font-weight: bold;
}

#sqozbjxjbk .gt_font_italic {
  font-style: italic;
}

#sqozbjxjbk .gt_super {
  font-size: 65%;
}

#sqozbjxjbk .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">studyName</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">plotArea</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">MaxNOHAV</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">plantsPerPlot</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">19.GS.C1.C2.C3.AYT.42.UB</td>
<td class="gt_row gt_right">16.0</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">10</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.IB</td>
<td class="gt_row gt_right">16.0</td>
<td class="gt_row gt_right">20</td>
<td class="gt_row gt_right">20</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setA.UB</td>
<td class="gt_row gt_right">16.0</td>
<td class="gt_row gt_right">18</td>
<td class="gt_row gt_right">20</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.IB</td>
<td class="gt_row gt_right">16.0</td>
<td class="gt_row gt_right">20</td>
<td class="gt_row gt_right">20</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C2.UYT.36.setB.UB</td>
<td class="gt_row gt_right">16.0</td>
<td class="gt_row gt_right">20</td>
<td class="gt_row gt_right">20</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.135.UB</td>
<td class="gt_row gt_right">8.0</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">10</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4B.PYT.140.IB</td>
<td class="gt_row gt_right">7.5</td>
<td class="gt_row gt_right">9</td>
<td class="gt_row gt_right">9</td></tr>
    <tr><td class="gt_row gt_left">19.GS.C4C.CE.822.IB</td>
<td class="gt_row gt_right">2.5</td>
<td class="gt_row gt_right">7</td>
<td class="gt_row gt_right">7</td></tr>
  </tbody>
  
  
</table>
</div>
```

For the sake of this example, it is 'ok' to make choices on the basis I have just done. 

***As a data generator, in-house at each breeding program, no reason not to get the correct answer and repair the metadata on the database!***


```r
dbdata %<>%
     # remove the empty genetic gain trial
     filter(studyName!="19geneticgainUB") %>% 
     select(-plantsPerPlot) %>% 
     # join plantsPerPlot_choices to the trial data
     left_join(plantsPerPlot_choices) %>% 
     # compute fresh root yield (FYLD) in tons per hectare
     mutate(PlantSpacing=plotArea/plantsPerPlot,
            FYLD=RTWT/(plantsPerPlot*PlantSpacing)*10)
#> Joining, by = "studyName"
dbdata %>% ggplot(.,aes(x=FYLD,fill=studyName)) + geom_density(alpha=0.75)
#> Warning: Removed 133 rows containing non-finite values
#> (stat_density).
```

<img src="03-clean_phenos_files/figure-html/unnamed-chunk-18-1.png" width="672" />

Additional things to compute:

1.  log-transform yield traits: this is a habit based on experience. Linear mixed-models should have normally distributed homoskedastic residuals, if they don't log-transform the response variable often helps. For FYLD and related traits, I always log-transform.


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
#> Warning: Removed 133 rows containing non-finite values
#> (stat_density).
```

<img src="03-clean_phenos_files/figure-html/unnamed-chunk-19-1.png" width="672" />

Debatable whether this is better. Let's not dwell on it. Onward!

**SUGGESTION:** For individuals working this manual, consider making different, or no transformations as you see fit, with your own data. Even better, set-up a direct comparison of results with- vs. without-transformation.*

## Save "cleaned" phenotypes {#save_cleaned_phenos}


```r
saveRDS(dbdata,file=here::here("output","phenotypes_cleaned.rds"))
```