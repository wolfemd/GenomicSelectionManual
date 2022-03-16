# Prepare phenotype data



-   **Context and Purpose:** In this step, we do quality control, clean and format training data for further analysis.
-   **Upstream:** Section \@ref(download-training-data) - training data download
-   **Downstream:** pretty much everything
-   **Inputs:** "Raw" field trial data
-   **Expected outputs:** "Cleaned" field trial data

## Process Map

![](images/prepare_phenos_process_map.png){width=100%}

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
#> 1             plot 1613
# table(dbdata$observationLevel)
```

Only plot-basis in this case.

**Checklist:** What experimental designs are present? How are they represented by the variables in the dataset? Are all designs consistent with your expectations, for example relative to the reported "`trialType`", "`studyName`" and/or "`studyDesign`"?

In this step, in the past, I have not been certain of the experimental designs of the trials I had downloaded. I was also not certain how the designs were represented in the column-names. For this reason, I developed an *ad hoc* custom code to "detect" the designs. I built the `genomicMateSelectR` function `detectExptDesigns()`. See an example [here](https://wolfemd.github.io/IITA_2021GS/01-cleanTPdata.html#Detect_experimental_designs).

**RECOMMENDATION:** Each analyst needs to use exploratory data anlaysis, making summary statistics and plots as necessary to determine how the data should be modelled downstream. If there are missing or incorrectly represented trial design variables, ***get it corrected on the database (contact breeding program data manager, as necessary)***.

Because I have a small example dataset, it is possible to look at 9 trials and evaluate.

Often, *many* more trials are part of a genomic prediction. This is why it is essential that trial designs be consistent, and clear to the analyst. You may need to derive a strategy similar to the `detectExptDesigns()` function to semi-automate the process.


```r
library(gt)
dbdata %>% 
     count(studyName,trialType, studyDesign, numberBlocks,numberReps,entryType) %>% 
     spread(entryType,n) %>% 
     gt()  %>% 
     tab_options(table.font.size = pct(75))
```

```{=html}
<div id="qjxdaxcqwk" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#qjxdaxcqwk .gt_table {
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

#qjxdaxcqwk .gt_heading {
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

#qjxdaxcqwk .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#qjxdaxcqwk .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#qjxdaxcqwk .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#qjxdaxcqwk .gt_col_headings {
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

#qjxdaxcqwk .gt_col_heading {
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

#qjxdaxcqwk .gt_column_spanner_outer {
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

#qjxdaxcqwk .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#qjxdaxcqwk .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#qjxdaxcqwk .gt_column_spanner {
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

#qjxdaxcqwk .gt_group_heading {
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

#qjxdaxcqwk .gt_empty_group_heading {
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

#qjxdaxcqwk .gt_from_md > :first-child {
  margin-top: 0;
}

#qjxdaxcqwk .gt_from_md > :last-child {
  margin-bottom: 0;
}

#qjxdaxcqwk .gt_row {
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

#qjxdaxcqwk .gt_stub {
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

#qjxdaxcqwk .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#qjxdaxcqwk .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#qjxdaxcqwk .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#qjxdaxcqwk .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#qjxdaxcqwk .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#qjxdaxcqwk .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#qjxdaxcqwk .gt_footnotes {
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

#qjxdaxcqwk .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#qjxdaxcqwk .gt_sourcenotes {
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

#qjxdaxcqwk .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#qjxdaxcqwk .gt_left {
  text-align: left;
}

#qjxdaxcqwk .gt_center {
  text-align: center;
}

#qjxdaxcqwk .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#qjxdaxcqwk .gt_font_normal {
  font-weight: normal;
}

#qjxdaxcqwk .gt_font_bold {
  font-weight: bold;
}

#qjxdaxcqwk .gt_font_italic {
  font-style: italic;
}

#qjxdaxcqwk .gt_super {
  font-size: 65%;
}

#qjxdaxcqwk .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_left">Clonal Evaluation</td>
<td class="gt_row gt_left">Augmented</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">25</td>
<td class="gt_row gt_right">839</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_left">Clonal Evaluation</td>
<td class="gt_row gt_left">Augmented</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">51</td>
<td class="gt_row gt_right">698</td></tr>
  </tbody>
  
  
</table>
</div>
```

Summary table above shows:

1)  **trialType** and **studyDesign** cannot be 100% relied upon, at least not here.
2)  The only trial actually listed as having `studyDesign=="Augmented"` does not have "check" vs. "test" distinguished in the "entryType".
3)  A `trialType=="Clonal Evaluation"` with `studyDesign=="RCBD"` but actually only 1 replication.

Next, I'll check if the **replicate** and **blockNumber** columns reliably distinguish complete and incomplete blocks in the data.


```r
dbdata %>% 
     group_by(studyName) %>% 
     summarize(N_replicate=length(unique(replicate)),
               N_blockNumber=length(unique(blockNumber))) %>% 
     gt() %>% tab_options(table.font.size = pct(75))
```

```{=html}
<div id="xsuzcdvpet" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#xsuzcdvpet .gt_table {
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

#xsuzcdvpet .gt_heading {
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

#xsuzcdvpet .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#xsuzcdvpet .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#xsuzcdvpet .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#xsuzcdvpet .gt_col_headings {
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

#xsuzcdvpet .gt_col_heading {
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

#xsuzcdvpet .gt_column_spanner_outer {
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

#xsuzcdvpet .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#xsuzcdvpet .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#xsuzcdvpet .gt_column_spanner {
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

#xsuzcdvpet .gt_group_heading {
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

#xsuzcdvpet .gt_empty_group_heading {
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

#xsuzcdvpet .gt_from_md > :first-child {
  margin-top: 0;
}

#xsuzcdvpet .gt_from_md > :last-child {
  margin-bottom: 0;
}

#xsuzcdvpet .gt_row {
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

#xsuzcdvpet .gt_stub {
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

#xsuzcdvpet .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#xsuzcdvpet .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#xsuzcdvpet .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#xsuzcdvpet .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#xsuzcdvpet .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#xsuzcdvpet .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#xsuzcdvpet .gt_footnotes {
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

#xsuzcdvpet .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#xsuzcdvpet .gt_sourcenotes {
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

#xsuzcdvpet .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#xsuzcdvpet .gt_left {
  text-align: left;
}

#xsuzcdvpet .gt_center {
  text-align: center;
}

#xsuzcdvpet .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#xsuzcdvpet .gt_font_normal {
  font-weight: normal;
}

#xsuzcdvpet .gt_font_bold {
  font-weight: bold;
}

#xsuzcdvpet .gt_font_italic {
  font-style: italic;
}

#xsuzcdvpet .gt_super {
  font-size: 65%;
}

#xsuzcdvpet .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td></tr>
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
<div id="zzrrvhnzju" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#zzrrvhnzju .gt_table {
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

#zzrrvhnzju .gt_heading {
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

#zzrrvhnzju .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#zzrrvhnzju .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#zzrrvhnzju .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zzrrvhnzju .gt_col_headings {
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

#zzrrvhnzju .gt_col_heading {
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

#zzrrvhnzju .gt_column_spanner_outer {
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

#zzrrvhnzju .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#zzrrvhnzju .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#zzrrvhnzju .gt_column_spanner {
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

#zzrrvhnzju .gt_group_heading {
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

#zzrrvhnzju .gt_empty_group_heading {
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

#zzrrvhnzju .gt_from_md > :first-child {
  margin-top: 0;
}

#zzrrvhnzju .gt_from_md > :last-child {
  margin-bottom: 0;
}

#zzrrvhnzju .gt_row {
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

#zzrrvhnzju .gt_stub {
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

#zzrrvhnzju .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#zzrrvhnzju .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#zzrrvhnzju .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#zzrrvhnzju .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#zzrrvhnzju .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#zzrrvhnzju .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zzrrvhnzju .gt_footnotes {
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

#zzrrvhnzju .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#zzrrvhnzju .gt_sourcenotes {
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

#zzrrvhnzju .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#zzrrvhnzju .gt_left {
  text-align: left;
}

#zzrrvhnzju .gt_center {
  text-align: center;
}

#zzrrvhnzju .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#zzrrvhnzju .gt_font_normal {
  font-weight: normal;
}

#zzrrvhnzju .gt_font_bold {
  font-weight: bold;
}

#zzrrvhnzju .gt_font_italic {
  font-style: italic;
}

#zzrrvhnzju .gt_super {
  font-size: 65%;
}

#zzrrvhnzju .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_center">TRUE</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_center">TRUE</td></tr>
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
<div id="eiyreqboxm" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#eiyreqboxm .gt_table {
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

#eiyreqboxm .gt_heading {
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

#eiyreqboxm .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#eiyreqboxm .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#eiyreqboxm .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#eiyreqboxm .gt_col_headings {
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

#eiyreqboxm .gt_col_heading {
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

#eiyreqboxm .gt_column_spanner_outer {
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

#eiyreqboxm .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#eiyreqboxm .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#eiyreqboxm .gt_column_spanner {
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

#eiyreqboxm .gt_group_heading {
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

#eiyreqboxm .gt_empty_group_heading {
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

#eiyreqboxm .gt_from_md > :first-child {
  margin-top: 0;
}

#eiyreqboxm .gt_from_md > :last-child {
  margin-bottom: 0;
}

#eiyreqboxm .gt_row {
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

#eiyreqboxm .gt_stub {
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

#eiyreqboxm .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#eiyreqboxm .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#eiyreqboxm .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#eiyreqboxm .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#eiyreqboxm .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#eiyreqboxm .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#eiyreqboxm .gt_footnotes {
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

#eiyreqboxm .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#eiyreqboxm .gt_sourcenotes {
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

#eiyreqboxm .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#eiyreqboxm .gt_left {
  text-align: left;
}

#eiyreqboxm .gt_center {
  text-align: center;
}

#eiyreqboxm .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#eiyreqboxm .gt_font_normal {
  font-weight: normal;
}

#eiyreqboxm .gt_font_bold {
  font-weight: bold;
}

#eiyreqboxm .gt_font_italic {
  font-style: italic;
}

#eiyreqboxm .gt_super {
  font-size: 65%;
}

#eiyreqboxm .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_right">828</td>
<td class="gt_row gt_right">828</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_right">692</td>
<td class="gt_row gt_right">692</td></tr>
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
<div id="pdlipzcqkr" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#pdlipzcqkr .gt_table {
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

#pdlipzcqkr .gt_heading {
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

#pdlipzcqkr .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#pdlipzcqkr .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#pdlipzcqkr .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#pdlipzcqkr .gt_col_headings {
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

#pdlipzcqkr .gt_col_heading {
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

#pdlipzcqkr .gt_column_spanner_outer {
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

#pdlipzcqkr .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#pdlipzcqkr .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#pdlipzcqkr .gt_column_spanner {
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

#pdlipzcqkr .gt_group_heading {
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

#pdlipzcqkr .gt_empty_group_heading {
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

#pdlipzcqkr .gt_from_md > :first-child {
  margin-top: 0;
}

#pdlipzcqkr .gt_from_md > :last-child {
  margin-bottom: 0;
}

#pdlipzcqkr .gt_row {
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

#pdlipzcqkr .gt_stub {
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

#pdlipzcqkr .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#pdlipzcqkr .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#pdlipzcqkr .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#pdlipzcqkr .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#pdlipzcqkr .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#pdlipzcqkr .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#pdlipzcqkr .gt_footnotes {
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

#pdlipzcqkr .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#pdlipzcqkr .gt_sourcenotes {
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

#pdlipzcqkr .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#pdlipzcqkr .gt_left {
  text-align: left;
}

#pdlipzcqkr .gt_center {
  text-align: center;
}

#pdlipzcqkr .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#pdlipzcqkr .gt_font_normal {
  font-weight: normal;
}

#pdlipzcqkr .gt_font_bold {
  font-weight: bold;
}

#pdlipzcqkr .gt_font_italic {
  font-style: italic;
}

#pdlipzcqkr .gt_super {
  font-size: 65%;
}

#pdlipzcqkr .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_right">828</td>
<td class="gt_row gt_right">828</td>
<td class="gt_row gt_right">828</td>
<td class="gt_row gt_right">828</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_right">692</td>
<td class="gt_row gt_right">692</td>
<td class="gt_row gt_right">692</td>
<td class="gt_row gt_right">692</td></tr>
  </tbody>
  
  
</table>
</div>
```

From this, you can see that except for `studyName=="19.GS.C1.C2.C3.AYT.42.UB"` the sub-blocks represented by **blockNumber** have only subsets of the total number of accessions in the trial, as expected.

Further, except for `studyName=="19geneticgainUB"` all trials have pretty consistently sized sub-blocks.

Now I will *ad hoc* create two variables (**CompleteBlocks** and **IncompleteBlocks**), indicating (`TRUE/FALSE`) whether to model using the **replicate** and/or **blockNumber** variable.

I *also* like to create *explicitly nested* design variables (**yearInLoc**, **trialInLocYr**, **repInTrial**, **blockInRep**).


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
<div id="ndwcvqmrvr" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#ndwcvqmrvr .gt_table {
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

#ndwcvqmrvr .gt_heading {
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

#ndwcvqmrvr .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ndwcvqmrvr .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ndwcvqmrvr .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ndwcvqmrvr .gt_col_headings {
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

#ndwcvqmrvr .gt_col_heading {
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

#ndwcvqmrvr .gt_column_spanner_outer {
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

#ndwcvqmrvr .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ndwcvqmrvr .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ndwcvqmrvr .gt_column_spanner {
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

#ndwcvqmrvr .gt_group_heading {
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

#ndwcvqmrvr .gt_empty_group_heading {
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

#ndwcvqmrvr .gt_from_md > :first-child {
  margin-top: 0;
}

#ndwcvqmrvr .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ndwcvqmrvr .gt_row {
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

#ndwcvqmrvr .gt_stub {
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

#ndwcvqmrvr .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ndwcvqmrvr .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#ndwcvqmrvr .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ndwcvqmrvr .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ndwcvqmrvr .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ndwcvqmrvr .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ndwcvqmrvr .gt_footnotes {
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

#ndwcvqmrvr .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#ndwcvqmrvr .gt_sourcenotes {
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

#ndwcvqmrvr .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#ndwcvqmrvr .gt_left {
  text-align: left;
}

#ndwcvqmrvr .gt_center {
  text-align: center;
}

#ndwcvqmrvr .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ndwcvqmrvr .gt_font_normal {
  font-weight: normal;
}

#ndwcvqmrvr .gt_font_bold {
  font-weight: bold;
}

#ndwcvqmrvr .gt_font_italic {
  font-style: italic;
}

#ndwcvqmrvr .gt_super {
  font-size: 65%;
}

#ndwcvqmrvr .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_right">864</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_center">FALSE</td>
<td class="gt_row gt_right">749</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">1</td></tr>
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
<div id="tcnxpktbxy" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#tcnxpktbxy .gt_table {
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

#tcnxpktbxy .gt_heading {
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

#tcnxpktbxy .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#tcnxpktbxy .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#tcnxpktbxy .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tcnxpktbxy .gt_col_headings {
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

#tcnxpktbxy .gt_col_heading {
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

#tcnxpktbxy .gt_column_spanner_outer {
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

#tcnxpktbxy .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#tcnxpktbxy .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#tcnxpktbxy .gt_column_spanner {
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

#tcnxpktbxy .gt_group_heading {
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

#tcnxpktbxy .gt_empty_group_heading {
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

#tcnxpktbxy .gt_from_md > :first-child {
  margin-top: 0;
}

#tcnxpktbxy .gt_from_md > :last-child {
  margin-bottom: 0;
}

#tcnxpktbxy .gt_row {
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

#tcnxpktbxy .gt_stub {
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

#tcnxpktbxy .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tcnxpktbxy .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#tcnxpktbxy .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#tcnxpktbxy .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#tcnxpktbxy .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#tcnxpktbxy .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#tcnxpktbxy .gt_footnotes {
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

#tcnxpktbxy .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#tcnxpktbxy .gt_sourcenotes {
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

#tcnxpktbxy .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#tcnxpktbxy .gt_left {
  text-align: left;
}

#tcnxpktbxy .gt_center {
  text-align: center;
}

#tcnxpktbxy .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#tcnxpktbxy .gt_font_normal {
  font-weight: normal;
}

#tcnxpktbxy .gt_font_bold {
  font-weight: bold;
}

#tcnxpktbxy .gt_font_italic {
  font-style: italic;
}

#tcnxpktbxy .gt_super {
  font-size: 65%;
}

#tcnxpktbxy .gt_footnote_marks {
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
#> Warning in mask$eval_all_mutate(quo): NAs introduced by
#> coercion
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

$$FYLD = \frac{RTWT_{kg / plot}}{MaxHarvestedPlantsPerPlot \times PlantSpacing}\times10$$ **NOTE:** *MaxHarvestedPlantsPerPlot* in formula above is to distinguish from the *plantsPerPlot* meta-data field, in case that a net-plot harvest is used. In other words, the value should be the total number of plants intended for harvest in a plot, assuming there were no missing plants in the plot.

*PlantSpacing* is the area in $m^2$ per plant.


```r
dbdata %>% 
     count(studyYear,studyName,studyDesign,plotWidth,plotLength,plantsPerPlot) %>% 
     mutate(plotArea=plotWidth*plotLength) %>% 
     gt() %>% tab_options(table.font.size = pct(67))
```

```{=html}
<div id="inursmalds" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#inursmalds .gt_table {
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

#inursmalds .gt_heading {
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

#inursmalds .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#inursmalds .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#inursmalds .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#inursmalds .gt_col_headings {
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

#inursmalds .gt_col_heading {
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

#inursmalds .gt_column_spanner_outer {
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

#inursmalds .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#inursmalds .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#inursmalds .gt_column_spanner {
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

#inursmalds .gt_group_heading {
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

#inursmalds .gt_empty_group_heading {
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

#inursmalds .gt_from_md > :first-child {
  margin-top: 0;
}

#inursmalds .gt_from_md > :last-child {
  margin-bottom: 0;
}

#inursmalds .gt_row {
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

#inursmalds .gt_stub {
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

#inursmalds .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#inursmalds .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#inursmalds .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#inursmalds .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#inursmalds .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#inursmalds .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#inursmalds .gt_footnotes {
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

#inursmalds .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#inursmalds .gt_sourcenotes {
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

#inursmalds .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#inursmalds .gt_left {
  text-align: left;
}

#inursmalds .gt_center {
  text-align: center;
}

#inursmalds .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#inursmalds .gt_font_normal {
  font-weight: normal;
}

#inursmalds .gt_font_bold {
  font-weight: bold;
}

#inursmalds .gt_font_italic {
  font-style: italic;
}

#inursmalds .gt_super {
  font-size: 65%;
}

#inursmalds .gt_footnote_marks {
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
    <tr><td class="gt_row gt_right">2013</td>
<td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_left">Augmented</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">2.5</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">864</td>
<td class="gt_row gt_right">2.5</td></tr>
    <tr><td class="gt_row gt_right">2013</td>
<td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_left">Augmented</td>
<td class="gt_row gt_right">1</td>
<td class="gt_row gt_right">8.0</td>
<td class="gt_row gt_center">NA</td>
<td class="gt_row gt_right">749</td>
<td class="gt_row gt_right">8.0</td></tr>
  </tbody>
  
  
</table>
</div>
```

In the example trial data, the **`plantsPerPlot`** meta-data field is empty. To my knowledge, no meta-data field is available in BreedBase to represent a net-plot harvest.

***RECOMMEND INPUTING plantsPerPlot meta-data to cassavabase for your*** **breeding program!**

Luckily, since there are only 9 trials and this is a tutorial, we will decisions manually.

Firstly noting that the trial **19geneticgainUB** actually *does not* have phenotypes (for any trait). It will be excluded downstream. (~~*I might find a substitute genetic gain trial, from an earlier year, for the sake of this example*~~)

To decide what the real **MaxHarvestedPlantsPerPlot** and **plantsPerPlot** were likely to have been, I make two plots below and also compute the maximum NOHAV for each trial.


```r
dbdata %>% 
     ggplot(.,aes(x=NOHAV, fill=studyName)) + geom_density(alpha=0.75)
#> Warning: Removed 71 rows containing non-finite values
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
#> Warning: Removed 71 rows containing non-finite values
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
#> Joining, by = "studyName"
plantsPerPlot_choices %>% gt() #%>% tab_options(table.font.size = pct(67))
```

```{=html}
<div id="guwmwvdkad" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#guwmwvdkad .gt_table {
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

#guwmwvdkad .gt_heading {
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

#guwmwvdkad .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#guwmwvdkad .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#guwmwvdkad .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#guwmwvdkad .gt_col_headings {
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

#guwmwvdkad .gt_col_heading {
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

#guwmwvdkad .gt_column_spanner_outer {
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

#guwmwvdkad .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#guwmwvdkad .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#guwmwvdkad .gt_column_spanner {
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

#guwmwvdkad .gt_group_heading {
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

#guwmwvdkad .gt_empty_group_heading {
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

#guwmwvdkad .gt_from_md > :first-child {
  margin-top: 0;
}

#guwmwvdkad .gt_from_md > :last-child {
  margin-bottom: 0;
}

#guwmwvdkad .gt_row {
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

#guwmwvdkad .gt_stub {
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

#guwmwvdkad .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#guwmwvdkad .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#guwmwvdkad .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#guwmwvdkad .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#guwmwvdkad .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#guwmwvdkad .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#guwmwvdkad .gt_footnotes {
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

#guwmwvdkad .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#guwmwvdkad .gt_sourcenotes {
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

#guwmwvdkad .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#guwmwvdkad .gt_left {
  text-align: left;
}

#guwmwvdkad .gt_center {
  text-align: center;
}

#guwmwvdkad .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#guwmwvdkad .gt_font_normal {
  font-weight: normal;
}

#guwmwvdkad .gt_font_bold {
  font-weight: bold;
}

#guwmwvdkad .gt_font_italic {
  font-style: italic;
}

#guwmwvdkad .gt_super {
  font-size: 65%;
}

#guwmwvdkad .gt_footnote_marks {
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
    <tr><td class="gt_row gt_left">13.GS.C1.CE.864.IB</td>
<td class="gt_row gt_right">2.5</td>
<td class="gt_row gt_right">13</td>
<td class="gt_row gt_right">13</td></tr>
    <tr><td class="gt_row gt_left">13geneticgainIB</td>
<td class="gt_row gt_right">8.0</td>
<td class="gt_row gt_right">10</td>
<td class="gt_row gt_right">10</td></tr>
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
#> Warning: Removed 89 rows containing non-finite values
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
#> Warning: Removed 89 rows containing non-finite values
#> (stat_density).
```

<img src="03-clean_phenos_files/figure-html/unnamed-chunk-19-1.png" width="672" />

Debatable whether this is better. Let's not dwell on it. Onward!

**SUGGESTION:** For individuals working this manual, consider making different, or no transformations as you see fit, with your own data. Even better, set-up a direct comparison of results with- vs. without-transformation.\*

## Save "cleaned" phenotypes {#save_cleaned_phenos}


```r
saveRDS(dbdata,file=here::here("output","phenotypes_cleaned.rds"))
```
