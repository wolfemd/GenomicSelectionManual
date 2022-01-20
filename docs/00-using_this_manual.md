# Using this manual

## Workflow notes

-   Create a `workflowR` repository for your genomic prediction analysis, following [instructions here](#create_project).

-   Follow along the following documents as templates and examples:

    -   **Genomic Selection Manual (THIS BOOK)**
    -   **GS Process Map (\~TO BE LINKED HERE, AND SHOWN THROUGHOUT THE MANUAL\~)**
    -   **GS Checklist (\~TO BE LINKED HERE\~)**

-   Use some variant of these documents and code examples to complete a genomic prediction analysis and develop a report on the results.

-   **Advice and best practices:**

    -   Choose your own data, traits, from cassavabase

    -   Work through what the example code actually does for yourself.

        -   Follow-up functions you don't know by going to the manual.
        -   I will strive to provide references to other tutorials, papers, etc. to give context and help you learn more detail where desired..

    -   Inevitably, we will want divergences, alterations, bells-and-whistles on top of the process documented. **SUGGEST** altering and developing your own **process maps** and **checklists** as you go.

    -   Use a combination of Rmarkdown (**.Rmd**) and Rscripts (**.R**) to document your analysis, as demonstrated.

    -   Take the time to write commentary throughout. In full sentences, what do you intend to do? How do you interpret the results? What is the next step? Etc.

    -   Take the time to think through the naming of datasets, files, folders, R objects, etc.

    -   Use Git version control, made easy with Rstudio.

    -   Publish your code to GitHub and a report on your results as a webpage using GitHub Pages. I will demonstrate using the package `workflowR` to manage these aspects.

## R sessions, packages to load

I will use the `tidyverse` and also `genomicMateSelectR` packages throughout the pipeline.

There are others that may appear.

I recommend, for each pipeline segment, starting with a new R session. Begin each segment, with a step to load these R packages:


```r
library(tidyverse)
library(genomicMateSelectR)
library(gt) # just for the nice looking tables
```

## Software to install

-   R

-   Rstudio

-   Packages

    -   tidyverse
    -   `genomicMateSelectR`

-   Bioinformatics tools

    -   `vcftools`
    -   `bcftools`

## High performance and remote computing