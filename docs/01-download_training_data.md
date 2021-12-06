# Download training data

## Cassavabase download

Go to [Cassavabase](https://www.cassavabase.org/) or your favorite alternative BreedBase.

Login.

Go to the **Search \> Wizard**

![](images/search_wizard.png){width="300"}

## Example dataset

For the sake of example, I will choose a small, but real dataset that I think will exemplify how data should be stored on the DB. Further, I choose data that have key features, including pedigree relationship and a high-proportion of genotyped accessions.

-   Use Wizard (Breeding Program: IITA > Years: 2013 > Locations: Ibadan >

-   Choose 2 trials `13geneticgainIB` and `13.GS.C1.CE.864.IB` --> create list **"ExampleIITAtrials"**

-   Clear Wizard panes

-   Start from new list **"ExampleIITAtrials" >** Accessions: "Select All" --> create list **"ExampleIITAaccesions"**

    ![](images/createExampleAccessionsList.png)

-   Genotyping Protocols > **"West Africa Clones Dart-GBS 2020"**

-   **Download Genotypes, Phenotypes and Metadata (**Bottom right on Search Wizard page)

    -   **"Related Genotype Data"**

        -   Download genotypes ("VCF File Format")

        <!-- -->

        -   This is going to take a while..... it timed out when my compute went to sleep, but when I came back and refreshed the site, the file was ready-to-go and downloaded right away.

    -   **"Related Trial Metadata"** --> exports a .csv named "phenotypes.csv", which is awkward....

    -   **"Related Trial Phenotypes"** --> exports a .csv also named "phenotypes.csv" so auto-changed to "phenotype (1).csv"

        -   **Manually rename and store here:** `data/phenotype.csv` and `data/metadata.csv`

-   **Download Pedigree** (on hold / blocked)

    -   Download > Download Pedigree

    -   "ExampleIITAaccessions" list > Direct parents only

        -   ... list failed to validate....

        -   Click "lists" in top right > Click the "Validate" checkmark on the row for the list of accessions..... "Working....." --> says list passed validation...

        -   Working with DB team...
