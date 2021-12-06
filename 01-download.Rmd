# Download training data

Go to [Cassavabase](https://www.cassavabase.org/) or your favorite alternative BreedBase.

Login.

Go to the **Search \> Wizard**

![](images/search_wizard.png){width="300"}

## Example dataset

For the sake of example, I will choose a small, but real dataset that I think will exemplify how data should be stored on the DB. Further, I choose data that have key features, including pedigree relationship and a high-proportion of genotyped accessions.

### Create trial list

Create a list of trials using the "Wizard"

![](images/wizard_create_trial_list.png)

IITA trials at Ibadan and Ubiaja locations, planted in 2019. Further chose key trial types and specific trials as seen in screenshot.

Create list: **IITA_ExampleGStrials_2021Dec04**

### Download related trial data

Clear Wizard panes

Start from the new list of trials created: **"IITA_ExampleGStrials_2021Dec04"**

Download **"Related Trial Metadata"** and **"Related Trial Phenotypes"**

![](images/wizard_dl_related_trial_metadata.png)

![](images/wizard_dl_related_trial_phenotypes.png)

Exports .csv files `phenotype.csv` and `metadata.csv`.

Store in `data/` sub-directory for current project.

### Make an accession list

**\[NEW + EXPTL\]** Choose:\
**Genotyping Protocol:** "IITA DArT-GBS 08 Aug 2021", *then*\
**Accessions:** "Select All"\

Create list **"IITA_ExampleGSaccessions_2021Dec05"**.

![](images/wizard_make_accession_list_genotypingprotocolfirst.png)

### Download related trial genotype data

Download **"Related Trial Genotype Data",** choosing both available formats: VCF and Dosage Matrix (.tsv).\
![](images/wizard_dl_related_genotype_data.png){width="485"}\
**NOTE:** This probably will take a while and usually times out. However, Cassavabase will complete preparation of the file and have it ready-to-go when you return, when ready, it should begin downloading immediately.

### Download Pedigree

**(on hold / blocked)**

Go to "Manage \> [Download](https://www.cassavabase.org/breeders/download) \> Download Pedigree"

![](images/manage_download_dropdown.png){width="226"}

![](images/download_pedigree.png)

-   Same issue currently has been brought up from MusaBase: <https://github.com/solgenomics/sgn/issues/3865>
