# Nails
This Repository provides the code for analysis of data to support the submission of the manuscript titled "Cell-free DNA from Nail Clippings as a Source of matched normal control for comprehensive Genomic Studies in hematologic malignancies. "

This repository is self-contained and includes all data to generate the figures being represented.  The data is deidentified data from the supplementary data provided with the manuscript and additional manifests required to generate figures and to complete statistical analysis. 

To acquire the repository, clone to machine with UNIX/Linux based OS. Local installation of git is required.

```bash
git clone https://github.com/chadvanderbilt/Nails.git
cd ./Nails
```

To Run the code please refer to install.R for proper software versioning and this can be ran as an Rscript to ensure installation of requisite packages.  All packages necessary are available on CRAN for R major version 4. 


```bash
Rscript install.R
cd scripts
```

Using editor change line 47 in Complete_Figure_with_Stats.R to set working director to the data folder from repository. Once saved, execute with Rscript. 

```bash
Rscript  Complete_Figure_with_Stats.R
```

Figures will be generated in the plot directory.  Please reach out to me (vanderbc@mskcc.org) with questions. 








