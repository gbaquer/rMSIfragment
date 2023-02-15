# rMSIfragment
rMSIfragment is an open-source R package for the automated annotation of in-source fragments for improved lipid annotation in MALDI-MSI.

More information can be found in the accompanying preprint: 
Gerard Baquer, Lluc Sementé, Pere Ràfols, Lucía Martín-Saiz, Christoph Bookmeyer, José A. Fernández,  Xavier Correig & María García-Altares. rMSIfragment: Automated Annotation of In-Source Fragments for Improved Lipid Annotation in MALDI-MSI. BioRxiv (2023). https://doi.org/xxxx

### 1. Installation

Install devtools:
```R
>  install.packages("devtools")
```
Install rMSI and rMSIproc
```R
>  devtools::install_github("prafols/rMSI", ref = "0.8")
>  devtools::install_github("prafols/rMSIproc", ref = "0.2")
```
Install rMSIfragment from github:
```R
> devtools::install_github("gbaquer/rMSIfragment")
```
This will install rMSIfragment package and all of its dependencies in your R environment. Then you can access its functions by loading the rMSIfragment package or through the `::` operator.

### 2. Basic Usage
```R
## 2.1. Load Data
pks<-rMSIproc::LoadPeakMatrix("[Full Path to Peak Matrix (.zip)]")

## 2.2. Annotate in-source fragments
[Update]
```
### 3. Downloading sample data
To easily try out the functionality of the package we provide a sample datasets available at https://doi.org/10.17632/53grw3ys6y.1

Baquer, Gerard (2022), “rMSIfragment datasets G1-G15”, Mendeley Data, V1, doi: 10.17632/53grw3ys6y.1

### 4. Processing your own imzML
rMSIfragment uses data in the rMSIproc format. To annotate your own data you will have to process the imzML using the following command:

```R
> rMSIproc::ProcessWizard()
```
A window will appear where you can set up all the processing parameters, the input data and the output directory to store the results.

Refer to the repository of rMSIproc for further details (https://github.com/prafols/rMSIproc)
