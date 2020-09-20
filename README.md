# DEVOLUTION
Phylogenetic reconstruction from multiregional sampling data

## Installation instructions

```{r, eval=FALSE, echo=TRUE}
install.packages(c("readxl","xlsxjars","rJava","xlsx","phangorn","ape",
                   "tree","writexl","tidyverse","colorspace","ggplot2","ggtree",
                   "phytools","bbmle","picante","ggforce","stringr","plotrix",
                   "reshape2","RColorBrewer","cluster","clValid","factoextra","NbClust","dplyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEVOLUTION")
```

If no error message has appeared, continue with the installation of the DEVOLUTION package. Install it in the folder in which your data files are (by replacint PATH with this) or make sure to set the path to that location after installation. You have to be in the folder of your files in order for the software to find them.

```install.packages("~PATH/DEVOLUTION_0.1.0.tar.gz", repos = NULL, type = "source")```

We can now load the package to the global environment.

```require(DEVOLUTION)```

Try loading dataset "test". You should not get any error message doing this.

## Usage

In order to illustrate the usage of the package, let's go through an example using the fabricated data set built into the package named "Segment". This is a made up data set representing a SNP-array output segment file.

**Load the data**

Let's load the data.

```test <- load_matrix(filename="Segment.xlsx",sheetname ="Tumor1")```

And take a look at the data itself

```head(test)```

What does the data set look like?

Describe the columns.

It contains data from two tumors. Let's extract the data from one of them.

```
x <- "Tumor1" #The tumor we want to look at.
datasegment <- splitdata(test,x)```

**Choose parameters**

The following parameters can be altered. If you do not want to change them, go with the default shown below.

```
event_co <- 1000000 #Cutoff for events.
sub_co <- 0.9 #Cutoff for stem events.
eventnumber <- 300 #The maximum number of events I think a particular sample or subclone will have.
TDS <- "No" #"Yes", "No" or "Only". Choose wether or not you want TDS-events to be included in the following computations. If you choose yes we will not remove anything.
```

**The event matrix**

Let us now make an event matrix using the fabricated data set.

```EM_test <- DEVOLUTION(datasegment,event_co,sub_co,eventnumber,TDS)```

Show an example of the event matrix and what it means.

**Let's produce the final event matrix**

What has the program done?

Table of sizes.

**Phylogenetic trees**

In the end one can use the event matrix in order to reconstruct a phylogenetic tree.

