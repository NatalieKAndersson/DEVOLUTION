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

Try loading dataset "test". You should not get any error message doing this.

We can now load the package to the global environment.

```require(DEVOLUTION)```

## Usage

In order to illustrate the usage of the package, let's go through an example using the fabricated data set built into the software named "Segmentfile".

**Load the data**

What does the data set look like?

Describe the columns.

**Choose parameters**

You have to determine the following parameters.


**The event matrix**

Let us now make an event matrix using the fabricated data set.

DEVOLUTION()

Show an example of the event matrix and what it means.

**Let's produce the final event matrix**

What has the program done?

Table of sizes.

**Phylogenetic trees**

In the end one can use the event matrix in order to reconstruct a phylogenetic tree.

