# DEVOLUTION
Devoution is a tool for phylogenetic reconstruction from multiregional sampling data that can incorporate information from SNP-array, WES, WGS, TDS etc. using the mutated sample fraction as input. The mutated sample fraction is the proportion of cancer cells in a praticular biopsy that harbor an alteration.

## Setting up DEVOLUTION
Download the R script denoted "DEVOLUTION" and put it in a folder on your computer. Double click on the script to open it in your R-environment.

Before starting your analysis, we must load the dependencies of the algorithm. These are found in the beginning of the code.

```
library("readxl") #Needed to load the data from the xlsx file.
library("stringr") #Needed for using the function "word".
library("ape")
library("phangorn") #Needed to transform the EM into phyDat and make trees.
library("ggplot2") #Needed to visualize the trees.
library("ggtree")
library("ggimage") #Needed to insert the pies in the tree.
library("dplyr") #Needed for the distinct function in pie.it.
```
If they are not installet you can install them by using the following command.

```R
install.packages(c("readxl","stringr","ape","phangorn","ggplot2",
                   "ggtree","ggimage","dplyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEVOLUTION")
```
Also load all functions by marking them and pressing “Run”.

#Ladda in en bild här.

If no error message has appeared we are good to go!

## Usage
If no error message has appeared, we are now ready to load some data and get going! At the top of the script you have to set the path i.e. the location of the files that will be analyzed.
```R
setwd("~/yourpath")
```
In order to illustrate the usage of the algorithm, let’s go through an example using the example data set that also can be downloaded from the github page. The data file is named “Segment.xlsx”. This is a fabricated data set representing a SNP-array output segment file for two neuroblastomas.

```R
data <- load_matrix(filename="Segment_NB.xlsx",sheetname ="Tumor1")
head(data)
```
What does the data set look like? Describe the columns. Make sure that you understand what each of the 11 columns represents in this data file. The log2 ratio and VAF are not essential for the analysis. If not supplied, these columns will be replaced by "NA".

- Tumor ID: The ID for the biopsies from the same patient.
- Sample ID: The biopsy names.
- Chr, Start, End: The position of the genetic alteration.
- Med.LogR: The median log2R obtained through chromosome analysis. Not needed for the DEVOLUTION algorithm, so you could skip it or set it to NA.
- VAF: Variant Allele Frequency. Not used here.
- Method: The method by which this genetic aberration has been identified.
- Cytoband: The location of the genetic alteration or gene name.
- Clone size: The mutated sample fraction i.e. the fraction of cancer cells in this particular biopsy that has this aberration.

#Sätt in bild på header.

It contains data from two tumors. The splitdata function determines the start and end position of these data sets. Here data is the whole file you loaded in the section above. The name is the tumor you would like to analyze now. Here you can choose between "Tumor1" and "Tumor2".

```R
datasegment <- splitdata(data,name)
```

**Choose parameters**

The following parameters can be altered. If you do not want to change them, go with the default shown below.

```R
event_co <- 1000000 #Cutoff for events.
sub_co <- 0.9 #Cutoff for stem events.
eventnumber <- 300 #The maximum number of events I think a particular sample or subclone will have.
TDS <- "No" #"Yes", "No" or "Only". Choose wether or not you want TDS-events to be included in the following computations. If you choose yes we will not remove anything.
```

**The event matrix**

Let us now make an event matrix using the fabricated data set.

```R
EM_test <- DEVOLUTION(test,Tumorname="NB1",event_co=1000000,samples,eventnumber=300,TDS="No") #Creating the event matrix.
```
We now have the event matrix illustrating the subclones and which events each incorporates.

Show an example of the event matrix and what it means.

**Let's produce the final event matrix**

```R
EM_test_newnames <- simplify.tree(file_samples_subclones,EM_test,sample_clone_matrix)
```

What has the program done?

Table of sizes.

**Phylogenetic trees**

In the end one can use the event matrix in order to reconstruct phylogenetic trees using the maximum likelihood and parsimony method.

```R
EM_test_phy <- phydatevent(EM_test_newnames) #Transforming the EM to phyDat format.
EM_test_mptree <- mp_tree(EM_test_phy,root) #Constructing the maximum parsimony tree.
EM_test_mltree <- ml_tree(EM_test_phy,root) #Constructing the maximum likelihood tree.
limitmp <- xlim(c(0, 20)) #Here you can determine the limits for the graph for mp.
limitml <- xlim(c(0, 0.5)) #Here you can determine the limits for the graph for ml.
Treemp <- MP_treeplot(EM_test_mptree,limitmp) #Illustrating the maximum parsimony tree.
Treeml <- ML_treeplot(EM_test_mltree,limitml) #Illustrating the maximum likelihood tree.

w = 10
h = 10
ggsave(Treemp,filename= "Tree_mp.png",width = w,height = h)
ggsave(Treeml,filename= "Tree_ml.png",width = w,height = h)
```

Create the pie charts. They are saved to the working directory. The function yields a matrix illustrating the sizes of each subclone throughout the samples.

```R
clone_size <- pie.it(clonenames_new_order,root) #Pie charts are created and saved. You also get a list of all subclones and their sizes in each sample.
```

This yields a phylogenetic tree looking like this.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/NB7_pie_ml.png" width="600">
