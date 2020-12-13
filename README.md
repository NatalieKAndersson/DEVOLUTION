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
library("RColorBrewer") #Needed to add the colored pie charts.
```
If they are not installet you can install them by using the following command.

```R
install.packages(c("readxl","stringr","ape","phangorn","ggplot2",
                   "ggtree","ggimage","dplyr","RColorBrewer"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEVOLUTION")
```
Also load all functions by marking them and pressing “Run”.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Functions.PNG" width="400">

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

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Head_test.PNG" width="600">

It contains data from two tumors. The splitdata function determines the start and end position of these data sets. Here data is the whole file you loaded in the section above. The name is the tumor you would like to analyze now. Here you can choose between "Tumor1" and "Tumor2".

```R
datasegment <- splitdata(data,name)
```

**Choose parameters**

The following parameters can be altered. If you do not want to change them, go with the default shown below.

```R
datatypes <- "All" 
event_co <- 1000000
root <- "Normal"
```
- datatypes: These are your data types such as SNP-array, TDS, WGS, WES etc that you want to include in the analysis. If you write "All", all events will be included in the input file. If you exclusively want to analyze particular alterations you can set an input vector indicating which events should be kept for analysis. The command "#c(unique(test[,9]))" gives you the unique methods in your dataset.
- event_co: The cutoff for the start and end positions of the events in your segment file.
- root: In what cell you want to root your tree. Choose between "Normal" or "Stem".

**DEVOLUTION**

Let us now use DEVOLUTION to create an event matrix using the fabricated data set!

```R
EM <- DEVOLUTION(data,event_co,datatypes)
EM_dev <- subclones(EM,file_samples_subclones)
View(EM_dev[[1]])
```
We now have the event matrix illustrating the subclones and which events each incorporates. Here each row represents an identified subclone. The columns represent the genetic alterations found across the biopsies. The presence of a particular alteration in the subclone is represented by the number 1.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/EM_R.png" width="600">

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
