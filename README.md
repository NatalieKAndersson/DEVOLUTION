# DEVOLUTION <img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Logga_devolution.PNG" align = "right" width="180"/>
Devolution is an algorithm for phylogenetic reconstruction from multiregional sampling data that can incorporate information from SNP-array, WES, WGS, TDS etc. in unison or separately. It uses the mutated sample fraction as input, which is defined as the proportion of cancer cells in each biopsy that harbor a particular genetic alteration. This can be calculated from the VAF or log2-ratio.

<a href="https://zenodo.org/badge/latestdoi/297145258"><img src="https://zenodo.org/badge/297145258.svg" alt="DOI"></a>

## Setting up DEVOLUTION
Download the R script denoted "DEVOLUTION" and double click on the script to open it in your R-environment.

Before starting the analysis, we must load the dependencies of the algorithm. These are found in the beginning of the code.

```
library("readxl") #Needed to load the data from the xlsx file.
library("xlsx") #Needed to save matrices into xlsx-files.
library("stringr") #Needed for using the function "word".
library("ape")
library("phangorn") #Needed to transform the EM into phyDat and make trees.
library("ggplot2") #Needed to visualize the trees.
library("ggtree")
library("ggimage") #Needed to insert the pies in the tree.
library("dplyr") #Needed for the distinct function in pie.it.
library("RColorBrewer") #Needed to add the colored pie charts.
```
If they are not installed you can install them by using the following command.

```R
install.packages(c("readxl","xlsx","stringr","ape","phangorn","ggplot2",
                   "ggtree","ggimage","dplyr","RColorBrewer"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEVOLUTION")
```
Also load all functions by marking them and pressing “Run”.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Functions.PNG" width="400">

If no error message has appeared we are good to go!

## Usage
We are now ready to load some data and get going! At the top of the script you have to set the path i.e. the location of the files that will be analyzed.
```R
setwd("~/yourpath")
```
In order to illustrate the usage of the algorithm, let’s go through an example using the example data set that also can be downloaded from the github page. The data file is named “Segment.xlsx”. This is a fabricated data set representing a SNP-array output segment file for two tumors.

```R
data <- load_matrix(filename="Segment_example.xlsx",sheetname ="Example_tumors")
head(data)
```
What does the data set look like? Describe the columns. Make sure that you understand what each of the 11 columns represents in this data file.

- Tumor ID: The ID for the biopsies from the same patient.
- Sample ID: The biopsy names.
- Chr, Start, End: The position of the genetic alteration.
- Med.LogR: The median log2R obtained through chromosome analysis. Not required.
- VAF: Variant Allele Frequency. Not required.
- Method: The method by which this genetic aberration has been identified.
- Cytoband: The location of the genetic alteration or gene name.
- Clone size: The mutated sample fraction i.e. the fraction of cancer cells in this particular biopsy that has this aberration.

It contains data from two tumors. The splitdata function determines the start and end position of these data sets. Here data is the whole file you loaded in the section above. The name is the tumor you would like to analyze now. Here you can choose between "Tumor1" and "Tumor2".

```R
datasegment <- splitdata(data,name)
head(datasegment)
```

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Head_test.PNG" width="600">

**Choose parameters**

The following parameters can be altered. If you do not want to change them, go with the default shown below.

```R
datatypes <- "All" 
event_co <- 10*10^6
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
Running DEVOLUTION on the dataset "Tumor1" yields an "Execution time" Time difference of 0.242177 secs. We now have the event matrix illustrating the subclones and which events each incorporates. Here each row represents an identified subclone. The columns represent the genetic alterations found across the biopsies. The presence of a particular alteration in the subclone is represented by the number 1.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/EM_Tumor1.PNG" width="500">

We can also look at the matrix named "Clustering". This illustrates to what clusters each alteration is estimated to belong. This matrix is perfect for you to go through the data set yorself to see if there are any contradictions in the data set. Make sure you understand the connection between this table and the final phylogenetic trees.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Tumor1_cluster.PNG" width="500"> <img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Tumor1_allocation.PNG" width="400">


We can also look at the distribution of genetic alterations across the biopsies by writing the command.
```R
DB <- distribution(overview)

w = 10
h = 10
ggsave(DB,filename= "Distribution.png",width = w,height = h)
```
This is the information DEVOLUTION uses to infer the most probable evolutionary trajectory of the tumor.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Distribution_Tumor1.png" width="400">


**Phylogenetic trees**

Now we will create the phylogenetic trees based on the event matrix above. The following commands creates the trees without pie charts.

```R
EM_test_phy <- phydatevent(EM_test_newnames) #Transforming the EM to phyDat format.

EM_test_mptree <- mp_tree(EM_test_phy,root) #Constructing the maximum parsimony tree.
limitmp <- xlim(c(0, 20)) #Here you can determine the limits for the graph for mp. Alter so that the tree fits in the plot viewer.
Treemp <- MP_treeplot(EM_test_mptree,limitmp) #Illustrating the maximum parsimony tree.

EM_test_mltree <- ml_tree(EM_test_phy,root) #Constructing the maximum likelihood tree.
limitml <- xlim(c(0, 0.5)) #Here you can determine the limits for the graph for ml. Alter so that the tree fits in the plot viewer.
Treeml <- ML_treeplot(EM_test_mltree,limitml) #Illustrating the maximum likelihood tree.

w = 10
h = 10
ggsave(Treemp,filename= "Tree_mp.png",width = w,height = h)
ggsave(Treeml,filename= "Tree_ml.png",width = w,height = h)
```

Create the pie charts. They are saved to the working directory. The function yields a matrix illustrating the sizes of each subclone throughout the samples. The pie charts are also saved to your computer as individual files. They will be replaced as you do a new tree if not saved elsewhere.

```R
pieData <- make_pie(EM_dev[[2]],root,samples,type="nocol") #Creates the pie charts.
pietree <- pie_it(Treemp,pieData,offset=1.5,size=0.21,col="no") #Adds pie charts to the tree.
```

- Type: You can choose which type of pies you want. Choose between
    - "col": Colored pie charts are used. Random order.
    - "nocol": Red pie charts.
    - "custom": You can use a vector of colors of your own as an input. The colors will be given from the color vector you give in the same order as in the samples vector also given to the pie.tree function.

The final trees look like this with the "nocol" and "col" setup. You can of course also use ggplot to alter the plot as you wish. The trees, pie charts and legend can all be saved separately.

Compare the trees to the event matrices that we discussed in the previous section!

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Tumor1_mp_pie_nocol.png" width="400"><img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Tumor1_mp_pie.png" width="400">

**Adding a rule**

The goal is to incorporate user-controlled rules for avoiding imposition of illicit biological trajectories Some genetic aberrations present in the data set might be known to never occur in the same cell for some well-known biological reason. Such constraints should optimally be supplied to the algorithm to ensure biologically plausible solutions. The user can therefore provide the DEVOLUTION algorithm with a matrix indicating which genetic aberrations in the data set that cannot be placed after one another. The subclonal deconvolution algorithm extracts a list for each genetic alteration containing information about in how many of the samples it can be allocated after a certain cluster. There might be multiple possible solutions, equally prevalent. In this instance the matrix containing information about illicit biological orders can aid the program in taking a decision regarding which of these allocations are less likely, subsequently discarding them. These rules will thus only be employed if the data set allows the genetic alterations to be placed in any other way. If the only possible way for the events to be allocated is to place them as descendants, the user will be advised to revise the original data set.

In the file "Segment.xlsx" there is an example of such a case in the sheet "Example_rule". In this sample 50 % of the cells have a loss of one copy of 17p13q21 which results in the allelic composition 1+0 (loss of heterozygozity = LOH) for this segment. In the biopsy 30 % of the cells have gained a copy of this segment hence having the allelic composition 2+1. A cell who has a LOH of a segment can never return to a heterozygous state. Hence, this evolutionary order of events is biologically unlikely.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Chromosomes_lossandgain.PNG" width="600">

The rule matrix should have the following structure where the first column is the mother event and the second one the daughter event it cannot have.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Rule_matrix.PNG" width="200">

Here the tree can be seen before (left) and after (right) using the rule.

**Alternative solutions**

There may be situations where multiple phylogenetic trees are able to explain the observed data. Therefore, an algorithm was constructed to obtain these alternative solutions. In the case where there is more than one solution, the user will be asked if the suggested solution by DEVOLUTION should be provided or an alternative solution.

In that case the clusters of genetic alterations unique for those subclones that have multiple solutions are removed from the tree structure and randomly reshuffled to produce a new phylogenetic tree, that does not confer any of the rules in any of the samples or rules provided by the user. The user is also provided with a matrix illustrating which subclones in the tree are reliable and which are uncertain due to multiple possible evolutionary trajectories.

An example. Imagine that you have the following data set.
<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Segment_file_multiple.PNG" width="600">

In this case the following event matrix is obtained from DEVOLUTION.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/EM_multiple_final.PNG" width="400">

The genetic alteration that distinguish subclone C from the other subclones have multiple solutions for its nesting.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Certainty_matrix.PNG" width="150">

The user is therefore asked if the suggested tree should be shown or an alternative solution. Here you can see the suggested tree (left) and an alternative tree (right). The user can also choose to color the subclone names for which only one solution is possible (red) or multiple solutions are possible (green). The colors can also easily be changed.

<img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Tumor_solution1_MP.png" width="400"><img src="https://github.com/NatalieKAndersson/DEVOLUTION/blob/master/Tumor_solution2_MP.png" width="400">
