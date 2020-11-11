# DEVOLUTION
Phylogenetic reconstruction from multiregional sampling data

## Installation instructions

```R
install.packages(c("readxl","stringr","ape","phangorn","ggplot2",
                   "ggtree","ggimage","dplyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEVOLUTION")
```

If no error message has appeared, continue with the installation of the DEVOLUTION package. Install it in the folder in which your data files are (by replacint PATH with this) or make sure to set the path to that location after installation. You have to be in the folder of your files in order for the software to find them.

```R
install.packages("~PATH/DEVOLUTION_0.1.0.tar.gz", repos = NULL, type = "source")
```

We can now load the package to the global environment.

```R
require(DEVOLUTION)
require(splitdata)
require(simplify.tree)
require(pie.it)
require(phydat)
require(mp_tree)
require(ml_tree)
require(MP_tree)
require(ML_tree)
```

Try loading dataset "test". You should not get any error message doing this.

Remember to load the packages that the software depends upon.

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
## Usage

In order to illustrate the usage of the package, let's go through an example using the fabricated data set built into the package named "Segment". This is a made up data set representing a SNP-array output segment file.

**Load the data**

Let's load the data and take a look at the data itself.

```R
test <- load_matrix(filename="Segment.xlsx",sheetname ="Tumor1")
head(test)
```

What does the data set look like? Describe the columns.

It contains data from two tumors. The splitdata function determines the start and end position of these data sets.

```R
samples <- splitdata(test)
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

```R

BAF <- 0.57
N_AN <- 2 #Background.
N_BN <- 2 #Background.
N_A1 <- 1 #Subclone 1.
N_B1 <- 1 #Subclone 1.
N_A2 <- 1 #Subclone 2.
N_B2 <- 2 #Subclone 2.

#Two subclones + background population.
left <- BAF*N_AN+BAF*N_BN-N_BN
right_C1 <- N_B1-N_BN + BAF*N_AN + BAF*N_BN - BAF*N_A1 - BAF*N_B1
right_C2 <- N_B2-N_BN + BAF*N_AN + BAF*N_BN - BAF*N_A2 - BAF*N_B2

solutions <- matrix(0,3,100)

i <- 1
j <- 1
k <- 1
for(i in 1:100){
  
  C1 <- i/100
  
  for(j in 1:100){
    
    C2 <- j/100
    
    tot <- C1*right_C1 + C2*right_C2
    
    diff <- left-tot
    
    tot_pop <- C1+C2
    
    if(abs(diff) < 0.01 && tot_pop <= 1){
      print("Here is the difference, C1 and C2")
      print(tot)
      print(diff)
      print(tot_pop)
      print(C1)
      print(C2)
      
      solutions[1,k] <- C1
      solutions[2,k] <- C2
      solutions[3,k] <- 1-C1-C2
      k <- k+1
      
    }
    
    j <- j+1
  }
  
  
  i <- i+1
}

solutions_original <- solutions[,solutions[1,]!="0"]

x <- c(1:ncol(solutions_original))

solutions <- t(solutions_original)

solutions <- c(c(solutions[,1]),c(solutions[,2]),c(solutions[,3]))
solutions <- cbind(solutions,c(c(rep("C1",ncol(solutions_original)))
                        ,c(rep("C2",ncol(solutions_original)))
                        ,c(rep("Background",ncol(solutions_original)))
                        ))
solutions <- cbind(solutions,c(x,x,x))

solutions[,1] <- round(as.numeric(solutions[,1]),2)
df <- as.data.frame(solutions)


ggplot(df,aes(x=as.numeric(V3),y=as.numeric(solutions)))+geom_point(aes(color=V2),size = 3)+
  ylab("CCF")+xlab("")+scale_y_continuous(breaks=c(seq(0,1,0.1)))+theme(legend.title = element_blank())


```
