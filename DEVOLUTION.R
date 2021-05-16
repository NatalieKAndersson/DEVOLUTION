###############################################################################
#-------------------------------DEVOLUTION------------------------------------#
###############################################################################
setwd("~/yourdirectory")
#Dependencies----
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
library("ggridges") #Used to plot the distribution.
library("cowplot")
library("dbscan") #Clustering

#Start with pressing the little triangle to the left here to collapse all functions!

#Functions----
#Function extracting the data
load_matrix <- function(filename, sheetname) {
  data <- as.data.frame(read_xlsx(filename, sheetname)) #Reading the xlsx file and saving it in the variable data.
  subdata <- data[ c(1:nrow(data)), c(1:ncol(data)) ] #Extracting the part of the file we are interested in.
  subdata <- subdata[is.na(subdata[,1])==FALSE,]
  return(subdata)
}

#Function creating the eventmatrix
DEVOLUTION <- function(file,eventcutoff,datatypes, rule, eps){
  start.time <- Sys.time()
  
  if(missing(eventcutoff)==TRUE){
    print("You have not chosen an event cutoff. Default 1 Mbp chosen.")
    eventcutoff <- 1000000
  }
  
  all_cols <- ncol(file)
  if(all_cols < 11){
    print("There are missing columns!")
    types <- c("Tumor ID","Samples","Chr","Start","End","Med LogR","VAF (TRS)","Type","Method","Cytoband/ Gene","Clone size (%)")
    thematch <- match(colnames(file),types)
    missing <- types[types%in%colnames(file)==FALSE]
    print("This is missing. The algorithm will add it.")
    print(missing)
    
    file_new <- matrix(0,nrow(file),11)
    colnames(file_new) <- types
    
    cols_file <- colnames(file)
    i <- 1
    for(i in 1:11){
      
      if(types[i]%in%missing){
        file_new[,i] <- "NA"
      }else{
        col <- match(types[i],cols_file)
        file_new[,i] <- file[,col]
      }
      i <- i+1
    }
    file <- file_new 
  }
  ################################
  #Treating NA for TC in the file#
  ################################
  i <- 1
  for(i in 1:as.numeric(nrow(file))){
    if(file[i,11] == "NA"){
      if(file[i,2] == "ALL"){
        file[i,11] <- "100"
      }
      if(file[i,2] != "ALL"){
        file[i,] <- "0"
      }
    }
    i <- i+1
  }
  
  ###################################################################
  #Removing events obtained with a method not specified by datatypes#
  ###################################################################
  #If you choose yes we will not remove anything.
  if(length(datatypes) > 1){
    print("Only the following datatypes are included in the analysis.")
    print(datatypes)
    i <- 1
    for(i in 1:as.numeric(nrow(file))){
      
      if(file[i,9] %in% datatypes == FALSE){
        file[i,] <- "0"
      }
      
      i <- i+1
    }
  }else{
    if(datatypes == "All"){
      print("All datatypes supplied in data are included in the analysis.")
    }else{
      print("The following datatype is included in the analysis")
      print(datatypes)
      i <- 1
      for(i in 1:as.numeric(nrow(file))){
        
        if(file[i,9] %in% datatypes == FALSE){
          file[i,] <- "0"
        }
        
        i <- i+1
      }
    }
  }
  
  file <- as.matrix(file[file[,2] != "0",])
  ###########################
  #Finding all unique events#
  ###########################
  #This loop defines all events that we have in the dataset for a particular tumor.
  #The events have to have the same name, be on the same chromosome and have breakpoints
  #within a certain cutoff that is set on beforehand.
  
  file_new <- file
  i <-  1
  versionnames <- c(paste( c("v"), 1:30, sep=""))
  k <- 1
  for(i in 1:nrow(file_new)){ #Choosing a row.
    j <- 1
    for(j in 1:nrow(file_new)){ #Choosing another row to compare it with.
      if(i != j){
        
        if(file_new[i,10] == file_new[j,10]){ #Comparing event names.
          if(file_new[i,8] == file_new[j,8]){ #Comparing the type of event.
            if(file_new[i,3]==file_new[j,3]){ #Are they on the same chromosome?
              if(file_new[i,2] == "ALL"){ #If the event is a part of the stem they shall always be separate.
                if(file_new[j,2] == "ALL"){
                  if(is.na(word(file[i,10],2)) == TRUE){ #If the one you compare with does not already have a version name.
                    file_new[j,10] <- paste(file_new[j,10],versionnames[k], sep = "") #Changing the name for the second version of the mutation.
                    k <- k+1
                  }else if(is.na(word(file[i,10],2)) == FALSE){
                    versionpos <- match(word(file[i,10],2),versionnames)
                    newversion <- versionpos+1
                    file_new[j,10] <- paste(file_new[j,10],versionnames[newversion], sep = "") #Changing the name for the second version of the mutation.
                  }
                }}
              if(file_new[i,2] != "ALL"){ #This part is only valid for non stem events.
                if(abs(as.numeric(file_new[i,4]) - as.numeric(file_new[j,4])) > eventcutoff){ #If the events differ too much in genetic distance they are seen as two separate events.
                  file_new[j,10] <- paste(file_new[j,10],"v1", sep = "")} #Changing the name for the second version of the mutation.
                else if(abs(as.numeric(file_new[i,5])-as.numeric(file_new[j,5])) > eventcutoff){ #The same but in the other direction.
                  file_new[j,10] <- paste(file_new[j,10],"v1", sep = "")}
              }
            }
          }
        }
      }
      j <- j+1
    }
    i <- i+1
  }
  
  for(i in 1:nrow(file_new)){ #Adding information to the events about which kind of alteration it is.
    file_new[i,10] <- paste(file_new[i,10],file_new[i,8],file_new[i,3], sep = " ")
    i <- i+1
  }
  
  un_s <- unique(file_new[,2])
  samples <- un_s[un_s!="ALL"]
  ###########################
  #Making an overview matrix#
  ###########################
  samples <- as.matrix(unique(file_new[,2])) #Extracting all unique samples.
  aberrations <- as.matrix(unique(file_new[,10])) #Extracting all unique events.
  #Constructing a matrix with all samples and their TC for each of the unique events.
  overview <- matrix(0,(length(aberrations)+1),(length(samples)+1))
  overview[1,2:as.numeric(ncol(overview))] <- samples
  overview[2:as.numeric(nrow(overview)),1] <- aberrations
  
  i <- 1
  for(i in 1:nrow(file_new)){ #Extracting all of the TC:s.
    
    samplepos <- match(file_new[i,2],overview[1,])
    aberrationpos <- match(file_new[i,10],overview[,1])
    overview[aberrationpos,samplepos] <- file_new[i,11]
    
    if(file_new[i,2] == "ALL"){
      overview[aberrationpos,2:ncol(overview)] <- 100 #All samples should have 100 on the "ALL" events.
    }
    
    i <- i+1
  }
  #Do we have any stem at all?
  if(overview[1,2] != "ALL"){
    print("You do not have any declared stem event denoted ALL.")
    
    allcolumn <- matrix(0,nrow(overview),1)
    allcolumn[1,1] <- "ALL"
    overview <- cbind(overview[,1],allcolumn,overview[,2:ncol(overview)])
    
  }
  
  #Treating cases where not all stem events have been declared.
  i <- 2
  firststem <- 1
  for(i in 2:nrow(overview)){
    
    stemornot <- (length(which(as.numeric(overview[i,2:ncol(overview)])>=90))-(as.numeric(ncol(overview))-2))
    
    if(stemornot == 0){
      #This is a stem event.
      overview[i,2:ncol(overview)] <- 100 #Declaring it as a stem event.
      
      #Now we have to declare it a stem in the file as well and remove it from the other ones.
      if(firststem == 1){
        row <- match(overview[i,1],file_new[,10])
        stemmatrix <- t(as.matrix(file_new[row,]))
        stemmatrix[1,2] <- "ALL"
        stemmatrix[1,11] <- "100"
        #stemmatrix[1,1] <- "Remove" #Temporary.
        firststem <- 2
        
        pos_stem <- as.numeric(which(overview[i,1]==file_new[,10])) #The positions in which the stem exists.
        file_new[pos_stem,1] <- "Remove"
        
      }else{
        row <- match(overview[i,1],file_new[,10])
        event <- t(as.matrix(file_new[row,]))
        event[1,2] <- "ALL"
        event[1,11] <- "100"
        #event[1,1] <- "Remove" #Temporary.
        stemmatrix <- rbind(stemmatrix,event)
        
        pos_stem <- as.numeric(which(overview[i,1]==file_new[,10])) #The positions in which the stem exists.
        file_new[pos_stem,1] <- "Remove"
      }
      
      
    }
    
    
    if(i == nrow(overview)){
      if(firststem == 1 && as.numeric(overview[2,2]) == 0){
        print("There is no stem events in the data. Adding a fabricated stem.")
        f_stem <- matrix(0,1,11)
        f_stem[1,] <- c(unique(file_new[file_new[,1]!="Remove",1]),"ALL",1,1,1,"NA","NA","Stem","WES","Stem","100")
        file_new <- rbind(f_stem,file_new)
      }else{
        if(firststem!=1){
          file_new <- rbind(stemmatrix,file_new)}
      }
    }
    
    
    i <- i+1
  }
  
  # #We want to order the overview a bit.
  overview_new <- matrix(0,nrow(overview),ncol(overview))
  overview_new[1,] <- overview[1,]
  sub <- overview[2:nrow(overview),]
  ov_stem <- sub[as.numeric(sub[,2])==100,]
  ov_notstem <- sub[as.numeric(sub[,2])!=100,]
  overview_new[2:nrow(overview_new),] <- rbind(ov_stem,ov_notstem)
  overview <- overview_new
  
  assign("file_new_stem", file_new, envir=globalenv())
  file_new <- file_new[file_new[,1]!="Remove",]
  assign("overview_stem", overview, envir=globalenv())
  assign("file_new_removed", file_new, envir=globalenv())
  
  if(firststem == 2){
    assign("stemmatrix", stemmatrix, envir=globalenv())}
  
  ###################################################################
  #Using Density-Based Spatial Clustering of Applications with Noise#
  ###################################################################
  #overview_truncated <- overview[2:nrow(overview),2:ncol(overview)]
  sub <- overview[2:nrow(overview),]
  overview_truncated <- sub[as.numeric(sub[,2])!=100,] #Removing the ALL events so that events are not clustered into the stem.
  oneevent <- 0
  assign("overview_test", overview, envir=globalenv())
  if(is.null(dim(overview_truncated))==FALSE){
    print("more")
    overview_truncated <- overview_truncated[,3:ncol(overview_truncated)]
    overview_truncated <- as.data.frame(overview_truncated)
  }else{
    print("only one")
    overview_truncated <- overview[2:nrow(overview),3:ncol(overview)]
    overview_truncated <- as.data.frame(overview_truncated)
    oneevent <- 1
  }
  overview_dfm <- data.matrix(overview_truncated, rownames.force = NA)
  
  if(missing(eps)==TRUE){
    eps <- 0.5
  }
  
  library(dbscan)
  x <- kNNdist(overview_dfm, k = 1)
  kNNdistplot(overview_dfm,k=1)
  abline(h=eps, col = "red", lty=2)
  
  myclusters <- dbscan(overview_dfm, eps = eps ,minPts = 1)
  #myclusters <- dbscan(overview_dfm, eps = 15 ,minPts = 1) #TRACERx
  
  print(myclusters)
  assign("myclusters", myclusters, envir=globalenv())
  
  
  # View(overview_dfm)
  # fviz_cluster(myclusters,data = overview_dfm, minPts = 1) #Plotting the clusters.
  #
  #  if(length(unique(overview_dfm)) >= 2){
  #  fviz_cluster(myclusters,data = overview_dfm, minPts = 1) #Plotting the clusters.
  #  }else{
  #    x = "Nope"
  #    print("Warning message: The input matrix does only contain one single subclone.")
  #    stopifnot(Datadimensions == "ok")
  #  }
  #assign("overview", overview, envir=globalenv())
  #######################################################################
  #Constructing a matrix indicating which events belong to which cluster#
  #######################################################################
  #If we use DBSCAN
  clusters <- as.matrix(myclusters$cluster)
  overview_new <- cbind(overview,matrix(0,nrow(overview),1)) #Adding the cluster belonging to the overview.
  
  if(oneevent == 0){
    if(is.null(nrow(ov_stem))==TRUE){
      t <- 1
    }else{
      t <- nrow(ov_stem)
    }
    overview_new[2:nrow(overview_new),ncol(overview_new)] <- c(c(rep("ALL",t)),clusters[,1])
    unique_clusters <- c(c("ALL"),c(unique(clusters[,1]))) #Constructing the matrix in which the events will be saved.
  }else{
    overview_new[2:nrow(overview_new),ncol(overview_new)] <- c(clusters[,1])
    unique_clusters <- c(unique(clusters[,1])) #Constructing the matrix in which the events will be saved.
  }
  assign("overview_cluster",overview_new,envir=globalenv())
  cluster_matrix <- matrix(0,as.numeric(length(unique_clusters)),200)
  print(unique_clusters)
  
  
  i <- 1
  for(i in 1:length(unique_clusters)){
    if(unique_clusters[i]== "ALL"){
      cluster_matrix[i,1] <- "ALL"
    }else{
      cluster_matrix[i,1] <- paste("Subclone",unique_clusters[i])}
    i <- i+1
  }
  
  i <- 1
  for(i in 2:nrow(overview_new)){ #Looping through the subclonal belonging.
    j <- 2
    
    Subclonerow <- match(overview_new[i,ncol(overview_new)],unique_clusters)
    for(j in 2:ncol(cluster_matrix)){ #Looping through the available spots for saving the event.
      if(cluster_matrix[Subclonerow,j] == "0"){
        cluster_matrix[Subclonerow,j] <- overview_new[i,1]
        break
      }
      j <- j+1
    }
    i <- i+1
  }
  
  addmatrix <- matrix(0,as.numeric(nrow(cluster_matrix)),3)
  addmatrix[,1] <- cluster_matrix[,1]
  
  clone_matrix_names <- cbind(addmatrix,cluster_matrix[,2:as.numeric(ncol(cluster_matrix))]) #Adding three columns in the beginning.
  #View(clone_matrix_names)
  #########################
  #Extracting the clusters#
  #########################
  #Extracting clusters, calculating the median TC for each one and assigning the
  #cluster names to the events in the segment file.
  clusterTC <- matrix(0,200,3)
  calculateTC <- matrix(0,400,1)
  i <- 1
  s <- 1
  t <- 1
  #View(file_new)
  for(i in 1:as.numeric(nrow(clone_matrix_names))){ #Looping through the rows in the clustermatrix.
    
    j <- 1
    for(j in 1:as.numeric(ncol(clone_matrix_names))){ #Looping through the mutations in a particular cluster i.
      
      if(clone_matrix_names[i,j] != 0 && clone_matrix_names[i,j] != "0"){
        
        k <- 1
        for(k in 1:nrow(file_new)){ #Looping through the mutations in our datafile.
          if(clone_matrix_names[i,j] == file_new[k,10]){
            calculateTC[s,1] <- file_new[k,11] #Saving the TC
            s <- s+1
          }
          k <- k+1
        }
        
      }
      j <- j+1
    }
    #We now have a vector with all of the TC:s for that cluster.
    clusterTC[t,1] <- paste("Cluster",t)
    clusterTC[t,2] <- median(as.numeric(calculateTC[calculateTC[,1] != 0,]))
    clone_matrix_names[i,2] <- median(as.numeric(calculateTC[calculateTC[,1] != 0,]))
    calculateTC <- matrix(0,400,1) #Resetting the matrix.
    s <- 1
    t <- t + 1
    
    i <- i+1
  }
  
  
  clusterTC_order <- clusterTC[order(as.numeric(clusterTC[,2]), decreasing = TRUE),]
  clone_matrix_names <- clone_matrix_names[order(as.numeric(clone_matrix_names[,2]), decreasing = TRUE),]
  
  i <- 1
  # #If we do not have any ALL-events we want to use another name vector.
  if("ALL" %in% file_new[,2] == FALSE){
    namevector <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z","ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZX","ZY","ZZ")
  }else{
    namevector <- c("ALL","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z","ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZX","ZY","ZZ")
  }
  
  for(i in 1:as.numeric(nrow(clusterTC_order))){
    
    if(as.numeric(clusterTC_order[i,2]) != 0){
      if(i == 1){
        if("ALL" %in% file_new[,2] == TRUE){ #If we have an ALL event we do not want to add the name "Subclone" to it.
          clusterTC_order[i,1] <- paste("Cluster",i)
          clusterTC_order[i,3] <- paste(namevector[i])
          clone_matrix_names[i,1] <- paste("Cluster",i)
          clone_matrix_names[i,3] <- paste(namevector[i])
        }else{
          clusterTC_order[i,1] <- paste("Cluster",i) #In the case where we do not have any ALL-events we want to add "Subclone" to the first subclone.
          clusterTC_order[i,3] <- paste("Subclone_",namevector[i])
          clone_matrix_names[i,1] <- paste("Cluster",i)
          clone_matrix_names[i,3] <- paste("Subclone_",namevector[i])
        }
      }else{
        clusterTC_order[i,1] <- paste("Cluster",i)
        clusterTC_order[i,3] <- paste("Subclone_",namevector[i])
        clone_matrix_names[i,1] <- paste("Cluster",i)
        clone_matrix_names[i,3] <- paste("Subclone_",namevector[i])
      }
    }
    i <- i+1
  }
  assign("clone_matrix_names_hej", t(clone_matrix_names), envir=globalenv()) #The subclone names and which mutations are included in each subclone are exported to the global environment.
  ##################################################################
  #Making a new file indicating which subclone each event belong to#
  ##################################################################
  file_original <- file_new
  file_subclones <- file_new
  
  i <- 1
  j <- 1
  k <- 1
  for(i in 1:as.numeric(nrow(clone_matrix_names))){ #Looping through a certain subclone.
    for(j in 4:as.numeric(ncol(clone_matrix_names))){ #Looping through the clone mutations.
      
      for(k in 1:as.numeric(nrow(file_new))){ #Looping through the file.
        
        if(clone_matrix_names[i,j] == file_new[k,10]){ #If we match a subclone event with an event in our data set.
          file_subclones[k,2] <- clone_matrix_names[i,3]
        }
        
        k <- k+1
      }
      
      j <- j+1
    }
    i <- i+1
  }
  
  ##################################################
  #Finding out which subclones exist in each sample#
  ##################################################
  samples <- as.matrix(unique(c(file_original[,2]))) #Extracting all unique samples.
  sample_clone_matrix <- matrix(0,40,(as.numeric(length(samples))*3)) #A matrix which will contain all of the samples and their clones. One column is used in order to asses the number of mutations that exist within the subclone.
  file_samples_subclones <- cbind(file_new,file_subclones[,2]) #Adding a column with all of the subclonal denotations of the alterations.
  file_samples_subclones <- cbind(file_samples_subclones,matrix(0,as.numeric(nrow(file_samples_subclones)),1))
  
  i <- 1
  m <- 1
  for(i in 1:as.numeric(nrow(file_samples_subclones))){ #Adding the subclonal names to the sample names.
    if(file_samples_subclones[i,2] != "ALL"){
      file_samples_subclones[i,13] <- paste(file_samples_subclones[i,2],file_samples_subclones[i,12])
    }else{
      file_samples_subclones[i,13] <- "ALL"
    }
    i <- i+1
  }
  
  subclones <- as.matrix(unique(c(file_samples_subclones[,13]))) #Extracting all unique subclones within samples.
  medianmatrix <- matrix(0,100,as.numeric(length(subclones))) #This matrix is to be used in order to calculate the median TC for each subclone within each sample.
  medianmatrix[1,] <- subclones #The first row consists of the subclone names.
  samples_unique <- as.matrix(unique(c(file_samples_subclones[,2]))) #Extracting all unique samples.
  
  i <- 1
  if("ALL" %in% file_new[,2] == TRUE){
    for(i in 1:(ncol(sample_clone_matrix)/3)){ #All samples have the subclone named "ALL".
      sample_clone_matrix[2,(3*i-2)] <- "ALL"
      sample_clone_matrix[2,(3*i-1)] <- "100"
      sample_clone_matrix[1,(3*i-2)] <- samples_unique[i,1]
      i <- i+1
      s <- 3
    }
  }else{
    s <- 2
  }
  
  i <- 1
  t <- 1
  for(i in 1:as.numeric(nrow(file_samples_subclones))){ #Looping through the dataset.
    
    if(i == 1){ #The first position will of course be our first sample in the matrix.
      
      if("ALL" %in% file_new[,2] == TRUE){
        sample_clone_matrix[1,1] <- file_samples_subclones[1,2] #Name.
      }else{
        sample_clone_matrix[1,1] <- file_samples_subclones[1,2] #When we do not have ALL-events in the tumor we need to extract the subclone name.
        sample_clone_matrix[s,1] <- file_samples_subclones[1,12] #Name.
      }
      
      if(file_samples_subclones[1,2] != "ALL"){
        sample_clone_matrix[s,2] <- file_samples_subclones[1,11] #TC.
        sample_clone_matrix[s,3] <- 1
        s <- s+1
        
        medianmatrix[1,1] <- file_samples_subclones[1,13] #Name.
        medianmatrix[2,1] <- file_samples_subclones[1,11] #TC.
      }
      if(file_samples_subclones[i,2] != file_samples_subclones[i+1,2]){
        t <-t+1
      }
      
    }
    
    if(i != 1){
      if(i < as.numeric(nrow(file_samples_subclones))){
        
        if(file_samples_subclones[i,2] == file_samples_subclones[i+1,2]){ #We are still within the same sample.
          if((file_samples_subclones[i,12] %in% sample_clone_matrix[,((3*t)-2)]) == FALSE){ #The event should not already be in that column. We just want the unique labels.
            sample_clone_matrix[1,(3*t-2)] <- file_samples_subclones[i,2] #Sample name.
            sample_clone_matrix[s,(3*t-2)] <- file_samples_subclones[i,12] #Saving the clone name.
            sample_clone_matrix[s,(3*t-1)] <- file_samples_subclones[i,11] #clone_matrix_names[as.numeric(match(file_samples_subclones[i,12], clone_matrix_names[,3])),2] #Saving the clone TC in that sample.
            sample_clone_matrix[s,(3*t)] <- (as.numeric(sample_clone_matrix[s,(3*t)])+1)
            #sample_clone_matrix[s,2*t] <- file_samples_subclones[i,11] #Saving the clone TC in that sample.
            s <- s+1
          }else{ #If the event already is in that column we add the TC.
            if(file_samples_subclones[i,2] != "ALL"){
              matchrow <- match(file_samples_subclones[i,12],sample_clone_matrix[,((3*t)-2)])
            }else{matchrow <- 2}
            sample_clone_matrix[matchrow,(3*t-1)] <- (as.numeric(sample_clone_matrix[matchrow,(3*t-1)]) + as.numeric(file_samples_subclones[i,11])) #Saving the clone TC in that sample.
            sample_clone_matrix[matchrow,(3*t)] <- (as.numeric(sample_clone_matrix[matchrow,(3*t)])+1)
          }
        }
        
        
        if(file_samples_subclones[i,2] != file_samples_subclones[i+1,2]){ #New sample.
          
          if(file_samples_subclones[i,12] %in% sample_clone_matrix[,((3*t)-2)] == FALSE){ #The event should not already be in that column. We just want the unique labels.
            sample_clone_matrix[1,(3*t-2)] <- file_samples_subclones[i,2] #Sample name.
            sample_clone_matrix[s,(3*t-2)] <- file_samples_subclones[i,12] #Saving the clone name.
            sample_clone_matrix[s,(3*t-1)] <- file_samples_subclones[i,11] #clone_matrix_names[as.numeric(match(file_samples_subclones[i,12], clone_matrix_names[,3])),2] #Saving the clone TC in that sample.
            sample_clone_matrix[s,(3*t)] <- (as.numeric(sample_clone_matrix[s,(3*t)]) + 1) #Counting events.
            #sample_clone_matrix[s,2*t] <- file_samples_subclones[i,11] #Saving the clone TC in that sample.
          }else{
            if(file_samples_subclones[i,2] != "ALL"){
              matchrow <- match(file_samples_subclones[i,12],sample_clone_matrix[,((3*t)-2)])
            }else{matchrow <- 2}
            sample_clone_matrix[matchrow,(3*t-1)] <- (as.numeric(sample_clone_matrix[matchrow,((3*t)-1)]) + as.numeric(file_samples_subclones[i,11])) #Saving the clone TC in that sample.
            sample_clone_matrix[matchrow,(3*t)] <- (as.numeric(sample_clone_matrix[as.numeric(matchrow),3*t])+1)
          }
          s <- 3 #Resetting the s.
          if(t != length(samples)){
            t <- t+1 #Going to the next triad of columns.
          }
        }
        
      }
    }
    
    if(i == as.numeric(nrow(file_samples_subclones))){ #We're at the end of the file.
      if(file_samples_subclones[i,2] != file_samples_subclones[i-1,2]){ #If the last row actually is a new sample.
        s <- 3 #Resetting the s.
        if(file_samples_subclones[i,12] %in% sample_clone_matrix[,((3*t)-2)] == FALSE){ #The event should not already be in that column. We just want the unique labels.
          sample_clone_matrix[1,(3*t-2)] <- file_samples_subclones[i,2] #Sample name
          sample_clone_matrix[s,(3*t-2)] <- file_samples_subclones[i,12] #Saving the clone name.
          sample_clone_matrix[s,(3*t-1)] <- file_samples_subclones[i,11] #clone_matrix_names[as.numeric(match(file_samples_subclones[i,12], clone_matrix_names[,3])),2] #Saving the clone TC in that sample.
          #sample_clone_matrix[s,2*t] <- file_samples_subclones[i,11] #Saving the clone TC in that sample.
        }
      }
      if(file_samples_subclones[i,2] == file_samples_subclones[i-1,2]){ #If the last row is the same sample.
        
        if(file_samples_subclones[i,12] %in% sample_clone_matrix[,((3*t)-2)] == FALSE){
          sample_clone_matrix[s,(3*t-2)] <- file_samples_subclones[i,12] #Saving the clone name.
          sample_clone_matrix[s,(3*t-1)] <- file_samples_subclones[i,11] #clone_matrix_names[as.numeric(match(file_samples_subclones[i,12], clone_matrix_names[,3])),2] #Saving the clone TC in that sample.
          sample_clone_matrix[s,(3*t)] <- (as.numeric(sample_clone_matrix[s,(3*t)]) + 1) #Counting events.
        }else{
          matchrow <- match(file_samples_subclones[i,12],sample_clone_matrix[,(3*t-2)])
          sample_clone_matrix[matchrow,(3*t-1)] <- (as.numeric(sample_clone_matrix[matchrow,(3*t-1)]) + as.numeric(file_samples_subclones[i,11])) #Saving the clone TC in that sample.
          sample_clone_matrix[matchrow,(3*t)] <- (as.numeric(sample_clone_matrix[as.numeric(matchrow),3*t])+1)
        }
      }
      
    }
    
    samplematch <- match(file_samples_subclones[i,13],medianmatrix[1,])
    
    k <- 2
    m <- 0
    for(k in 2:as.numeric(nrow(medianmatrix))){ #Making a matrix for calculating the median.
      
      if(medianmatrix[k,samplematch] == "0"){
        if(m == 0){
          medianmatrix[k,samplematch] <- file_samples_subclones[i,11]
          m <- 1
        }
      }
      
      k <- k+1
    }
    
    
    i <- i+1
  }
  
  i <- 1
  for(i in 1:as.numeric(ncol(medianmatrix))){
    column <- as.matrix(medianmatrix[2:nrow(medianmatrix),i])
    column <- column[column[,1] != "0",1]
    medianmatrix[as.numeric(nrow(medianmatrix)),i] <- median(as.numeric(column))
    i <- i+1
  }
  
  #Adding the TC to the matrix illustrating the subclonal architecture within a sample.
  i <- 1
  for(i in 1:as.numeric(ncol(medianmatrix))){ #Looping through the subclones.
    
    columnsample <- match(word(medianmatrix[1,i],1),sample_clone_matrix[1,]) #Locating the sample. We get the column for the sample in sample_clone_matrix.
    
    if(medianmatrix[1,i] != "ALL"){
      rowsample <- match(word(medianmatrix[1,i],2,3),sample_clone_matrix[,columnsample]) #Locating the subclone in the row for the sample.
    }else{rowsample <- match(word(medianmatrix[1,i],1),sample_clone_matrix[,columnsample])}
    
    if(is.na(rowsample) == FALSE){
      if(columnsample != 1){
        sample_clone_matrix[rowsample,(columnsample+1)] <- medianmatrix[nrow(medianmatrix),i] #Adding the median TC.
      }else{sample_clone_matrix[2,2] <- "100"}
    }
    i <- i+1
  }
  
  assign("sample_clone_matrix", sample_clone_matrix, envir=globalenv()) #The matrix which tells us which mutations belong to which subclone is transferred to the global environment.
  
  ####################################
  #Building the event matrix skeleton#
  ####################################
  subclones <- as.matrix(unique(c(file_samples_subclones[,13]))) #Extracting all unique subclones.
  samples <- as.matrix(unique(c(file_new[,2]))) #Extracting all unique samples.
  events <- as.matrix(unique(c(file_new[,10]))) #Extracting all unique events.
  
  EMc <- nrow(subclones)+1
  EMr <- nrow(events)+1
  eventmatrix <- matrix(0,EMr,EMc) #Creating an empty event matrix.
  eventmatrix[1,2:EMc] <- subclones #The subclone names are placed on the firs row of the event matrix.
  eventmatrix[2:EMr,1] <- events #The event names are placed in the first column of the event matrix.
  eventnumber <- nrow(file_new) #The unpper bound of events we think the tumor will have.
  events <- matrix(0,eventnumber,as.numeric(nrow(subclones))) #Creating an empty matrix for the events belonging to each subclone.
  events[1,] <- subclones
  #########################################################################################################
  #Allocating the events to the samples/subclones. All subclones should have the events belonging to "ALL"#
  #########################################################################################################
  i <- 1
  for(i in 1:ncol(events)){ #Looping through every subclone separately.
    j = 1
    s = 2
    for(j in 1:nrow(file_samples_subclones)){ #Going through all of the events for the data set.
      if(file_samples_subclones[j,13] == "ALL"){ #If we find an "ALL"-event the sample should always have this one.
        events[s,i] <- file_samples_subclones[j,10]
        s <- s+1
      }
      else if(events[1,i] == file_samples_subclones[j,13]){ #If we find an event belonging to the subclone we add it to the EM.
        if((file_samples_subclones[j,10] %in% events[,i]) == FALSE){
          events[s,i] <- file_samples_subclones[j,10]
          s <- s+1
        }
      }
      j <- j+1
    }
    i <- i+1
  }
  #############################
  #Adding the events to the EM#
  #############################
  i <- 1
  for(i in 2:nrow(eventmatrix)){ #Events in the EM.
    j <- 1
    for(j in 1:ncol(events)){ #Cells.
      if(eventmatrix[i,1] %in% events[,j] == TRUE){ #Check if the events exist in this cell.
        eventmatrix[i,j+1] <- 1
      }
      else(eventmatrix[i,j+1] <- 0)
    }
    i <- i+1
  }
  ###############################################################
  #The subclones should have the events that its motherclone has#
  ###############################################################
  i <- 2
  j <- 1
  s <- 3
  space <- matrix(0,50,1)
  
  #The events of the mother clones are allocated within a single sample.
  i <- 1
  j <- 1
  s <- 2
  t <- 1
  
  space <- matrix(0,50,2) #Spaces within a sample. Dynamic.
  totalspace <- matrix(0,(as.numeric(nrow(space)+1)),((2*as.numeric(ncol(sample_clone_matrix))/3)+1)) #A matrix used for calculating the spaces available.
  possible_mothers <- matrix(0,(as.numeric(nrow(space)+1)),((as.numeric(nrow(subclones))-1)*2)) #A matrix used for saving the possible motherclones.
  rowsofhundred <- 40
  hundredpercentclones <- matrix(0,rowsofhundred,length(samples_unique)) #Matrix that is to be used in the cases where we have > 2 clones in a sample that have 100 %.
  hundredpercentclones[1,] <- samples_unique
  hpc <- 2
  cl <- 1
  nr_eq <- 6
  
  equalclones <- matrix(0,rowsofhundred,(length(samples_unique)*nr_eq)) #Matrix that is to be used in the cases where we have > 2 clones in a sample that have 100 %.
  
  equalclones[1,] <- rep(samples_unique,nr_eq)
  ec <- 2 #Count.
  ecl <- 1 #Column number.
  
  k <- 1
  for(k in 1:(ncol(sample_clone_matrix)/3)){ #Constructing a matrix were every two columns represent a sample. The first one tells us which subclone harbors the space and the second the remaining space on top of this sample.
    totalspace[1,(2*k-1)] <- sample_clone_matrix[1,(3*k-2)]
    k <- k+1
  }
  
  k <- 1
  for(k in 2:as.numeric(nrow(subclones))){ #Constructing a matrix were every two columns represent a subclone within a sample. The first one tells us which the chosen motherclone is and the other which other possible solutions there are.
    possible_mothers[1,(2*k-3)] <- subclones[k,1]
    k <- k+1
  }
  
  #Mother-daughter clone matrix.
  allocation_samples <- matrix(0,(as.numeric(nrow(clone_matrix_names))+1),(as.numeric(nrow(samples))+1)) #Matrix which is to be used for comparing the subclonal designation within each sample.
  allocation_samples[2:(as.numeric(nrow(clone_matrix_names))+1),1] <- clone_matrix_names[,3] #The subclone names are in the first row.
  allocation_samples[1,2:(as.numeric(nrow(samples))+1)] <- samples #The sample names are in the first column.
  
  
  subcloneswithinsample <- matrix(0,(as.numeric(nrow(sample_clone_matrix))-1),2)
  
  i <- 1
  for(i in 1:(as.numeric(ncol(sample_clone_matrix))/3)){ #Looping through all of the samples.
    #for(i in 1:5){
    
    subcloneswithinsample <- sample_clone_matrix[2:as.numeric(nrow(sample_clone_matrix)),(3*i-2):(3*i-1)] #Extraxting the subclonal architecture and TC for a certain sample.
    subcloneswithinsample_order <- subcloneswithinsample[order(as.numeric(subcloneswithinsample[,2]),decreasing = TRUE),] #Ordering the subclones from highest to lowest TC.
    
    #Ordering the subclones.
    ord <- 2
    subcloneswithinsample_order_old <- matrix(0,(as.numeric(nrow(sample_clone_matrix))-1),2)
    subcloneswithinsample_order_new <- subcloneswithinsample_order
    
    while(all(subcloneswithinsample_order_new == subcloneswithinsample_order_old) == FALSE){
      subcloneswithinsample_order_old <- subcloneswithinsample_order_new
      for(ord in 2:(as.numeric(nrow(subcloneswithinsample_order_old))-1)){ #Writing a function/loop that orders the subclones of the same size according to their median size.
        if(subcloneswithinsample_order_old[ord,2] != "0"){
          if(subcloneswithinsample_order_old[ord,2] == subcloneswithinsample_order_old[ord+1,2]){
            
            orderpos1 <- match(word(subcloneswithinsample_order_old[ord,1],2),namevector)
            orderpos2 <- match(word(subcloneswithinsample_order_old[ord+1,1],2),namevector)
            
            if(as.numeric(orderpos2) < as.numeric(orderpos1)){
              subcloneswithinsample_order_new <- subcloneswithinsample_order_old[c(1:(ord-1),ord+1,ord,(ord+2):nrow(subcloneswithinsample_order_old)), ]
            }
            
          }
        }
        ord <- ord+1
      }
    }
    
    subcloneswithinsample_order <- subcloneswithinsample_order_new
    #print(subcloneswithinsample_order)
    
    j <- 1
    ecl_original <- ecl
    equal <- 1
    for(j in 2:as.numeric(nrow(sample_clone_matrix))){ #Looping through the subclones within the sample.
      #for(j in 2:4){
      
      if(j == 2){ #We're in the first position. This is the ALL-event.
        space[1,1] <- subcloneswithinsample_order[j-1,1] #The name.
        space[1,2] <- subcloneswithinsample_order[j-1,2] #The TC.
      }
      
      
      if(j != 2){
        if(subcloneswithinsample_order[j-1,1] != "0"){
          if(subcloneswithinsample_order[j-1,1] != "ALL"){ #We should not add it again.
            maxspace <- which.max(space[,2]) #Finding the largest available space.
            newname <- subcloneswithinsample_order[j-1,1] #The name of the new subclone.
            newspace <- subcloneswithinsample_order[j-1,2] #The space of the new subclone.
            maxname <- space[maxspace,1]
            #Adding a loop that checks whether or not this is the only possible solution for the subclone to be placed as a daughter to.
            c <- 1
            so <- 0
            for(c in 1:nrow(space)){
              if(as.numeric(space[c,2]) != 0){ #The new test space should not be zero. You cannot put anything there.
                if((as.numeric(space[c,2])-as.numeric(newspace)) >= -0.1){ #Added this due to the simulation. Too many decimals. Rounding makes events not being placed in parallel.
                  
                  daughter_pos <- match(paste(sample_clone_matrix[1,(3*i-2)],newname), possible_mothers[1,])
                  
                  if(c == maxspace){
                    possible_mothers[2,daughter_pos] <- space[c,1] #Adding the original solution.
                  }
                  if(c != maxspace){
                    #print("There are other solutions")
                    if(space[c,1] %in% possible_mothers[2,(as.numeric(daughter_pos))] == FALSE){
                      #print("Now we will add it")
                      daughter_pos <- match(paste(sample_clone_matrix[1,(3*i-2)],newname), possible_mothers[1,])
                      possible_mothers[(2+so),(as.numeric(daughter_pos)+1)] <- space[c,1]
                      
                      #Tystade detta 200720 samt 200820. Made some rules disappear.
                      # if(space[c,2] == newspace){ #Beh?ver ju dock inte inneb?ra att de faktiskt ?r equalclones och kan placeras i varandra.
                      #   if(space[c,1] != "ALL"){ #Varf?r inte ALL? Man kan f? fel.
                      #     #print(space)
                      #     #print(c)
                      #
                      #     mothername <- paste(sample_clone_matrix[1,(3*i-2)],space[c,1])
                      #     mothercolumn <- match(mothername,possible_mothers[1,])
                      #     possible_mothers[1,mothercolumn]
                      #     #print("h?r ?r mothername and newname")
                      #     #print(mothername)
                      #     #print(newname)
                      #
                      #     n <- 1
                      #     for(n in 2:nrow(possible_mothers)){
                      #       #print("H?r ?r ett n")
                      #       #print(n)
                      #       #print(mothercolumn)
                      #       #print(possible_mothers[n,(mothercolumn+1)])
                      #       if(possible_mothers[n,(mothercolumn+1)] == "0"){
                      #         possible_mothers[n,(mothercolumn+1)] <- newname
                      #         break
                      #       }
                      #       n <- n+1
                      #     }
                      #   }
                      # }
                      
                      so <- so+1
                    }}
                }
              }
              c <- c+1
            }
            
            space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
            space[s,1] <- newname #Adding the new spacename and space size to the spacematrix.
            space[s,2] <- newspace
            
            #Treating the case when space[maxspace,2] = 100 % and the newspace as well. Then the motherclone and the daughterclone are both part of the base.
            if(subcloneswithinsample_order[j-1,2] == "100"){
              if(subcloneswithinsample_order[j-1,1] != "ALL"){
                if(subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),2] == "100"){
                  if(subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1] != "ALL"){
                    
                    newspace_name <- subcloneswithinsample_order[j-1,1]
                    maxspace_name <- subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1]
                    
                    if(hpc == 2){
                      hundredpercentclones[2,cl] <- subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1]
                      hundredpercentclones[3,cl] <- subcloneswithinsample_order[j-1,1]
                      hpc <- 4
                    }else{
                      if(subcloneswithinsample_order[j-1,1] %in% hundredpercentclones[,cl] == FALSE){
                        hundredpercentclones[hpc,cl] <- subcloneswithinsample_order[j-1,1]
                        hpc <- hpc + 1
                      }
                    }
                    
                    allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- maxspace_name #Annoting the mother clone of each of the subclones within each sample.
                    allocation_samples[match(maxspace_name,allocation_samples[,1]),(i+1)] <- newspace_name #Annoting the mother clone of each of the subclones within each sample.
                    
                  }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
                  }
                }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
                }
              }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
              }
            }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
            }
            
            
            #Treating the case when we have multiple clones of equal size that have to be placed inside each other.
            if(subcloneswithinsample_order[j-1,1] != "ALL"){
              if(subcloneswithinsample_order[j-1,2] != "100"){
                if(as.numeric(subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),2]) == as.numeric(subcloneswithinsample_order[j-1,2])){ #It should be equal in size to the other cluster.
                  
                  newspace_name <- subcloneswithinsample_order[j-1,1]
                  maxspace_name <- subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1]
                  
                  if(ec == 2){ #We have not yet added any events to the equalcolumn for this sample.
                    equalclones[2,ecl] <- subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1]
                    equalclones[3,ecl] <- subcloneswithinsample_order[j-1,1]
                    ec <- 4
                  }else{ #We have added events earlier. We now want to add even one more to this column.
                    if(subcloneswithinsample_order[j-1,1] %in% equalclones[,ecl] == FALSE){
                      equalclones[ec,ecl] <- subcloneswithinsample_order[j-1,1]
                      ec <- ec + 1
                    }
                  }
                  
                  #This part adds the names such that they get each others names in the allocation_samples matrix.
                  #Silenced this one 200308 since it sometimes made events be allocated in a weird way. 100 % events got each others and 50 % each others but the 50 did not get the 100.
                  #allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- maxspace_name #Annoting the mother clone of each of the subclones within each sample.
                  #allocation_samples[match(maxspace_name,allocation_samples[,1]),(i+1)] <- newspace_name #Annoting the mother clone of each of the subclones within each sample.
                  equal <- 2
                  
                }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
                if(equal == 2){
                  ecl <- ecl+length(samples_unique)
                  ec <- 2
                  equal <- 1}
                }
              }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
              }
            }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
            }
            
          } #!= "ALL".
        } #!= "0".
        
        if(j == as.numeric(nrow(sample_clone_matrix))){ #We're at the end of a sample.
          totalspace[2:as.numeric(nrow(totalspace)),((t*2)-1):((t*2))] <- space
          t <- t+1
          
          s <- 2 #Resetting s and space.
          space <- matrix(0,50,2)
          
        }else{s <- s+1}
        
      } #j != 2.
      
      
      j <- j+1
    } #Subclones within a sample.
    
    hpc <- 2
    cl <- cl + 1
    ec <- 2
    ecl <- ecl_original + 1
    i <- i+1
    
  } #Samples.
  
  
  i <- 1
  Clustering <- file_samples_subclones
  Clustering[,12] <- "No"
  Clustering <- Clustering[,Clustering[1,]!="No"]
  colnames(Clustering)[12] <- "Cluster"
  for(i in 1:nrow(Clustering)){
    
    w1 <- word(Clustering[i,12],1)
    w3 <- word(Clustering[i,12],3)
    w2 <- "Cluster_"
    
    if(w1 != "ALL"){
      Clustering[i,12] <- paste(w1,w2,w3,sep=" ")}
    
    i <- i+1
  }
  
  assign("Clustering", Clustering, envir=globalenv()) #This is a matrix illustrating the events and their subclonal belonging.
  assign("file_samples_subclones", file_samples_subclones, envir=globalenv()) #This is a matrix illustrating the events and their subclonal belonging.
  assign("possible_mothers", possible_mothers, envir=globalenv()) #This is a matrix illustrating the chosen mother clone as well as other possible mothers.
  assign("allocation_samples", allocation_samples, envir=globalenv()) #The mother-daughter division are exported to the global environment.
  assign("equalclones", equalclones, envir=globalenv()) #The equal clones.
  assign("hundredpercentclones", hundredpercentclones, envir=globalenv()) #The equal clones.
  
  #Fusing the equalclones and the hundredpercentclones.
  i <- 1
  
  for(i in 1:ncol(hundredpercentclones)){
    
    if(hundredpercentclones[2,i] != "0"){ #We have some hundredpercentlones in this sample.
      
      if(equalclones[2,i] == "0"){ #We do not have any other equalclones in this sample.
        equalclones[2:nrow(equalclones),i] <- hundredpercentclones[2:nrow(hundredpercentclones),i] #We paste it here.
        
      }else if(equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){ #We have something but not in the next one.
        equalclones[2:nrow(equalclones),(i+as.numeric(ncol(hundredpercentclones)))] <- hundredpercentclones[2:nrow(hundredpercentclones),i] #We paste it here.
      }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){ #We have something but not in the nextnext one.
        equalclones[2:nrow(equalclones),(i+as.numeric(ncol(hundredpercentclones)))] <- hundredpercentclones[2:nrow(hundredpercentclones),i] #We paste it here.
      }
      
    }
    
    i <- i+1
  }
  
  #The equalclones should have each other's mothers. Created 200804.
  i <- 1
  save_order <- 1
  unique_mothers <- list()
  for(i in 1:ncol(equalclones)){ #Looping through the clones of equal size.
    removed <- 0
    if(equalclones[2,i] != "0"){
      
      equal_mothers <- as.vector(equalclones[equalclones[,i] != "0",i])
      order <- matrix(0,length(equal_mothers),2)
      overview_sub <- matrix(0,length(equal_mothers),ncol(overview))
      order[,1] <- as.matrix(equal_mothers)
      
      j <- 2
      for(j in 2:length(equal_mothers)){
        
        daughter <- paste(equal_mothers[1],equal_mothers[j])
        daughter_column <- match(daughter,possible_mothers[1,])
        all <- rbind(as.matrix(possible_mothers[,daughter_column]),as.matrix(possible_mothers[,daughter_column+1]))
        all_nozero <- as.matrix(all[all!="0",])
        
        #if(length(all_nozero)>1){ #Added this for the cases where there is none.
        unique_mothers[[(j-1)]] <- all_nozero[2:nrow(all_nozero),] #Saving all unique mothers. 2 because we do not want the daughter.
        #}
        
        file_samples_subclones_row <- match(daughter,file_samples_subclones[,13])
        overview_row <- match(file_samples_subclones[file_samples_subclones_row,10],overview[,1])
        overview_sub[j,] <- overview[overview_row,]
        j <- j+1
      }
      
      unique_mother_tot <- as.matrix(unique(unlist(unique_mothers)))
      
      j <- 2
      for(j in 3:nrow(overview_sub)){
        overview_sub <- overview_sub[,overview_sub[j,]!="0"] #Removing columns where not all of them are present simultaneously.
        j <- j+1
      }
      
      overview_sub <- as.matrix(overview_sub)
      overview_sub_n <- overview_sub[2:nrow(overview_sub),2:as.numeric(ncol(overview_sub))]
      class(overview_sub_n) <- "numeric"
      
      #Checking if the equalclones really are allowed to be placed inside each other in all samples.
      m <- 1
      again <- 0
      
      if(is.null(nrow(overview_sub_n))==FALSE){
        for(m in 1:nrow(overview_sub_n)){ #Looping through the clones in this equalclones.
          
          n <- 1
          for(n in 1:nrow(overview_sub_n)){
            
            thesigns <- unique(sign(overview_sub_n[m,]-overview_sub_n[n,])) #Subtracting the rows.
            
            if("-1" %in% thesigns && "1" %in% thesigns && again == 0 && length(equalclones[equalclones[,i]%in% hundredpercentclones,i])!=nrow(equalclones)){ #the last part was added due to handle cases where we have contradictions in the data so that we do not remove it.
              # print("These should not be in equalclones.")
              again <- 1 #We should not remove it multiple times.
              removed <- 1
              eq_segment <- ceiling(as.numeric(i)/as.numeric(length(samples_unique))) #Telling us which segment of equalclones we are in.
              
              segments_left <- as.numeric(nr_eq)-as.numeric(eq_segment)
              therest <- as.numeric(i)%%length(samples_unique) #Calculating in which biopsy we are.
              if(as.numeric(therest) == 0){
                biopsy <- nr_eq
              }else{
                biopsy <- therest
              }
              if(as.numeric(segments_left) != 0){ #There is at least one that should be moved.
                
                p <- as.numeric(eq_segment)+1
                
                for(p in (as.numeric(eq_segment)+1):as.numeric(nr_eq)){
                  if(as.numeric(biopsy)+(as.numeric(length(samples_unique))*(as.numeric(p)-1)) <= as.numeric(ncol(equalclones))){
                    equalclones[,as.numeric(biopsy)+(as.numeric(length(samples_unique))*(as.numeric(p)-2))] <- equalclones[,as.numeric(biopsy)+(as.numeric(length(samples_unique))*(as.numeric(p)-1))]
                  }
                  p <- p+1
                }
                
              }else{
                equalclones[2:nrow(equalclones),i] <- "0"
              }
            }
            n <- n+1
          }
          m <- m+1
        }
      }
      
      if(is.null(nrow(overview_sub_n))==FALSE){
        if(removed == 0){
          
          order[2:nrow(order),2] <- as.matrix(rowSums(overview_sub_n))
          assign("order",order,envir=globalenv())
          order[2:nrow(order),2] <- as.numeric(order[2:nrow(order),2])
          order_new <- order[order(as.numeric(as.matrix(order[,2])),decreasing=TRUE),] #Ordering the matrix after size.
          
          order[1,] <- order_new[nrow(order_new),]
          order[2:nrow(order),] <- order_new[1:(nrow(order_new)-1),]
          
          # if(save_order == 1){ #Silenced it 210316 since equalclones_order and order had differing numbers of rows.
          #   equalclones_order <- order
          #   save_order <- 2
          # }else{
          #   print(equalclones_order)
          #   print(order)
          #   print(nrow(equalclones_order))
          #   print(nrow(order))
          #   equalclones_order <- cbind(equalclones_order,order)
          # }
          
          j <- 2
          for(j in 2:nrow(order)){
            
            daughter_column <- match(paste(order[1,1],order[j,1]),possible_mothers[1,])
            
            if(j == 2){
              
              events <- as.matrix(unique_mother_tot[unique_mother_tot %in% order == FALSE,]) #Finding which possible mothers are not one of the equal clones. Why?
              #events <- as.matrix(unique_mother_tot[unique_mother_tot != order[j,1],]) #Finding which possible mothers are left.
              
              if(nrow(events)==1){
                possible_mothers[2,daughter_column] <- events
              }else{
                #print(possible_mothers[1,daughter_column])
                possible_mothers[2,daughter_column] <- events[1,]
                possible_mothers[2:nrow(events),(daughter_column+1)] <- events[2:nrow(events),]
              }
              
            }else{
              mother_column <- match(paste(order[1],order[j-1,1]),possible_mothers[1,])
              #print(paste(order[1],order[j-1,1]))
              #print(mother_column)
              possible_mothers[2,daughter_column] <- word(possible_mothers[1,mother_column],2,3)
              
              events <- as.matrix(unique_mother_tot[unique_mother_tot != order[j,1],]) #Finding which possible mothers are not one of the equal clones.
              as.matrix(events[word(possible_mothers[1,mother_column],2,3) != order[j,1],])
              
              possible_mothers[2:nrow(events),(daughter_column+1)] <- "0" #events[2:nrow(events),1] #"0" was added 200821.
              
            }
            j <- j+1
          }
        }#Removed.
      }
    }
    i <- i+1
  }
  
  assign("possible_mothers_2",possible_mothers,envir=globalenv())
  assign("equalclones_new", equalclones, envir=globalenv()) #The equal clones.
  #assign("equalclones_order",equalclones_order,envir=globalenv())
  #########################################################################################
  #We want to add all equalclones to possible mothers if one of them are a possible mother#
  #########################################################################################
  # i <- 1
  # for(i in 1:(ncol(possible_mothers)/2)){
  #   
  #   daughter <- possible_mothers[1,(2*i-1)]
  #   mother <- possible_mothers[2,(2*i-1)] #Right now I only check the "primary mother".
  #   eq_col <- which(word(daughter,1)==equalclones[1,])
  #   eq_sample <- equalclones[,eq_col] #All equalclones in this sample.
  # 
  #   if(mother %in% eq_sample){ #Checking if the mother is in equalclones.
  # 
  #     j <- 1
  #     for(j in 1:length(eq_col)){
  #       
  #       if(mother %in% eq_sample[,j]){
  #        
  #        eq_sub <- eq_sample[eq_sample[,j]!=mother,j] 
  #        eq_sub <- eq_sub[eq_sub!="0"]
  # 
  #         othermothers <- possible_mothers[possible_mothers[,(2*i)]!="0",(2*i)]
  #         pos_m_rows <- length(othermothers)
  #         eq_sub <- eq_sub[eq_sub!=othermothers]
  #         
  #         if(length(eq_sub)!=0 && is.na(length(pos_m_rows))==FALSE){
  #         possible_mothers[(pos_m_rows+2):(pos_m_rows+length(eq_sub)),(2*i)] <- eq_sub[2:length(eq_sub)]}
  #        
  #       }
  #       
  #       j <- j+1
  #     }
  #   }
  #   i <- i+1
  # }
  #   
  ###################################################################################
  #Looking for discrepancies between the subclonal architecture of different samples#
  ###################################################################################
  #Taking rules into consideration.
  if(missing(rule)==FALSE){
    
    i <- 1
    for(i in 1:nrow(rule)){
      row_m <- match(rule[i,1],file_samples_subclones[,10]) 
      row_d <- match(rule[i,2],file_samples_subclones[,10])
      
      rule[i,1] <- file_samples_subclones[row_m,12]
      rule[i,2] <- file_samples_subclones[row_d,12]
      rule <- as.matrix(rule)
    }
    
  }else{
    rule <- matrix(0,1,2)
  }
  
  assign("rule_new",rule,envir=globalenv())
  
  theonlymothers <- matrix(0,as.numeric(nrow(possible_mothers)),as.numeric(ncol(possible_mothers)))
  if(as.numeric(length(unique(file_samples_subclones[,2]))) > 2){ #If we only have one biopsy we do not have to compare stuff.
    i <- 2
    x <- matrix(0,2,as.numeric(ncol(allocation_samples)))
    x[1,] <- allocation_samples[1,]
    tom <- 1 #We will start to save down the data in column 1.
    Event_rule_removed <- 0
    not_again <- 0
    removed <- 0
    stop_while <- 0
    again <- 1
    while(i <= as.numeric(nrow(allocation_samples))){ #Looping through the subclones.
      #while(i <= 6){
      # print("Here is i - Subclones")
      # print(i)
      only <- 0
      x[2,] <- allocation_samples[i,1:as.numeric(ncol(allocation_samples))] #Extracting information about the motherclones in all samples.
      y <- as.data.frame(allocation_samples[i,2:as.numeric(ncol(allocation_samples))])
      
      if(length(allocation_samples[i,allocation_samples[i,2:as.numeric(ncol(allocation_samples))]!="0"]) > 1 || x[2,1] %in% equalclones){ #If we only have this event in one sample we want to solely go on largest space.
        y <- y[y[]!=0] #Removing the data points containing zeros. We now have all the motherclones.
        
        mother_all_biopsies <- matrix(0,as.numeric(nrow(possible_mothers))+1,ncol(allocation_samples))
        mother_all_biopsies[1,] <- allocation_samples[1,]
        
        k <- 2
        for(k in 2:ncol(x)){ #Looping through the daughter clones.
          daughtersubclone <- paste(x[1,k],x[2,1]) #Finding the name of the daughter subclone.
          daughterposition <- match(daughtersubclone,possible_mothers[1,]) #Finding the position for the daughter subclone in the possible_mothers matrix.
          mother_all_biopsies[2,k] <- possible_mothers[2,as.numeric(daughterposition)]
          mother_all_biopsies[3:nrow(mother_all_biopsies),k] <- possible_mothers[2:nrow(possible_mothers),(as.numeric(daughterposition)+1)] #This matrix will contain all mothers in all biopsies.
          k <- k+1
        }
        
        distribution <- table(mother_all_biopsies[2:nrow(mother_all_biopsies),3:ncol(mother_all_biopsies)])
        mother_not_all <- distribution[distribution!=length(y)]
        mother_all <- distribution[distribution==length(y)] #This table illustrates the mothers that can be given in all biopsies.
        
        if(length(mother_all) > 1){
          mother_all <- mother_all[mother_all>=length(y)]
        }
        
        mother_almost_all <- distribution[distribution==(length(y)-1)]
        mother_not_all <- mother_not_all[rownames(mother_not_all)!="0"] #This table illustrates all mothers that cannot be given in each biopsy.
        
        j <- 2
        mp <- 0
        count_replace <- 1 #Used when we get a rule where we only have one mother in one sample and have to change earlier clones.
        for(j in 2:ncol(x)){ #Looping through the motherclones.
          if(x[2,j] != "0"){ #We do not want to analyze situations where we do not even have the daughter subclone in question.
            daughtersubclone <- paste(x[1,j],x[2,1]) #Finding the name of the daughter subclone.
            daughterposition <- match(daughtersubclone,possible_mothers[1,]) #Finding the position for the daughter subclones in the possible_mothers matrix.
            justzeros <- table(possible_mothers[,(as.numeric(daughterposition)+1)]) #Extracting the column for the other possible mothers.
            
            #print("Before")
            if((as.numeric(justzeros[1])/(as.numeric(nrow(possible_mothers)))) == 1){ #If it is one, then every position is a zero. This means that this is the only solution for this mother-daughter allocation in this sample.
              daughterrowdata <- match(daughtersubclone,file_samples_subclones[,13]) #Finding the row of the subclone in the file. This will be used in order to obatin the clone size.
              only <- 1
              onlymother <- possible_mothers[2,as.numeric(daughterposition)]
              
              #print("Now we are changing things backwards!")
              column_equal <- which(x[1,j]==equalclones[1,])
              equalclones_biopsy <- equalclones[,column_equal]
              
              if(x[1,j] != "B1" && onlymother %in% mother_all %in% equalclones == FALSE){ #This does not matter if we only have one biopsy. why is it not allowed to be in equalclones?
                if(count_replace!=1 && tom > 1){
                  sub <- which(word(theonlymothers[1,1:(tom-1)],2,3)==word(theonlymothers[1,tom],2,3))
                  howmany <- length(sub)
                  if(length(sub)==0){howmany <- 0}
                  if(howmany!=0){
                    if(theonlymothers[3,tom-1] == "0"){
                      theonlymothers[3,(tom-as.numeric(howmany)):(tom-1)] <- theonlymothers[2,(tom-as.numeric(howmany)):(tom-1)] #Changed from j to count_replace. Otherwise we get problems when certain events are not present in all biopsies.
                    }
                    theonlymothers[2,(tom-as.numeric(howmany)):(tom-1)] <- onlymother
                  }
                }
              }
              
              if(possible_mothers[2,as.numeric(daughterposition)] != "ALL"){ #Extracting the name of the mother it has to have.
                mothername <- paste(word(possible_mothers[1,as.numeric(daughterposition)],1),possible_mothers[2,as.numeric(daughterposition)]) #The name of the only possible mother.
                motherrowdata <- match(mothername,file_samples_subclones[,13]) #Finding its row in the file in order to obtain the clone size.
              }else{mothername <- "ALL"
              motherrowdata <- 1}
              
              if(mothername != "ALL" && as.numeric(file_samples_subclones[daughterrowdata,11]) + as.numeric(file_samples_subclones[motherrowdata,11]) != 200){ #If they were they they should be each other's motherclones.
                theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                theonlymothers[2,tom] <- possible_mothers[2,as.numeric(daughterposition)] #The mother clone which it has to have.
                tom <- tom + 1
                norules <- 0
              }else if(mothername == "ALL"){ #The only possible mother for this clone are the "ALL" events.
                theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                theonlymothers[2,tom] <- possible_mothers[2,as.numeric(daughterposition)] #The mother clone which it has to have.
                tom <- tom + 1
                norules <- 0
              }
            }else{ #There are multiple solutions.
              #print("Multiple solutions")
              
              if(is.na(match(x[2,1],hundredpercentclones)) == FALSE){ #The subclone is present in a hundredpercentclones.
                columnhundred <- round(match(x[2,1],hundredpercentclones)/nrow(hundredpercentclones))+1 #The column. This column contains all of the alterations in the hundredpercentclone.
                rowhundred <- match(x[2,1],hundredpercentclones)- columnhundred*nrow(hundredpercentclones)
                if(is.na(match(x[2,j],hundredpercentclones[,columnhundred])) == FALSE){ #The mother exist in the same hundredpercentclone.
                  theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                  theonlymothers[2,tom] <- possible_mothers[2,as.numeric(daughterposition)] #The mother clone which it has to have.
                  tom <- tom + 1
                }else{
                  k <- 2
                  for(k in 2:nrow(possible_mothers)){ #Looping through the other solutions.
                    if(possible_mothers[k,(as.numeric(daughterposition)+1)]!= "0"){
                      if(is.na(match(possible_mothers[k,(as.numeric(daughterposition)+1)],hundredpercentclones[,columnhundred])) == FALSE){
                        theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                        theonlymothers[2,tom] <- possible_mothers[k,as.numeric(daughterposition)+1] #The mother clone which it has to have.
                        tom <- tom + 1
                      }
                    }
                    k <- k+1
                  }
                }
              }else if(length(names(mother_all)) < 1){
                #print("There is no mother that is possible in all samples.")
                if(only == 1){
                  #print("Onlymother")
                  theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                  theonlymothers[2,tom] <- onlymother
                  tom <- tom+1
                }else if(length(mother_almost_all) != 0){ #There is a possibility that this almost event is the true one.
                  
                  p <- 1
                  for(p in 1:length(mother_almost_all)){ #Looping through the mothers that are possible in almost all samples.
                    
                    columns <- which(theonlymothers[2,]==names(mother_almost_all)[p]) #Finding all places where this mother is present.
                    
                    biopsy_nr <- x[1,j] #The biopsy we are looking at.
                    sample_columns <- which(word(theonlymothers[1,],1)== biopsy_nr) #Finding all positions belonging to this biopsy.
                    match_columns <- intersect(columns,sample_columns) #All the rules for events being placed in this mother.
                    
                    #Information about the mother in this sample.
                    if(names(mother_almost_all)[p]=="ALL"){
                      mother_rule <- "ALL"
                    }else{
                      mother_rule <- paste(biopsy_nr,names(mother_almost_all)[p]) #Name
                    }
                    row_TC_mother_rule <- match(mother_rule,file_samples_subclones[,13]) #Position.
                    mother_rule_size <- file_samples_subclones[row_TC_mother_rule,11] #Size in the sample.
                    
                    if(length(match_columns) >= 2){ #Changed to 2.
                      #print("There is a rule for this mother.")
                      
                      our_new_daughter <- daughtersubclone
                      row_TC_our_new_daughter_rule <- match(our_new_daughter,file_samples_subclones[,13])
                      our_new_daughter_rule_size <- file_samples_subclones[row_TC_our_new_daughter_rule,11]
                      
                      #Calculating if there is any room left.
                      #print("Calculating if there is any room left")
                      r <- 1
                      for(r in 1:length(match_columns)){
                        daughter_rule <- theonlymothers[1,as.numeric(match_columns[r])]
                        row_TC_daughter_rule <- match(daughter_rule,file_samples_subclones[,13])
                        daughter_rule_size <- file_samples_subclones[row_TC_daughter_rule,11]
                        
                        if((as.numeric(daughter_rule_size)-as.numeric(our_new_daughter_rule_size)) > 0){ #If our new daughter is larger than the ones we're comparing with now, it is not interesting to subtract them since this new alteration will have f?retr?de.
                          mother_rule_size <- as.numeric(mother_rule_size)-as.numeric(daughter_rule_size)}
                        r <- r+1
                      }
                      
                      if(is.na(as.numeric(mother_rule_size))==FALSE){
                        if(is.na(as.numeric(our_new_daughter_rule_size))==FALSE){
                          if((as.numeric(mother_rule_size)+0.1) < as.numeric(our_new_daughter_rule_size)){ #Added 0.1 because otherwise you might get rounding errors.
                            #print("There is no longer room.")
                            
                            pos_rem1 <- match(names(mother_almost_all)[p],possible_mothers[,daughterposition])
                            
                            if(is.na(pos_rem1) == FALSE){
                              possible_mothers[pos_rem1,daughterposition] <- "0"
                            }
                            
                            if(daughterposition < ncol(possible_mothers)){
                              pos_rem2 <- match(names(mother_almost_all)[p],possible_mothers[,(daughterposition+1)])
                              if(is.na(pos_rem2) == FALSE){
                                possible_mothers[pos_rem2,(daughterposition+1)] <- "0"
                              }
                            }
                            
                            removed <- 1
                            Event_rule_removed <- 1 #Indicating that an event has been removed.
                          }}}
                    }
                    p < p+1
                  }
                  
                  theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                  theonlymothers[2,tom] <- names(mother_almost_all)[1] #The mother clone which is possible in almost all samples.
                  tom <- tom + 1
                }else{
                  theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)]
                  theonlymothers[2,tom] <- names(distribution)[2]
                  tom <- tom+1 #Added this since we do not get the tom-count during the second round.
                }
              }else if(length(mother_all) == 1){
                #Multiple solutions
                #Adding a rule in order to make all the daughters originate after the same mother.
                #print("The only solution now.")
                daughtersubclone <- paste(x[1,j],x[2,1]) #Finding the name of the daughter subclone.
                daughterposition <- match(daughtersubclone,possible_mothers[1,]) #Finding the position for the daughter subclones in the possible_mothers matrix.
                
                theonlymothers[2:nrow(theonlymothers),tom] <- "0"
                theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                theonlymothers[2,tom] <- names(mother_all)[1] #The mother clone which it has to have.
                
                if(length(mother_almost_all) != 0){
                  theonlymothers[3,tom] <- names(mother_almost_all)[1]
                }
                
                tom <- tom + 1
                norules <- 0
              }else if(length(mother_all) > 1){ #There are multiple solutions that are possible in all samples.
                #print("We have multiple possible allocations that are equally probable")
                daughtersubclone <- paste(x[1,j],x[2,1]) #Finding the name of the daughter subclone.
                daughterposition <- match(daughtersubclone,possible_mothers[1,]) #Finding the position for the daughter subclones in the possible_mothers matrix.
                theonlymothers[2:nrow(theonlymothers),tom] <- "0"
                
                #Adding a rule algorithm.
                rule_d <- which(rule[,2]==x[2,1])
                r <- 1
                rule_applied <- 0
                for(r in 1:length(mother_all)){
                  
                  rule_m <- which(rule[,1]==names(mother_all)[r])
                  rule_both <- intersect(rule_d,rule_m)
                  if(length(rule_both)>0){
                    #print("They are not allowed to be placed after one another.")
                    mother_all <- mother_all[names(mother_all)!=names(mother_all)[r]] #Removing this mother entirely.
                    rule_applied <- 1
                  }
                  
                  r <- r+1
                }
                
                #Algorithm for choosing between many mothers that are possible in all samples.
                if(ncol(overview) > 3){ #If it is 3 we do only have one sample. No reason to look for patterns.
                  mother_all_mtrx <- as.matrix(names(mother_all)) #The mother names.
                  if(length(names(mother_all))==1){
                    mother_all_type <- matrix(0,2,ncol(overview))
                  }else{
                    mother_all_type <- matrix(0,(nrow(mother_all)+1),ncol(overview))}
                  
                  n <- 1
                  
                  for(n in 1:(nrow(mother_all_mtrx)+1)){
                    if(n == 1){
                      mother_all_mtrx_row<- match(possible_mothers[1,as.numeric(daughterposition)],file_samples_subclones[,13]) #Position.
                      mother_all_type[n,1] <- file_samples_subclones[mother_all_mtrx_row,12]
                      
                      type <- word(file_samples_subclones[mother_all_mtrx_row,10],1:3)
                      mother_all_type[n,2] <- paste(type[1],type[2],type[3],sep=" ")
                      
                      overview_row <- match(file_samples_subclones[mother_all_mtrx_row,10], overview[,1])
                      mother_all_type[n,3:ncol(mother_all_type)] <- overview[overview_row,3:ncol(overview)]
                    }else{
                      
                      if(mother_all_mtrx[n-1,1] != "ALL"){
                        
                        mother_all_type[n,1] <- mother_all_mtrx[n-1,1]
                        
                        mother_all_mtrx_row <- match(mother_all_mtrx[n-1,1],word(file_samples_subclones[,13],2,3)) #Position.
                        type <- word(file_samples_subclones[mother_all_mtrx_row,10],1:3)
                        mother_all_type[n,2] <- paste(type[1],type[2],type[3],sep=" ")
                        
                        overview_row <- match(file_samples_subclones[mother_all_mtrx_row,10], overview[,1])
                        mother_all_type[n,3:ncol(mother_all_type)] <- overview[overview_row,3:ncol(overview)]
                        
                      }else{
                        mother_all_type[n,1] <- "ALL"
                        mother_all_type[n,2] <- "ALL"
                        mother_all_type[n,3:ncol(mother_all_type)] <- "100"
                      }
                    }
                    
                    n <- n+1
                  }
                  m <- 4
                  for(m in 4:ncol(mother_all_type)){
                    if(is.null(nrow(mother_all_type))==FALSE){
                      sign_vector <- (as.numeric(mother_all_type[,(as.numeric(m)-1)]) - as.numeric(mother_all_type[,m]))
                      same_sign <- as.matrix(sign(sign_vector[1]) == sign(sign_vector))
                      mother_all_type <- mother_all_type[same_sign,]
                    }
                    
                    m <- m+1
                  }
                  
                  if(is.null(nrow(mother_all_type))==FALSE){
                    mother_all_type <- mother_all_type[!is.na(mother_all_type[,1]),]
                    if(is.null(nrow(mother_all_type))==FALSE){
                      if(nrow(mother_all_type) > 1){
                        #print("There are mothers that follow the same pattern.")
                        mother_preferred <- mother_all_type[2,1] #We prefer to choose a mother that shows a similar pattern.
                        mp <- 1
                      }}
                  }
                  #print("Final")
                }
                
                #Adding the events.
                if(not_again != 1){
                  if(mp != 1){
                    theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                    if(possible_mothers[2,as.numeric(daughterposition)]%in% names(mother_all)){ #If the largest space mother is in this, we will choose it.
                      theonlymothers[2,tom] <- possible_mothers[2,as.numeric(daughterposition)]
                    }else{
                      theonlymothers[2,tom] <- names(mother_all)[1] #The mother clone.
                      theonlymothers[3:(length(mother_all)+1),tom] <- names(mother_all)[2:length(mother_all)]}
                  }else{
                    theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                    theonlymothers[2,tom] <- mother_preferred
                    pos_mp <- match(mother_preferred,names(mother_all))
                    mother_all[pos_mp] <- mother_all[1] #Replacing the first event with the preferred one.
                    theonlymothers[3:(length(mother_all)+1),tom] <- names(mother_all)[2:length(mother_all)]
                  }
                }else{
                  if(mp != 1){
                    
                    if(possible_mothers[2,as.numeric(daughterposition)]%in% names(mother_all)){ #If the largest space mother is in this, we will choose it.
                      theonlymothers[2,tom] <- possible_mothers[2,as.numeric(daughterposition)]
                    }else{
                      theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                      theonlymothers[2,tom] <- names(mother_all)[1] #The mother clone which it has to have.
                      theonlymothers[3,tom] <- names(mother_all)[2]
                      if(length(mother_all) >= 3){
                        theonlymothers[4:(length(mother_all)+1),tom] <- names(mother_all)[3:length(mother_all)]}
                    }
                  }else{
                    theonlymothers[1,tom] <- possible_mothers[1,as.numeric(daughterposition)] #The clone which is to be allocated.
                    pos_mp <- match(mother_preferred,names(mother_all))
                    theonlymothers[2,tom] <- mother_preferred
                    mother_all[pos_mp] <- mother_all[1] #Replacing the first event with the preferred one.
                    theonlymothers[3:(length(mother_all)+1),tom] <- names(mother_all)[2:length(mother_all)]
                  }
                }
                
                #We have to see if any of these mothers are not possible any more since we've gotten some rules for the allocations.
                if(not_again != 1){
                  p <- 1
                  for(p in 2:(length(mother_all)+1)){ #Looping through the mothers that are possible in all samples.
                    
                    columns <- which(theonlymothers[2,]==theonlymothers[p,tom]) #Finding all places where this mother is present.
                    
                    biopsy_nr <- x[1,j] #The biopsy we are looking at.
                    sample_columns <- which(word(theonlymothers[1,],1)== biopsy_nr) #Finding all positions belonging to this biopsy.
                    match_columns <- intersect(columns,sample_columns) #All the rules for events being placed in this mother.
                    
                    
                    #Information about the mother in this sample.
                    if(theonlymothers[p,tom]=="ALL"){
                      mother_rule <- "ALL"
                    }else{
                      mother_rule <- paste(biopsy_nr,theonlymothers[p,tom]) #Name
                    }
                    row_TC_mother_rule <- match(mother_rule,file_samples_subclones[,13]) #Position.
                    mother_rule_size <- file_samples_subclones[row_TC_mother_rule,11] #Size in the sample.
                    
                    if(length(match_columns) > 1){
                      #print("There is a rule for this mother.")
                      our_new_daughter <- daughtersubclone
                      row_TC_our_new_daughter_rule <- match(our_new_daughter,file_samples_subclones[,13])
                      our_new_daughter_rule_size <- file_samples_subclones[row_TC_our_new_daughter_rule,11]
                      
                      #Calculating if there is any room left.
                      r <- 1
                      #print("Calculating if there is any room left")
                      for(r in 1:length(match_columns)){
                        daughter_rule <- theonlymothers[1,as.numeric(match_columns[r])]
                        row_TC_daughter_rule <- match(daughter_rule,file_samples_subclones[,13])
                        daughter_rule_size <- file_samples_subclones[row_TC_daughter_rule,11]
                        
                        if((as.numeric(daughter_rule_size)-as.numeric(our_new_daughter_rule_size)) > 0){ #If our new daughter is larger than the ones we're comparing with now, it is not interesting to subtract them since this new alteration will have f?retr?de.
                          mother_rule_size <- as.numeric(mother_rule_size)-as.numeric(daughter_rule_size)}
                        r <- r+1
                      }
                      
                      if(is.na(as.numeric(mother_rule_size))==FALSE){
                        if(is.na(as.numeric(our_new_daughter_rule_size))==FALSE){
                          if((as.numeric(mother_rule_size)+0.1) < as.numeric(our_new_daughter_rule_size)){ #Added 0.1 because otherwise you might get rounding errors.
                            #print("There is no longer room.")
                            
                            #changed from theonlymothers[p,tom].
                            pos_rem1 <- match(theonlymothers[p,tom],possible_mothers[,daughterposition])
                            pos_rem2 <- match(theonlymothers[p,tom],possible_mothers[,(daughterposition+1)])
                            
                            if(is.na(pos_rem1) == FALSE){
                              possible_mothers[pos_rem1,daughterposition] <- "0"
                            }
                            if(is.na(pos_rem2) == FALSE){
                              possible_mothers[pos_rem2,(daughterposition+1)] <- "0"
                            }
                            
                            removed <- 1
                            Event_rule_removed <- 1 #Indicating that an event has been removed.
                            
                            #mother_all <- mother_all[names(mother_all) != theonlymothers[p,tom]] #Removing the mother from the possible ones.
                            
                            #} While loop for error searching.
                          }}}
                      
                    }
                    
                    p < p+1
                  }}
                norules <- 0
                tom <- tom + 1
              }
            }
            count_replace <- count_replace+1 #Increasing this one which calculates how many columns we are into theonlymothers for this subclone. Needed when we will replace earlier chosen mothers since we later got a definitive rule.
            
          }
          
          j <- j+1
        } #Looping through motherclones.
        
        if(Event_rule_removed == 1 && not_again!=1){
          # print("We will now redo the loop.")
          Event_rule_removed <- 0
          not_again <- 1
          removed <- 0
          tom <- (tom-(count_replace-1)) #We have to reset this. changed from allocation samples to count_replace since some events are not present in all samples.
        }else{
          not_again <- 0
          i <- i+1
        }
        
        stop_while <- stop_while+1
        if(stop_while > 2*ncol(theonlymothers)){
          i <- (as.numeric(nrow(allocation_samples))+1)
          break
        }
        
        if(i==as.numeric(nrow(allocation_samples))&&again==1){ #We will redo it all in order to minimize discrepancies.
          #print("We will now redo it all.")
          tom <- 1
          again <- 0
          i <- 2
        }
        
      }else{i <- i+1} #If the event is only preent in one sample.
      
    }
    
  }
  assign("theonlymothers", theonlymothers, envir=globalenv()) #The equal clones.
  assign("equalclones_before", equalclones, envir=globalenv()) #The equal clones.
  assign("possible_mothers_new", possible_mothers, envir=globalenv()) #The equal clones.
  ###############################################
  #Updating the mother-daughter-clone allocation#
  ###############################################
  i <- 1
  j <- 1
  s <- 2
  t <- 1
  
  space <- matrix(0,50,2) #Spaces within a sample. Dynamic.
  totalspace <- matrix(0,(as.numeric(nrow(space)+1)),((2*as.numeric(ncol(sample_clone_matrix))/3)+1)) #A matrix used for calculating the spaces available.
  possible_mothers <- matrix(0,(as.numeric(nrow(space)+1)),((as.numeric(nrow(subclones))-1)*2)) #A matrix used for saving the possible motherclones.
  Not_allocated_correctly <- matrix(0,ncol(theonlymothers),3)
  nac <- 1
  
  k <- 1
  for(k in 1:(ncol(sample_clone_matrix)/3)){ #Constructing a matrix where every two columns represent a sample. The first one tells us which subclone harbors the space and the second the remaining space on top of this sample.
    totalspace[1,(2*k-1)] <- sample_clone_matrix[1,(3*k-2)]
    k <- k+1
  }
  
  k <- 1
  for(k in 2:as.numeric(nrow(subclones))){ #Constructing a matrix were every two columns represent a subclone within a sample. The first one tells us which the chosen motherclone is and the other which other possible solutions there are.
    possible_mothers[1,(2*k-3)] <- subclones[k,1]
    k <- k+1
  }
  
  subcloneswithinsample <- matrix(0,(as.numeric(nrow(sample_clone_matrix))-1),2)
  
  #SAMPLE LOOP
  i <- 1
  for(i in 1:(as.numeric(ncol(sample_clone_matrix))/3)){ #Looping through all of the samples.
    #Change the loop number
    #for(i in 1:5){
    #for(i in 3:3){
    #print("Here is i")
    #print(i)
    subcloneswithinsample <- sample_clone_matrix[2:as.numeric(nrow(sample_clone_matrix)),(3*i-2):(3*i-1)] #Extraxting the subclonal architecture and TC for a certain sample.
    subcloneswithinsample_order <- subcloneswithinsample[order(as.numeric(subcloneswithinsample[,2]),decreasing = TRUE),] #Ordering the subclones from highest to lowest TC.
    sameclones <- 0
    
    current_sample <- sample_clone_matrix[1,(3*i-2)]
    #print(c)
    
    or <- 1
    subcloneswithinsample_order <- cbind(subcloneswithinsample_order,matrix(0,nrow(subcloneswithinsample_order),1))
    for(or in 1:nrow(subcloneswithinsample_order)){
      col <- match(subcloneswithinsample_order[or,1],clone_matrix_names[,3]) #The third column contains the subclones.
      subcloneswithinsample_order[or,3] <- clone_matrix_names[2,as.numeric(col)]
      or <- or+1
    }
    
    #Arranging the subclones.
    subcloneswithinsample_order_old <- matrix(0,(as.numeric(nrow(sample_clone_matrix))-1),2)
    subcloneswithinsample_order_new <- subcloneswithinsample_order
    
    ord <- 2
    while(all(subcloneswithinsample_order_new[,1] == subcloneswithinsample_order_old[,1]) == FALSE){
      subcloneswithinsample_order_old <- subcloneswithinsample_order_new
      ord <- 2
      for(ord in 2:(as.numeric(nrow(subcloneswithinsample_order_old))-1)){ #Writing a function/loop that orders the subclones of the same size according to their median size.
        if(subcloneswithinsample_order_old[ord,2] != "0"){
          if(subcloneswithinsample_order_old[ord,2] == subcloneswithinsample_order_old[ord+1,2]){
            
            orderpos1 <- match(word(subcloneswithinsample_order_old[ord,1],2),namevector)
            orderpos2 <- match(word(subcloneswithinsample_order_old[ord+1,1],2),namevector)
            
            if(as.numeric(orderpos2) < as.numeric(orderpos1)){ 
              subcloneswithinsample_order_new <- subcloneswithinsample_order_old[c(1:(ord-1),ord+1,ord,(ord+2):nrow(subcloneswithinsample_order_old)), ]
            }
            
          }
        }
        ord <- ord+1
      }
    }
    
    subcloneswithinsample_order <- subcloneswithinsample_order_new
    #print(subcloneswithinsample_order)
    
    #SUBCLONE LOOP
    j <- 1
    for(j in 2:as.numeric(nrow(sample_clone_matrix))){ #Looping through the subclones within the sample.
      #print("Here is j")
      #print(j)
      tick <- 0
      if(j == 2){ #We're in the first position. This is the ALL-event.
        space[1,1] <- subcloneswithinsample_order[j-1,1] #The name.
        space[1,2] <- subcloneswithinsample_order[j-1,2] #The TC.
      }
      
      if(j != 2){
        if(subcloneswithinsample_order[j-1,1] != "0"){
          if(subcloneswithinsample_order[j-1,1] != "ALL"){ #We should not add it again.
            
            maxspace <- which.max(space[,2]) #Finding the largest available space.
            newname <- subcloneswithinsample_order[j-1,1] #The name of the new subclone.
            newspace <- subcloneswithinsample_order[j-1,2] #The space of the new subclone.
            full_newname <- paste(sample_clone_matrix[1,(3*i-2)],newname)
            
            #print("Precisely before conditioned.")
            if(newspace != "100"){
              if(newname %in% word(theonlymothers[1,],2,3) == TRUE){ #The clone in question has a condition on it.
                
                ######################
                #CONDITIONED SUBCLONE#
                ######################
                #print("Conditioned")
                newnamecolumn <- match(newname,word(theonlymothers[1,],2,3)) #Finding the column in theonlymothers that the daughter has.
                maxname <- theonlymothers[2,newnamecolumn] #This is the name of the mother that it has to have.
                
                #print("The maxname and the newname")
                #print(maxname)
                #print(newname)
                
                if(maxname %in% space[,1] == TRUE){ #Does the mother exist in the sample in question?
                  #print("The mother exist in the sample")
                  maxspace <- match(maxname,space[,1]) #This is the new maxspace's row position in the space matrix.
                  maxspaceTC <- space[maxspace,2] #Added this since maxspace in the row above only gives the row position and not the actual TC.
                  space[s,1] <- newname #Adding the new spacename and space size to the spacematrix.
                  space[s,2] <- newspace
                  
                  #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                  equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                  
                  e <- 0
                  for(e in 0:(equalclones_multiples-1)){
                    
                    if(e == 0){
                      e_x <- NA #The column multiple in which it exist.
                      e_y <- NA
                      x <- NA
                      y <- NA
                    }
                    
                    x_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the maxname exist.
                    y_test <- match(newname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the newname exist.
                    
                    if(is.na(x_test) == FALSE){ #The maxname exist in an equalclones column for this sample.
                      x <- (x_test+1)
                      e_x <- (i+length(samples_unique)*e) #The column where the maxname exist.
                    }
                    if(is.na(y_test) == FALSE){ #The newname exist in an equalclones column for this sample.
                      y <- (y_test+1)
                      e_y <- (i+length(samples_unique)*e) #The column where the newname exist.
                    }
                    e <- e+1
                  }
                  
                  #Added TC to the row below since we should compare TC:s. maxspace is just a row position.
                  #Added 0.1 here 200721 since rounding problems can occur when handling simulated data.
                  #print(maxspaceTC)
                  #print(newspace)
                  if(as.numeric(maxspaceTC)+0.1 >= as.numeric(newspace)){ #There must be enough room left.
                    #print("There is room")
                    
                    if(is.na(x) == FALSE){ #The maxname is in equalclones.
                      #print("Maxname is in equalclones")
                      if(is.na(y) == TRUE){ #The newname is not in equalclones.
                        #print("Newname is not")
                        
                        true_maxsize <- file_samples_subclones[match(paste(current_sample,maxname),file_samples_subclones[,13]),11]
                        true_newsize <- file_samples_subclones[match(paste(current_sample,newname),file_samples_subclones[,13]),11]
                        
                        if(true_maxsize!=true_newsize){ #Added 200920 since I got weird equalclones.
                          #if(as.numeric(newspace) != space[maxspace,2]){ #They are not of the same size.
                          #print("They are not of the same size.")
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                            if(equalclones[eq,e_x] != "0"){ #I changed the i to i+e_x
                              eq_row <- match(equalclones[eq,e_x],space[,1])
                              if(is.na(eq_row) == FALSE){
                                space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                              }
                            }
                            eq <- eq+1
                          }
                          
                        }else{ #This event should belong to the equalclones-subclone since they are of the same size.
                          #print("The new event should belong to equalclones for maxclone")
                          
                          # if(newname %in% equalclones[2:nrow(equalclones),i] == FALSE){ #If it does not already belong to the equalclones of this sample.
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){
                            if(equalclones[eq,e_x] == "0"){
                              if(newname %in% equalclones[2:nrow(equalclones),e_x] == FALSE){
                                equalclones[eq,e_x] <- newname}
                              eq <- nrow(equalclones)
                            }
                            eq <- eq+1
                          }
                        }
                      }else if(as.numeric(newspace) == as.numeric(space[maxspace,2])){ #Newname does belong to equalclones as well as maxname.
                        #print("Newname is in equalclones as well")
                        #Is it in the same equalclones?
                        
                        if(e_x == e_y){
                          # if(newname %in% equalclones[2:nrow(equalclones),i] == FALSE){ #If it does not already belong to the equalclones of this sample.
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){
                            if(equalclones[eq,e_x] == "0"){
                              if(newname %in% equalclones[2:nrow(equalclones),e_x] == FALSE){
                                equalclones[eq,e_x] <- newname}
                              eq <- nrow(equalclones)
                            }
                            eq <- eq+1
                          }
                        }
                        
                      }else{ #Maxclone and new name is in equaclones. They are not of the same size.
                        #space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
                        #Added 200202.
                        #print("Newname is in equalclones as well")
                        #print("They are not of the same size")
                        eq <- 2
                        for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                          if(equalclones[eq,e_x] != "0"){ #I changed the i to i+e_x
                            eq_row <- match(equalclones[eq,e_x],space[,1])
                            if(is.na(eq_row) == FALSE){
                              space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                            }
                          }
                          eq <- eq+1
                        }
                      }
                    }else{ #Maxname is not in equalclones.
                      #print("Maxname is not in equalclones")
                      thename <- space[maxspace,1]
                      therow <- match(thename,subcloneswithinsample_order_new[,1])
                      maxTC <- subcloneswithinsample_order_new[therow,2]
                      
                      if(as.numeric(maxTC) == as.numeric(newspace)){ #The maxspace and the newclone are of equal size.
                        #print("Maxname and newname is of equal size.")
                        if(is.na(y) == FALSE){
                          #print("Newname is in equalclones")
                          
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){
                            if(equalclones[eq,e_y] == "0"){
                              if(maxname %in% equalclones[2:nrow(equalclones),e_y] == FALSE){ #Maxname is not in this equalclone column.
                                equalclones[eq,e_y] <- maxname} #Adding the maxname to the equalclones.
                              break
                            }
                            eq <- eq+1
                          }
                          
                        }else{ #The newname does not exist in the equalclones for this sample.
                          #print("Newname is not in equalclones")
                          if(equalclones[2,i] == "0"){ #Adding them to the equalclones matrix.
                            equalclones[2,i] <- newname
                            equalclones[3,i] <- thename
                          }else if (equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){
                            #print("Adding them")
                            equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] <- newname
                            equalclones[3,(i+as.numeric(ncol(hundredpercentclones)))] <- thename
                          }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){
                            #print("Adding them")
                            equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] <- newname
                            equalclones[3,(i+2*as.numeric(ncol(hundredpercentclones)))] <- thename
                          }
                        }
                      }
                      #print(space)
                      if(as.numeric(space[maxspace,2]) >= as.numeric(newspace)){ #Added = 210412.
                        space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
                      }
                      #print(space)
                      
                    }
                  }else{ #This is the case where the conditioned subclone does not have room for the new event.
                    # print("We could not allocate it to the conditioned place.")
                    # print(maxname)
                    # print(maxspaceTC)
                    # print(newname)
                    # print(newspace)
                    
                    Not_allocated_correctly[nac,1] <- current_sample
                    Not_allocated_correctly[nac,2] <- newname
                    Not_allocated_correctly[nac,3] <- newspace
                    nac <- nac+1
                    if(theonlymothers[3,newnamecolumn]!= "0" && theonlymothers[3,newnamecolumn] %in% space[,1]){ #New 200720.
                      #This a second conditioned clone.
                      #print("Second conditioned")
                      maxspace <- match(theonlymothers[3,newnamecolumn],space[,1])
                      maxname <- space[maxspace,1]
                      maxspaceTC <- space[maxspace,2]
                      
                      #Actually we do not want the other subclones to be allocated to the event since it is not possible in all samples any more.
                      
                    }else{
                      maxspace <- which.max(space[,2]) #Finding the largest available space.
                      maxname <- space[maxspace,1]
                      maxspaceTC <- space[maxspace,2]}
                    
                    if(maxspaceTC == newspace && maxname!="ALL"){ # They are of equal size.
                      #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                      equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                      
                      e <- 0
                      for(e in 0:(equalclones_multiples-1)){
                        
                        if(e == 0){
                          e_x <- NA #The column multiple in which it exist.
                          e_y <- NA
                          x <- NA
                          y <- NA
                        }
                        
                        x_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the maxname exist.
                        y_test <- match(newname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the newname exist.
                        
                        if(is.na(x_test) == FALSE){ #The maxname exist in an equalclones column for this sample.
                          x <- (x_test+1)
                          e_x <- (i+length(samples_unique)*e) #The column where the maxname exist.
                        }
                        if(is.na(y_test) == FALSE){ #The newname exist in an equalclones column for this sample.
                          y <- (y_test+1)
                          e_y <- (i+length(samples_unique)*e) #The column where the newname exist.
                        }
                        e <- e+1
                      }
                      
                      if(is.na(x) ==TRUE){
                        #print("The maxname is not in equalclones.")
                        if(is.na(y) == TRUE){
                          #print("The newname is not in equalclones either.")
                          if(equalclones[2,i] == "0"){ #Adding them to the equalclones matris.
                            equalclones[2,i] <- newname
                            equalclones[3,i] <- maxname #Changed "thename" to "maxname" 200720.
                          }else if (equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){
                            #print("Adding them")
                            equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] <- newname
                            equalclones[3,(i+as.numeric(ncol(hundredpercentclones)))] <- maxname
                          }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){
                            #print("Adding them")
                            equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] <- newname
                            equalclones[3,(i+2*as.numeric(ncol(hundredpercentclones)))] <- maxname
                          }
                        }else{
                          #print("Newname is in equalclones.")
                          #print(newname)
                          #print(maxname)
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){
                            if(equalclones[eq,e_y] == "0"){
                              equalclones[eq,e_y] <- maxname #Adding the maxname to the equalclones.
                              break
                            }
                            eq <- eq+1
                          }
                          
                        }
                        
                      }else{
                        #print("Maxname is in equalclones.")
                        if(is.na(y) == TRUE){
                          #print("Newname is not in equalclones.")
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){
                            if(equalclones[eq,e_x] == "0"){
                              equalclones[eq,e_x] <- newname} #Adding the newname to the equalclones.
                            break
                            eq <- eq+1
                          }
                          
                        }else{
                          #print("Newname is in equalclones as well.")
                          if(e_x != e_y){
                            #print("They are not in the same equalclones.")
                            eq <- 2
                            for(eq in 2:nrow(equalclones)){
                              if(equalclones[eq,e_x] == "0"){
                                if(newname %in% equalclones[2:nrow(equalclones),e_x] == FALSE){ #Maxname is not in this equalclone column.
                                  equalclones[eq,e_x] <- newname} #Adding the maxname to the equalclones.
                                break
                              }
                              eq <- eq+1
                            }
                            equalclones[y,e_y] <- "0" #Removing the newname from its old position.
                          }
                        }
                      }
                      
                    }else{
                      #space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace))
                      eq <- 2
                      for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                        if(equalclones[eq,i] != "0"){
                          eq_row <- match(equalclones[eq,i],space[,1])
                          space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                        }
                        eq <- eq+1
                      }
                    }
                  }
                  
                }else{
                  #print("It is not allocated yet.")
                  if(maxname %in% sample_clone_matrix[,(3*i-2)] == TRUE){ #The clone exist in the sample but has not been allocated yet. This only happens if they are equal in size.
                    #print("The mother has not been allocated yet but it exist in the biopsy.")
                    
                    ###################################################################
                    #The mother has not been allocated yet but it exist in the biopsy.#
                    ###################################################################
                    
                    rowofthemother <- match(maxname,subcloneswithinsample_order[,1]) #The row in which the mother exist in the subcloneswithinsample matrix.
                    themothername <- maxname
                    thedaughtername <- newname
                    
                    if(themothername %in% word(theonlymothers[1,],2,3) == TRUE){ #The mother is conditioned.
                      #print("The mother is conditioned")
                      #print(subcloneswithinsample_order)
                      mothernamecolumn <- match(paste(current_sample,themothername),theonlymothers[1,]) #Finding the column in theonlymothers that the daughter has.
                      themothermothername <- theonlymothers[2,mothernamecolumn]
                      
                      rowofthemothermother <- match(themothermothername,subcloneswithinsample_order[,1])
                      
                      if(is.na(rowofthemothermother)==TRUE){
                        #print("The mothermother has not been allocated yet.")
                        #Added this since we do not get a tree otherwise.
                        maxspace <- which.max(space[,2])
                        rowofthemothermother <- match(space[maxspace,1],subcloneswithinsample_order[,1])
                        themothermothername <- space[maxspace,1]
                      }
                      
                    }else{
                      #print("The mother is not conditioned.")
                      rowofthemothermother <- match(space[maxspace,1],subcloneswithinsample_order[,1])
                      themothermothername <- space[maxspace,1]
                      
                    }
                    
                    tick <- 1 #Just so that we know that we've been in this loop and that space[maxspace,1] outside the loop will be the mother to the mother.
                    
                    #It may happen that the mother it should have here is smaller than the daughter. Switch positions.
                    if(as.numeric(subcloneswithinsample_order[rowofthemother,2]) < as.numeric(newspace)){
                      #print("The mother is smaller than the daughter.")
                      newspace <- subcloneswithinsample_order[as.numeric(rowofthemother),2] #Finding the new newspace.
                      rowofthemother <- match(newname,subcloneswithinsample_order[,1]) #The new row of the mother.
                      
                      temp1 <- thedaughtername
                      temp2 <- themothername
                      themothername <- temp1 #Changing the name.
                      newname <- temp2
                      
                    }
                    
                    #Mother
                    #space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(subcloneswithinsample_order[rowofthemother,2])) #The mother is allocated to its mother.
                    space[s,1] <- themothername #Adding the new spacename and space size to the spacematrix belonging to the mother.
                    space[s,2] <- subcloneswithinsample_order[rowofthemother,2]
                    maxspaceTC <- subcloneswithinsample_order[rowofthemother,2]
                    s <- s+1
                    
                    if(as.numeric(subcloneswithinsample_order[as.numeric(rowofthemothermother),2]) > as.numeric(subcloneswithinsample_order[as.numeric(rowofthemother),2])){
                      #print("The mothermother is larger than the mother.")
                      space[as.numeric(maxspace),2] <- (as.numeric(subcloneswithinsample_order[as.numeric(rowofthemothermother),2])-as.numeric(subcloneswithinsample_order[as.numeric(rowofthemother),2])) #The mother is allocated to its mother.
                      allocation_samples[match(maxname,allocation_samples[,1]),(i+1)] <- subcloneswithinsample_order[as.numeric(rowofthemothermother),1] #We have to add information about the mother's mother to the allocation matrix.
                    }else{
                      #print("The mother's mother and the mother are of equal size and should be equalclones together")
                      allocation_samples[match(maxname,allocation_samples[,1]),(i+1)] <- subcloneswithinsample_order[rowofthemothermother,1] #We have to add information about the mother's mother to the allocation matrix.
                      
                      #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                      equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                      
                      e <- 0
                      for(e in 0:(equalclones_multiples-1)){
                        
                        if(e == 0){
                          e_x <- NA #The column multiple in which it exist.
                          e_y <- NA
                          x <- NA
                          y <- NA
                        }
                        
                        x_test <- match(themothermothername,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the maxname exist.
                        y_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the newname exist.
                        
                        if(is.na(x_test) == FALSE){ #The mothermother exist in an equalclones column for this sample.
                          x <- (x_test+1)
                          e_x <- (i+length(samples_unique)*e) #The column where the maxname exist.
                        }
                        if(is.na(y_test) == FALSE){ #The mother exist in an equalclones column for this sample.
                          y <- (y_test+1)
                          e_y <- (i+length(samples_unique)*e) #The column where the newname exist.
                        }
                        e <- e+1
                      }
                      
                      
                      if(is.na(x) == FALSE){
                        #print("The mothermother is in equalclones")
                        
                        eq <- 2
                        for(eq in 2:nrow(equalclones)){
                          if(equalclones[eq,e_x] == "0"){
                            if(themothername %in% equalclones[2:nrow(equalclones),e_x] == FALSE){ #Maxname is not in this equalclone column.
                              equalclones[eq,e_x] <- themothername} #Adding the mothername to the equalclones.
                            break
                            #eq <- nrow(equalclones)
                          }
                          eq <- eq+1
                        }
                        
                      }else{ #The mothermother does not exist in the equalclones for this sample.
                        #print("The mothermother is not in equalclones yet. We now add it.")
                        if(equalclones[2,i] == "0"){ #Adding them to the equalclones matris.
                          equalclones[2,i] <- themothermothername
                          equalclones[3,i] <- themothername
                        }else if (equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] <- themothermothername
                          equalclones[3,(i+as.numeric(ncol(hundredpercentclones)))] <- themothername
                        }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] <- themothermothername
                          equalclones[3,(i+2*as.numeric(ncol(hundredpercentclones)))] <- themothername
                        }
                      }
                      
                      
                    }
                    
                    #Daughter
                    #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                    equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                    
                    e <- 0
                    for(e in 0:(equalclones_multiples-1)){
                      
                      if(e == 0){
                        e_x <- NA #The column multiple in which it exist.
                        e_y <- NA
                        x <- NA
                        y <- NA
                      }
                      
                      x_test <- match(themothermothername,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the maxname exist.
                      y_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the newname exist.
                      
                      if(is.na(x_test) == FALSE){ #The mothermother exist in an equalclones column for this sample.
                        x <- (x_test+1)
                        e_x <- (i+length(samples_unique)*e) #The column where the maxname exist.
                      }
                      if(is.na(y_test) == FALSE){ #The mother exist in an equalclones column for this sample.
                        y <- (y_test+1)
                        e_y <- (i+length(samples_unique)*e) #The column where the newname exist.
                      }
                      e <- e+1
                    }
                    
                    if(as.numeric(subcloneswithinsample_order[rowofthemother,2]) != as.numeric(newspace)){
                      #print("They are not of equal size.")
                      
                      if(is.na(x) == TRUE){
                        #print("The mother is not in equalclones")
                        space[s-1,2] <- (as.numeric(space[s-1,2])-as.numeric(newspace)) #Replacing the old maxspace.
                      }else{
                        #print("The mother is in equalclones")
                        eq <- 2
                        for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                          if(equalclones[eq,e_x] != "0"){
                            eq_row <- match(equalclones[eq,e_x],space[,1])
                            space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                          }
                          eq <- eq+1
                        }
                        
                      }
                      space[s,1] <- newname #Adding the new spacename and space size to the spacematrix belonging to the mother.
                      space[s,2] <- newspace
                      
                    }else{
                      #print("The mother and the newclone are of equal size")
                      
                      #The newclone and maxname should be in equalclones.
                      #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                      equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                      
                      e <- 0
                      for(e in 0:(equalclones_multiples-1)){
                        
                        if(e == 0){
                          e_x <- NA #The column multiple in which it exist.
                          e_y <- NA
                          x <- NA
                          y <- NA
                        }
                        
                        x_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the maxname exist.
                        y_test <- match(newname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the newname exist.
                        
                        if(is.na(x_test) == FALSE){ #The maxname exist in an equalclones column for this sample.
                          x <- (x_test+1)
                          e_x <- (i+length(samples_unique)*e) #The column where the maxname exist.
                        }
                        if(is.na(y_test) == FALSE){ #The newname exist in an equalclones column for this sample.
                          y <- (y_test+1)
                          e_y <- (i+length(samples_unique)*e) #The column where the newname exist.
                        }
                        e <- e+1
                      }
                      
                      if(is.na(y) == FALSE){
                        #print("Newname is in equalclones")
                        eq <- 2
                        for(eq in 2:nrow(equalclones)){
                          if(equalclones[eq,e_y] == "0"){
                            if(maxname %in% equalclones[2:nrow(equalclones),e_y] == FALSE){ #Maxname is not in this equalclone column.
                              equalclones[eq,e_y] <- maxname} #Adding the maxname to the equalclones.
                            eq <- nrow(equalclones)
                          }
                          eq <- eq+1
                        }
                        
                      }else{ #The newname does not exist in the equalclones for this sample.
                        #print("The newname does not exist in the equalclones for this sample.")
                        thename <- maxname #200516 ???? F?rs?kte fixa till en sak.
                        if(equalclones[2,i] == "0"){ #Adding them to the equalclones matris.
                          #print("Adding them")
                          equalclones[2,i] <- newname
                          equalclones[3,i] <- thename
                        }else if (equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] <- newname
                          equalclones[3,(i+as.numeric(ncol(hundredpercentclones)))] <- thename
                        }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] <- newname
                          equalclones[3,(i+2*as.numeric(ncol(hundredpercentclones)))] <- thename
                        }
                      }
                      
                      space[s,1] <- newname #Adding the new spacename and space size to the spacematrix belonging to the mother.
                      space[s,2] <- newspace
                    }
                    #Removing the mother from the matrix so that it is not added again.
                    rowtoremove <- match(maxname,subcloneswithinsample_order[,1])
                    #print("Here a row is removed")
                    subcloneswithinsample_order[rowtoremove,] <- "0" #210203. Testar att tysta.
                    
                    sameclones <- 1 #This is just a way to illustrate the fact that we have done this.
                  }else{ #Conditioned clone but the mother it has to have in another sample does not exist in this sample at all.
                    
                    maxspace <- which.max(space[,2]) #Finding the largest available space.
                    maxspaceTC <- as.numeric(space[maxspace,2])
                    maxname <- space[maxspace,1]
                    space[s,1] <- newname #Adding the new spacename and space size to the spacematrix.
                    space[s,2] <- as.numeric(newspace)
                    
                    # print("The conditioned mother is not present in the sample.")
                    #Looking if there is other possible places for this event to be placed
                    #which better corresponds to earlier samples.
                    other <- 1
                    nej <- 0
                    for(other in 1:j){ #j is the latest event to be placed.
                      if(as.numeric(space[other,2]) >= as.numeric(newspace)){
                        if(space[other,1] %in% allocation_samples[match(space[s,1],allocation_samples[,1]),2:i]){
                          maxspace <- other
                          maxname <- space[maxspace,1]
                          maxspaceTC <- space[maxspace,2]
                          nej <- 1
                          break
                        }
                      }
                      if(nej == 0){
                        maxname <- space[maxspace,1]
                        maxspaceTC <- space[maxspace,2]
                      }
                      other <- other+1
                    }
                    
                    #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                    equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                    
                    e <- 0
                    for(e in 0:(equalclones_multiples-1)){
                      
                      if(e == 0){
                        e_x <- NA #The column multiple in which it exist.
                        e_y <- NA
                        x <- NA
                        y <- NA
                      }
                      
                      x_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the maxname exist.
                      y_test <- match(newname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)]) #The row (minus 1) where the newname exist.
                      
                      if(is.na(x_test) == FALSE){ #The maxname exist in an equalclones column for this sample.
                        x <- (x_test+1)
                        e_x <- (i+length(samples_unique)*e) #The column where the maxname exist.
                      }
                      if(is.na(y_test) == FALSE){ #The newname exist in an equalclones column for this sample.
                        y <- (y_test+1)
                        e_y <- (i+length(samples_unique)*e) #The column where the newname exist.
                      }
                      e <- e+1
                    }
                    
                    # x <- match(space[maxspace,1],equalclones[2:nrow(equalclones),i])
                    # y <- match(newname,equalclones[2:nrow(equalclones),i])
                    
                    if(maxname != "ALL"){
                      if(as.numeric(maxspaceTC) == as.numeric(newspace)){
                        if(is.na(x) == FALSE){ #The maxname is in equalclones.
                          #print("The maxname is in equalclones")
                          if(is.na(y) == TRUE){ #The newname is not.
                            #print("The newname is not in equalclones")
                            
                            eq <- 2
                            for(eq in 2:nrow(equalclones)){
                              if(equalclones[eq,e_x] == "0"){
                                equalclones[eq,e_x] <- newname #Adding the newname to the equalclones for maxname.
                                break
                              }
                              eq <- eq+1
                            }
                            
                          }else{
                            #print("The newname is in equalclones as well")
                            if(as.numeric(newspace) != space[maxspace,2]){
                              eq <- 2
                              for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                                if(equalclones[eq,e_y] != "0"){ #Changed i to e_y.
                                  eq_row <- match(equalclones[eq,e_y],space[,1])
                                  space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                                }
                                eq <- eq+1
                              }
                            }
                          }
                        }
                      }else{
                        if(is.na(x) == FALSE){
                          #print("The maxname is in equalclones")
                          eq <- 2
                          for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                            if(equalclones[eq,e_x] != "0"){ #Changed i to e_x
                              eq_row <- match(equalclones[eq,e_x],space[,1])
                              space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                            }
                            eq <- eq+1
                          }
                        }else{
                          space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
                        }
                      }
                      
                    }else{
                      #print("The maxname is ALL")
                      space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace))
                    }
                    
                  }
                }
                ############################
                #NOT A CONDITIONED SUBCLONE#
                ############################
              }else{
                #print("Not conditioned")
                space[s,1] <- newname #Adding the new spacename and space size to the spacematrix.
                space[s,2] <- newspace
                
                #Seeing if there is other possible places for this event to be placed
                #which better corresponds to earlier samples.
                other <- 1
                nej <- 0
                for(other in 1:j){ #j is the latest event to be placed.
                  if(as.numeric(space[other,2]) >= as.numeric(newspace)){
                    if(space[other,1] %in% allocation_samples[match(space[s,1],allocation_samples[,1]),2:i]){
                      maxspace <- other
                      maxname <- space[maxspace,1]
                      maxspaceTC <- space[maxspace,2]
                      nej <- 1
                      break
                    }
                  }
                  if(nej == 0){
                    maxname <- space[maxspace,1]
                    maxspaceTC <- space[maxspace,2]
                  }
                  other <- other+1
                }
                
                #I try to take into account that we actually have multiple columns of equalclones for one single sample.
                equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
                
                e <- 0
                for(e in 0:(equalclones_multiples-1)){
                  
                  if(e == 0){
                    e_x <- NA #The column multiple in which it exist.
                    e_y <- NA
                    x <- NA
                    y <- NA
                  }
                  
                  x_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)])
                  y_test <- match(newname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)])
                  
                  if(is.na(x_test) == FALSE){ #The maxname exist in an equalclones column for this sample.
                    x <- (x_test+1) #Row
                    e_x <- (i+length(samples_unique)*e) #Column
                  }
                  if(is.na(y_test) == FALSE){ #The newname exist in an equalclones column for this sample.
                    y <- (y_test+1)
                    e_y <- (i+length(samples_unique)*e)
                  }
                  e <- e+1
                }
                
                #x <- match(maxname,equalclones[2:nrow(equalclones),i])
                #y <- match(newname,equalclones[2:nrow(equalclones),i])
                
                if(is.na(x) == FALSE){ #The maxname is in equalclones.
                  #print("Maxname is in equalclones")
                  if(is.na(y) == TRUE){ #The newname is not.
                    #print("Newname is not in equalclones")
                    
                    if(as.numeric(newspace) == as.numeric(subcloneswithinsample_order[maxspace,2])){
                      #print("They are are of equal size")
                      eq <- 2
                      for(eq in 2:nrow(equalclones)){
                        if(equalclones[eq,e_x] == "0"){
                          if(newname %in% equalclones[2:nrow(equalclones),e_x] == FALSE){
                            equalclones[eq,e_x] <- newname}
                          eq <- nrow(equalclones)
                        }
                        eq <- eq+1
                      }
                      
                    }else{
                      #print("They are not equal in size.")
                      #Removing the TC from all of the events in the equalclones.
                      eq <- 2
                      breaking <- 1
                      for(eq in 2:nrow(equalclones)){
                        if(equalclones[eq,e_x] != "0"){
                          eq_row <- match(equalclones[eq,e_x],space[,1])
                          space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                        }
                        eq <- eq+1
                      }}
                    
                  }else{
                    #print("Newname is in equalclones as well")
                    
                    #if(as.numeric(newspace) == as.numeric(space[maxspace,2])){ #Here we just compare the newspace with the space the maxspace has right now. They might not be equal.
                    #We cannot go around adding events to equalclones just based on this above. Changed it to the line below instead 200308.
                    #if(as.numeric(newspace) == as.numeric(subcloneswithinsample_order[maxspace,2])){
                    if(as.numeric(newspace) == as.numeric(subcloneswithinsample_order[maxspace,2])){
                      #print("They are are of equal size")
                      if(newname %in% equalclones[2:nrow(equalclones),e_x] == FALSE){ #If it does not already belong to the equalclones of this sample.
                        eq <- 2
                        
                        for(eq in 2:nrow(equalclones)){
                          if(equalclones[eq,e_x] == "0"){
                            if(newname %in% equalclones[2:nrow(equalclones),e_x] == FALSE){
                              equalclones[eq,e_x] <- newname}
                            eq <- nrow(equalclones)
                          }
                          eq <- eq+1
                        }
                      }
                    }else{
                      #space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
                      equalclones[y,e_y] <- "0" #Removing the newname from equalclones.
                      
                      eq <- 2
                      for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                        if(equalclones[eq,e_x] != "0"){
                          eq_row <- match(equalclones[eq,e_x],space[,1])
                          space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                        }
                        eq <- eq+1
                      }
                    }
                  }
                }else{
                  #print("Maxname is not in equalclones")
                  
                  if(maxname != "ALL"){
                    maxnamerow <- match(maxname,subcloneswithinsample_order[,1])
                    maxspaceTC <- subcloneswithinsample_order[maxnamerow,2]
                    # print(subcloneswithinsample_order)
                    # print(maxname)
                    # print(maxspaceTC)
                    # print(newname)
                    # print(newspace)
                    if(as.numeric(newspace) == as.numeric(maxspaceTC)){
                      #print("Newclone and maxclone are of the same size.")
                      
                      if(is.na(y) == FALSE){
                        #print("Newclone is in equalclones.")
                        
                        equalclones[y,e_y] <- "0" #Removing it from its old place and placing it on a new one.
                        
                        if(equalclones[2,i] == "0"){ #Adding them to the equalclones matris.
                          equalclones[2,i] <- newname
                          equalclones[3,i] <- maxname
                        }else if (equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] <- newname
                          equalclones[3,(i+as.numeric(ncol(hundredpercentclones)))] <- maxname
                        }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] <- newname
                          equalclones[3,(i+2*as.numeric(ncol(hundredpercentclones)))] <- maxname
                        }
                        
                      }else{
                        #print("Neither the maxclone nor the newclone is in equalclones.")
                        
                        if(equalclones[2,i] == "0"){ #Adding them to the equalclones matris.
                          equalclones[2,i] <- newname
                          equalclones[3,i] <- maxname
                        }else if (equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+as.numeric(ncol(hundredpercentclones)))] <- newname
                          equalclones[3,(i+as.numeric(ncol(hundredpercentclones)))] <- maxname
                        }else if(equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] == "0"){
                          #print("Adding them")
                          equalclones[2,(i+2*as.numeric(ncol(hundredpercentclones)))] <- newname
                          equalclones[3,(i+2*as.numeric(ncol(hundredpercentclones)))] <- maxname
                        }
                      }
                    }else{
                      space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
                      if(is.na(y) == FALSE){
                        #print("Newname is in equalclones")
                        equalclones[y,e_y] <- "0" #Removing it.
                      }
                      
                    }
                  }else{
                    #print("Maxname is ALL.")
                    if(is.na(y) == FALSE){
                      #print("The daughter exist in an equalclones.")
                      equalclones[y,e_y] <- "0" #Removing it
                    }
                    
                    space[maxspace,2] <- (as.numeric(space[maxspace,2])-as.numeric(newspace)) #Replacing the old maxspace.
                  }
                }
              }
              
            }else{ #The newspace is 100.
              if(space[1,1] == "ALL"){ #The ALL-space is now occupied by the new subclone.
                space[1,2] <- "0"
              }
              space[s,1] <- newname #Adding the new spacename and space size to the spacematrix. We do not have to alter anything when dealing with hundredpercentclones.
              space[s,2] <- newspace
            } #100 %
            
          }else{ #ALL
          }
        }else{ #"0"
        }
        
        #The point of these loops are to take into account cases where the newname is in an equalclone situation with
        #an event in another sample and this event is also present in this sample but they are not equal here.
        if(newname %in% equalclones == TRUE){ #The clone exists in equalcones.
          
          #I try to take into account that we actually have multiple columns of equalclones for one single sample.
          equalclones_multiples <- as.numeric(ncol(equalclones))/length(samples_unique)
          
          e <- 0
          for(e in 0:(equalclones_multiples-1)){
            
            if(e == 0){
              e_x <- NA #The column multiple in which it exist.
              e_y <- NA
              x <- NA
              y <- NA
            }
            
            x_test <- match(maxname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)])
            y_test <- match(newname,equalclones[2:nrow(equalclones),i+(length(samples_unique)*e)])
            
            if(is.na(x_test) == FALSE){ #The maxname exist in an equalclones column for this sample.
              x <- (x_test+1) #Row
              e_x <- (i+length(samples_unique)*e) #Column
            }
            if(is.na(y_test) == FALSE){ #The newname exist in an equalclones column for this sample.
              y <- (y_test+1)
              e_y <- (i+length(samples_unique)*e)
            }
            e <- e+1
          }
          
          if(is.na(y) == TRUE && maxname != "ALL"){ #But not in this sample. Added ALL 201108.
            if(equalclones[2,i] != "0"){ #There exists equalclones in this sample.
              equalclonename <- equalclones[2,i]
              equalclonepos <- match(equalclonename,subcloneswithinsample_order_new) #The position for the equalclone.
              theTCforequalpos <- subcloneswithinsample_order_new[equalclonepos,2] #Changed to new.
              
              if(as.numeric(theTCforequalpos) == as.numeric(newspace) && subcloneswithinsample_order_new[match(maxname,subcloneswithinsample_order_new[,1]),2]==subcloneswithinsample_order_new[match(newname,subcloneswithinsample_order_new[,1]),2]){ #The new clone and the clone in equalclones are equal in size.
                #print("Den kom in hit")
                h <- 1
                for(h in 1:ncol(equalclones)){ #Looping through the columns.
                  n <- 1
                  u <- 1
                  t <- 2
                  for(n in 1:nrow(equalclones)){ #Looping through the rows.
                    if(equalclones[n,h] == equalclonename){ #We have found the equalclone in a particular sample.
                      u <- h
                    }
                    if(equalclones[n,h] == newname){ #We have found the newclone.
                      t <- h
                    }
                    if(t == u){ #The equalclone and the newclone actually exists together as equalclones in another sample.
                      o <- 1
                      
                      for(o in 1:length(equalclones[,i])){ #We add the newclone to the equalclone.
                        if(o == 1){
                          #print("They exist together in another sample!")
                        }
                        
                        if(equalclones[o,i] == "0"){
                          if(newname %in% equalclones[,i] == FALSE){
                            equalclones[o,i] <- newname
                            
                            eq <- 2
                            for(eq in 2:nrow(equalclones)){ #I have to reduce the space for all of the events belonging to this subclone.
                              if(equalclones[eq,i] != "0"){ #Changed it to i instead of e_x since we did not have e_x for RMS6 B2. 200329.
                                eq_row <- match(equalclones[eq,i],space[,1])
                                space[eq_row,2] <- (as.numeric(space[eq_row,2])-as.numeric(newspace))
                              }
                              eq <- eq+1
                            }
                            
                            #space[maxspace,2] <- (as.numeric(space[maxspace,2])+as.numeric(newspace)) #Resetting the space.
                            maxspace <- match(equalclonename, space[,1])
                            
                            o <- nrow(equalclones)
                          }
                        }
                        o <- o+1
                      }
                      
                    }
                    n <- n+1
                  }
                  h <- h+1
                }
              }
            }
          }
        }
        
        #print(space[maxspace,1])
        #print(newname)
        #print(allocation_samples[match(newname,allocation_samples[,1]),(i+1)])
        if(tick == 0){
          if(space[s,1] != "0"){
            if(sameclones != 1){
              #Treating the case when space[maxspace,2] = 100 % and the newspace as well. Then the motherclone and the daughterclone are both part of the base.
              if(subcloneswithinsample_order[j-1,2] == "100"){
                
                if(subcloneswithinsample_order[j-1,1] != "ALL"){
                  if(subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),2] == "100"){
                    if(subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1] != "ALL"){
                      
                      newspace_name <- subcloneswithinsample_order[j-1,1]
                      maxspace_name <- subcloneswithinsample_order[match(space[maxspace,1],subcloneswithinsample_order[,1]),1]
                      
                      allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- maxspace_name #Annoting the mother clone of each of the subclones within each sample.
                      allocation_samples[match(maxspace_name,allocation_samples[,1]),(i+1)] <- newspace_name #Annoting the mother clone of each of the subclones within each sample.
                      
                      
                    }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
                    }
                  }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
                  }
                }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
                }
              }else{allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- space[maxspace,1] #Annoting the mother clone of each of the subclones within each sample.
              }
              
            }else{
              #if(themothername != "ALL"){
              #allocation_samples[match(thedaughtername,allocation_samples[,1]),(i+1)] <- themothername}else{
              allocation_samples[match(newname,allocation_samples[,1]),(i+1)] <- space[maxspace,1] #200330
              #}
            }
          } #"0"
        }else{
          allocation_samples[match(space[s,1],allocation_samples[,1]),(i+1)] <- themothername
        }
      } #j != 1.#}
      
      #print(space)
      if(j == as.numeric(nrow(sample_clone_matrix))){ #We're at the end of a sample.
        #print("Totalspace is added.")
        totalspace[2:as.numeric(nrow(totalspace)),((i*2)-1):((i*2))] <- space
        #t <- t+1
        
        s <- 2 #Resetting s and space.
        space <- matrix(0,50,2)
        
      }else{s <- s+1}
      
      j <- j+1
    }
    
    i <- i+1
    
  }
  
  assign("totalspace", totalspace, envir=globalenv()) #The equal clones.
  assign("equalclones_after", equalclones, envir=globalenv()) #The equal clones.
  assign("allocation_samples_updated", allocation_samples, envir=globalenv()) #The mother-daughter division are exported to the global environment.
  assign("Not_allocated_correctly",Not_allocated_correctly,envir=globalenv())
  
  #Treating cases where the alteration did not get placed in its correct mother.
  nac <- Not_allocated_correctly[as.numeric(Not_allocated_correctly[,3])<10,]
  i <- 1
  for(i in 1:nrow(nac)){
    
    col <- match(nac[i,1],allocation_samples[1,])
    row <- match(nac[i,2],allocation_samples[,1])
    
    mothercol <- match(paste(nac[i,1],nac[i,2]),theonlymothers[1,])
    mothercond <- theonlymothers[2,mothercol]
    
    allocation_samples[row,col] <- mothercond
    
    i <- i+1
  }
  
  assign("allocation_samples_revised", allocation_samples, envir=globalenv()) 
  
  #Cleaning the equalclones matrix.
  if(length(equalclones[1:10,equalclones[2,]!="0"]) > 1){
    equalclones_eq <- equalclones[1:20,equalclones[2,]!="0"]
    if(is.null(ncol(equalclones_eq))==FALSE){
      i <- 1
      for(i in 1:ncol(equalclones_eq)){
        sample <- equalclones_eq[1,i]
        all_col <- match(sample,allocation_samples_updated[1,])
        
        j <- 2
        for(j in 2:nrow(equalclones_eq)){
          yes <- 0
          subclone <- equalclones_eq[j,i]
          
          if(subclone!="0"){
            subclone_row <- match(subclone,allocation_samples_updated[,1])
            
            mother <- allocation_samples_updated[subclone_row,all_col]
            if(mother %in% equalclones_eq[,i] == FALSE){
              #print("This subclone did not get placed in any of the equalclones.")
              
              k <- 3
              for(k in 3:nrow(allocation_samples_updated)){
                if(allocation_samples_updated[k,all_col] == subclone){
                  subclone_mother <- allocation_samples_updated[k,1]
                  if(subclone_mother != "0"){
                    if(subclone_mother %in% equalclones_eq[,i] == TRUE){
                      #print("This subclone got an equalclone as daughter.")
                      yes <- 1
                    }
                  }
                  
                }
                k <- k+1
              }
              
            }else{
              #print("The mother is in equalclones.")
              yes <- 1
            }
          }
          
          if(yes == 0 && equalclones_eq[j,i] != "0"){
            #print("This subclone did not get placed in any of these equalclones nor got any of them as daughter.")
            rm <- match(paste(equalclones_eq[1,i],equalclones_eq[j,i]),file_samples_subclones[,13])
            rs <- file_samples_subclones[rm,11]
            if(as.numeric(rs) <= 50){
              equalclones_eq[j,i] <- "0"
              if(j != nrow(equalclones_eq) && as.numeric(rs) <= 50){ #If under 50 they have to be inside each other.
                equalclones_eq[j:(nrow(equalclones_eq)-1),i] <- equalclones_eq[(j+1):nrow(equalclones_eq),i]
              }
              yes <- 0}
          }
          
          j <- j+1
        }
        i <- i+1
      }
      assign("equalclones_cleaned", equalclones_eq, envir=globalenv()) #The equal clones.
      
      equalclones <- equalclones_eq
      
    }
  }
  
  ###################################################################
  #Treating the case when we have many clones which have the same TC#
  ###################################################################
  i <- 2
  j <- 1
  therows <- nrow(equalclones)
  for(j in 1:as.numeric(ncol(equalclones))){ #Looping through the samples.
    v <- equalclones[,j] #Extracting the column.
    v <- cbind(v,matrix(0,nrow(equalclones),1)) #Adding a new column to it. Changed from "rowsofhundred" to nrow(equalclones).
    clonenumber <- length(v[v!="0"])
    v <- v[v[,1] !="0",]
    if(clonenumber > 2){ #If we have more than 1 subclone within the sample with equal TC % who are placed inside each other.
      k <- 2
      for(k in 2:clonenumber){ #Giving alterations in all combinations.
        mother <- paste(v[1,1],v[k,1])
        l <- 2
        for(l in 2:clonenumber){ #Looping through the clones.
          if(k != l){
            daughter <- paste(v[1,1], v[l,1])
            mothercolumn <- match(mother,eventmatrix[1,])
            daughtercolumn <- match(daughter,eventmatrix[1,])
            m <- 2
            for(m in 2:as.numeric(nrow(eventmatrix))){ #Looping through the events.
              if(eventmatrix[m,mothercolumn] == "1"){
                eventmatrix[m,daughtercolumn] <- "1"
              }else if(eventmatrix[m,daughtercolumn] == "1"){ #Added 200821.
                eventmatrix[m,mothercolumn] <- "1"
              }
              m <- m+1
            }
            
          }
          l <- l+1
        }
        k <- k+1
      }
      
    }
    j <- j+1
  }
  
  ##########################################################
  #Giving the daughter subclones their motherclone's events#
  ##########################################################
  motherdaughterevent <- matrix(0,as.numeric(nrow(allocation_samples)),3)
  motherdaughterevent[,1] <- allocation_samples[,1]
  motherdaughterevent_new <- matrix(0,as.numeric(nrow(allocation_samples)),3) #The one used for the order.
  
  eventmatrix_original <- eventmatrix
  motherdaughter <- matrix(0,50,2)
  themedian <- as.numeric(nrow(medianmatrix))
  
  j <- 2
  s <- 1
  for(j in 2:as.numeric(ncol(allocation_samples))){ #Looping through the columns.
    motherdaughterevent[,2] <- allocation_samples[,j] #Extracting that particular column.
    motherdaughterevent[1,3] <- "100" #Setting the subclone name to "100" just so that it will stay where it is.
    
    i <- 2
    for(i in 2:as.numeric(nrow(allocation_samples))){ #Looping through the rows.
      
      if(motherdaughterevent[i,2] != "ALL"){
        if(motherdaughterevent[1,2] != "ALL"){
          mothername <- paste(motherdaughterevent[1,2],motherdaughterevent[i,2])
          daughtername <- paste(motherdaughterevent[1,2],motherdaughterevent[i,1])
        }else{
          mothername <- "0"
          daughtername <- "0"}
      }else{mothername <- "ALL"
      daughtername <- "ALL"}
      
      if(motherdaughterevent[i,2] != "0"){ #Adding the median TC of the mother clone to the matrix.
        samplecolumn <- match(mothername,medianmatrix[1,])
        theTC <- medianmatrix[themedian,samplecolumn]
        samplecolumn_daughter <- match(daughtername,medianmatrix[1,])
        theTC_daughter <- medianmatrix[themedian,samplecolumn_daughter]
      }else{theTC <- "0"
      theTC_daughter <- "0"}
      
      motherdaughterevent[i,3] <- theTC_daughter #The median TC for that subclone in that particular sample.
      
      i <- i+1
    }
    
    motherdaughter_order_before <- motherdaughterevent[order(as.numeric(motherdaughterevent[1:nrow(motherdaughterevent),3]), decreasing = TRUE),]
    motherdaughter_totalspace <- as.matrix(totalspace[,(j-1)*2-1]) #The order obtained from totalspace.
    motherdaughter_totalspace <- as.matrix(motherdaughter_totalspace[motherdaughter_totalspace != "0",]) #Removing the zero rows.
    
    a <- 2
    for(a in 2:as.numeric(nrow(motherdaughter_totalspace))){ #Looping through the order we should have.
      
      clone <- motherdaughter_totalspace[a,1]
      clonerow <- match(clone,motherdaughterevent)
      
      motherdaughterevent_new[a,] <- motherdaughterevent[clonerow,] #New matrix in which the correct order is saved.
      
      a <- a+1
    }
    
    motherdaughter_order <- as.matrix(motherdaughterevent_new)
    motherdaughter_order[1,2] <- motherdaughterevent[1,2] #In this position we want to to have the sample name.
    #print(motherdaughter_order)
    i <- 2
    for(i in 2:as.numeric(nrow(allocation_samples))){ #Looping through the subclones.
      eq <- 0
      
      if(motherdaughter_order[i,2] != "0"){ #The subclone does not exist in that particular sample.
        if(motherdaughter_order[i,2] != "ALL"){ #The subclone has already gotten these alterations.
          
          daughter_name <- paste(motherdaughter_order[1,2],motherdaughter_order[i,1]) #The name of the subclone in a particular sample.
          daughter_column <- match(daughter_name,eventmatrix[1,]) #The corresponding column in the eventmatrix.
          
          mother_name <- paste(motherdaughter_order[1,2],motherdaughter_order[i,2]) #The motherclone.
          mother_column <- match(mother_name,eventmatrix[1,]) #The column in the eventmatrix corresponding to the motherclone.
          
          col <- match(word(daughter_name,1),equalclones_new[1,])
          
          if(is.na(col) == FALSE){
            if(word(daughter_name,2,3) %in% equalclones_new[,col]){
              #It is part of equalclones in this sample. Then the other ones here should also have this mother.
              #print("It is in equalclones.")
              eq <- 1
            }
          }
          
          k <- 2
          for(k in 2:as.numeric(nrow(eventmatrix))){
            if(eventmatrix[k,mother_column] == "1"){
              eventmatrix[k,daughter_column] <- "1"
              
              if(eq == 1){#Equalclones.
                l <- 2
                sub <- equalclones_new[equalclones_new[,col]!="0",col]
                for(l in 2:length(sub)){
                  eq_column <- match(paste(sub[1],sub[l]),eventmatrix[1,])
                  eventmatrix[k,eq_column] <- "1"
                  l <- l+1
                }
              }
              
              
            }
            k <- k+1
          }
        }
      }
      i <- i+1
    }
    motherdaughterevent <- matrix(0,as.numeric(nrow(allocation_samples)),3)
    motherdaughterevent[,1] <- allocation_samples[,1]
    motherdaughterevent_new <- matrix(0,as.numeric(nrow(allocation_samples)),3) #The one used for the order.
    
    j <- j+1
  }
  
  ###################################################################
  #Treating the case when we have many clones which have the same TC#
  ###################################################################
  i <- 2
  j <- 1
  therows <- nrow(equalclones)
  for(j in 1:as.numeric(ncol(equalclones))){ #Looping through the samples.
    v <- equalclones[,j] #Extracting the column.
    v <- cbind(v,matrix(0,nrow(equalclones),1)) #Adding a new column to it. Changed from "rowsofhundred" to nrow(equalclones).
    clonenumber <- length(v[v!="0"])
    v <- v[v[,1] !="0",]
    if(clonenumber > 2){ #If we have more than 1 subclone within the sample with equal TC % who are placed inside each other.
      k <- 2
      for(k in 2:clonenumber){ #Giving alterations in all combinations.
        mother <- paste(v[1,1],v[k,1])
        l <- 2
        for(l in 2:clonenumber){ #Looping through the clones.
          if(k != l){
            daughter <- paste(v[1,1], v[l,1])
            mothercolumn <- match(mother,eventmatrix[1,])
            daughtercolumn <- match(daughter,eventmatrix[1,])
            m <- 2
            for(m in 2:as.numeric(nrow(eventmatrix))){ #Looping through the events.
              if(eventmatrix[m,mothercolumn] == "1"){
                eventmatrix[m,daughtercolumn] <- "1"
              }else if(eventmatrix[m,daughtercolumn] == "1"){ #Added 200821.
                eventmatrix[m,mothercolumn] <- "1"
              }
              m <- m+1
            }
            
          }
          l <- l+1
        }
        k <- k+1
      }
      
    }
    j <- j+1
  }
  
  eventmatrix_new <- matrix(0,(as.numeric(nrow(eventmatrix))-1), (as.numeric(ncol(eventmatrix))-1)) #Skapar en ny h?ndelsematris d?r vi bara har med 1:orna och 0:orna.
  eventmatrix_new <- eventmatrix[2:as.numeric(nrow(eventmatrix)),2:as.numeric(ncol(eventmatrix))]
  eventmatrix_new <- as.matrix(eventmatrix_new)
  rownames(eventmatrix_new) <- eventmatrix[2:as.numeric(nrow(eventmatrix)),1] #L?gger till radnamnen och kolumnnamnen till den nya matrisen.
  colnames(eventmatrix_new) <- eventmatrix[1,2:as.numeric(ncol(eventmatrix))]
  eventmatrix_new <- t(eventmatrix_new)
  
  stop.time <- Sys.time()
  print("Execution time")
  print(stop.time-start.time)
  return(eventmatrix_new)
}

#Splitting the input file.
splitdata <- function(file,name){
  k <- 1
  s <- 1
  samples <- matrix(0,40,2)
  rownames(samples) <- c(1:40)
  file <- file[is.na(file[,1])==FALSE,]
  for(k in 1:as.numeric(nrow(file))){ #Looping over all samples.
    if(k == 1){ #The first position.
      samples[s,1] <- k
      rownames(samples)[s] <- file[k,1]
    }
    if(k != 1){ #Every other position.
      if(file[k-1,1] != file[k,1]){
        if(k != nrow(file)){
          samples[s,2] <- k-1 #End position.
          s <- s+1
          samples[s,1] <- k
          rownames(samples)[s] <- file[k,1]}}
    }
    if(k == nrow(file)){ #Last row.
      if(file[k-1,1] != file[k,1]){
        samples[s,2] <- k-1 #End position.
        s <- s+1
        samples[s,1] <- k
        samples[s,2] <- k
        rownames(samples)[s] <- file[k,1]
      }else{
        samples[s,2] <- k}
    }
    k <- k+1
  }
  
  i <- 1
  for(i in 1:nrow(samples)){ #Localizing that particular tumor in the sample file.
    if(samples[i,1] != 0){
      tumorname <- match(name,rownames(samples))
    }
    i <- i+1
  }
  datasegment <- file[samples[tumorname,1]:samples[tumorname,2],] #Extracting the data for that particular tumor from the large segment file.
  
  return(datasegment)
}

#Adding a stem to the data.
stem <- function(eventmatrix,co,root){
  i = 1
  j = 1
  s = 1
  class(eventmatrix) <- "numeric"
  eventmatrix_new <- eventmatrix
  
  if(root == "Stem"){
    stemroot <- matrix(0, 1, as.numeric(ncol(eventmatrix_new)))
    i = 1
    for(i in 1:as.numeric(ncol(eventmatrix_new))){
      if(sum(as.numeric(eventmatrix_new[,i]))/as.numeric(nrow(eventmatrix_new)) == 1){
        stemroot[1,i] <- 1
        i <- i+1
      }
    }
    eventmatrix_new <- rbind(eventmatrix_new,stemroot)
    rownames(eventmatrix_new)[nrow(eventmatrix_new)] <- "Stem"}
  
  if(root == "Normal"){
    M <- matrix(0, 1, as.numeric(ncol(eventmatrix_new)))
    eventmatrix_new <- rbind(eventmatrix_new,M)
    rownames(eventmatrix_new)[as.numeric(nrow(eventmatrix_new))] <- "Normal"
    i <- i+1}
  
  if(root == "None"){
    eventmatrix_new <- eventmatrix_new
  }
  return(eventmatrix_new)
}

#Transform the file into phyDat format.
phydatevent <- function(excelfil){
  patient.matrix <- as.matrix(excelfil)
  patient.phydat <- phyDat(patient.matrix,type="USER",levels=c(0,1),ambiguity='0')
  return(patient.phydat)
}

#Maximum Likelihood
ml_tree <- function(Eventmatrix,root) {
  dm_h <- dist.hamming(Eventmatrix)
  starting_tree <- NJ(dm_h)
  starting_tree <- root(starting_tree, outgroup = root,resolve.root = TRUE)
  Lf <- pml(starting_tree, Eventmatrix) #Obtaining an object of class pml
  Lf_JC <- optim.pml(Lf, model = "JC", optEdge = TRUE)
  return(Lf_JC)
}

#Maximum parsimony
mp_tree <- function(Eventmatrix,root){
  MP_tree_pratchet <- pratchet(Eventmatrix, start = NULL, method = "fitch", maxit = 2000, k = 10, #Funktionen anv?nder sig av Fithchs algoritm. Pratchet = p ratchett (1999).
                               trace = 1, all = FALSE, rearrangements = "TBR",
                               perturbation = "ratchet") #Den ger det b?sta tr?det den funnit. Minimerar parsimony score.
  
  MP_tree_pratchet <- root(MP_tree_pratchet, outgroup = root,resolve.root = TRUE)
  
  treeRatchet <- acctran(MP_tree_pratchet, Eventmatrix) #Gives us the tree with an edge length fulfilling the acctran criterion.
  return(treeRatchet)
}

#Visualising the MP-tree.
MP_treeplot <- function(MP_tree,limitmp){
  EM_testmp <- ggplot(MP_tree) + geom_tree() + geom_tiplab(size=4, color = "black") #+ geom_treescale(width = 1)
  EM_testmp <- EM_testmp +  theme_tree() + limitmp
  mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
  print(EM_testmp+mytheme)
  return(EM_testmp)}

#Visualising the ML-tree.
ML_treeplot <- function(ML_tree,limitml){
  EM_test_mltree <- ML_tree$tree
  EM_test_mltree <- ggplot(EM_test_mltree) + geom_tree() + geom_tiplab(size=4, color = "black") #+ geom_treescale(width = 1)
  EM_test_mltree <- EM_test_mltree +  theme_tree() + limitml
  mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
  print(EM_test_mltree+mytheme)
  return(EM_test_mltree)}

#Making new subclones.
subclones <- function(EM_test,file_samples_subclones,root){
  if(missing(root)==TRUE){root <- "Normal"} #The default is to root the tree in a normal cell.
  EM_newnames <- unique(EM_test) #Finding all unique rows in the EM i.e. all subclones that have different sets of mutations.
  clonenames_new <- matrix(0,(as.numeric(nrow(EM_newnames))*2),50) #Creating a new matrix that will contain the new subclone names and which former subclones it includes.
  
  samples_all <- t(as.matrix(unique(file_samples_subclones[,2]))) #A matrix containing all unique samples.
  samples <- t(as.matrix(samples_all[samples_all != "ALL"]))
  sampleTC <- matrix(0,1,ncol(samples))
  sampleTC[1,1:ncol(samples)] <- "100"
  
  i <- 1
  l <- 2
  k <- 1
  for(i in 1:as.numeric(nrow(EM_newnames))){ #Looping through each of the unique subclones.
    uniquenamerow <- match(EM_newnames[i,1],file_samples_subclones[,13]) #Finding the position of the subclone.
    uniquenameTC <- file_samples_subclones[uniquenamerow,11] #Finding the TC of the subclone.
    
    j <- 1
    for(j in 1:as.numeric(nrow(EM_test))){ #Every unique subclone is to be compared with the others.
      
      if(all(EM_newnames[i,] == EM_test[j,]) == TRUE){ #They have to include the same events.
        clonenames_new[k,l] <- rownames(EM_test)[j] #Saving the subclone name to the matrix.
        
        if(rownames(EM_test)[j] != "ALL"){
          column <- match(word(rownames(EM_test)[j],1),sample_clone_matrix[1,])
          row <- match(word(rownames(EM_test)[j],2,3),sample_clone_matrix[,column])
          theTC <- sample_clone_matrix[row,(column+1)] #Finding the TC for the subclone.
        }else{theTC <- "100"}
        
        clonenames_new[(k+1),l] <- theTC #Saving the TC below its subclone name in the matrix.
        l <- l+1
      }
      
      j <- j+1
    }
    
    l <- 2
    k <- k + 2
    i <- i+1
  }
  
  m <- 1
  for(m in 1:(as.numeric(nrow(clonenames_new))/2)){ #Calculating the mean of all of the subclones within the new subclones.
    #print(sum(as.numeric(clonenames_new[2*m,2:ncol(clonenames_new)])))
    clonenames_new[2*m,1] <- mean(as.numeric(clonenames_new[2*m,clonenames_new[2*m,] != 0]))
    clonenames_new[(2*m-1),1] <- mean(as.numeric(clonenames_new[2*m,clonenames_new[2*m,] != 0]))
    m <- m+1
  }
  
  #Giving the new subclone names. The order is determined based on the subclones' mean TC:s.
  clonenames_new_order <- clonenames_new[order(as.numeric(clonenames_new[,1]), decreasing = TRUE),] #Ordering the subclones based on their TC.
  newnames <- c("Subclone A", "Subclone B","Subclone C","Subclone D","Subclone E","Subclone F","Subclone G","Subclone H","Subclone I","Subclone J","Subclone K","Subclone L","Subclone M","Subclone N","Subclone O","Subclone P", "Subclone Q","Subclone R","Subclone S","Subclone T","Subclone U","Subclone V","Subclone X","Subclone Y","Subclone Z",
                "Subclone ZA","Subclone ZB","Subclone ZC", "Subclone ZD","Subclone ZE", "Subclone ZF", "Subclone ZG", "Subclone ZH", "Subclone ZI", "Subclone ZJ", "Subclone ZK", "Subclone ZL", "Subclone ZM", "Subclone ZN", "Subclone ZO", "Subclone ZP", "Subclone ZQ", "Subclone ZR", "Subclone ZS", "Subclone ZT", "Subclone ZU", "Subclone ZV", "Subclone ZX","Subclone ZY", "Subclone ZZ",
                "Subclone ZZA","Subclone ZZB","Subclone ZZC", "Subclone ZZD","Subclone ZZE", "Subclone ZZF", "Subclone ZZG", "Subclone ZZH", "Subclone ZZI", "Subclone ZZJ", "Subclone ZZK", "Subclone ZZL", "Subclone ZZM", "Subclone ZZN", "Subclone ZZO", "Subclone ZZP", "Subclone ZZQ", "Subclone ZZR", "Subclone ZZS", "Subclone ZZT", "Subclone ZZU", "Subclone ZZV", "Subclone ZZX","Subclone ZZY", "Subclone ZZZ",
                "Subclone ZZZA","Subclone ZZZB","Subclone ZZZC", "Subclone ZZZD","Subclone ZZZE", "Subclone ZZZF", "Subclone ZZZG", "Subclone ZZZH", "Subclone ZZZI", "Subclone ZZZJ", "Subclone ZZZK", "Subclone ZZZL", "Subclone ZZZM", "Subclone ZZZN", "Subclone ZZZO", "Subclone ZZZP", "Subclone ZZZQ", "Subclone ZZZR", "Subclone ZZZS", "Subclone ZZZT", "Subclone ZZZU", "Subclone ZZZV", "Subclone ZZZX","Subclone ZZZY", "Subclone ZZZZ")
  
  i <- 1
  s <- 1
  for(i in 1:(nrow(clonenames_new_order)/2)){ #Looping through all of the new subclones and giving them their new names.
    #print(i)
    if(clonenames_new_order[2*i-1,2] != "ALL"){
      clonenames_new_order[2*i-1,1] <- newnames[s]
      s <- s+1
    }else{
      clonenames_new_order[2*i-1,1] <- "Stem"
      clonenames_new_order[2*i-1,2:(ncol(samples)+1)] <- samples
      clonenames_new_order[2*i,1] <- "100"
      clonenames_new_order[2*i,2:(ncol(samples)+1)] <- sampleTC
      rowofall <- match("ALL",rownames(EM_newnames))
      rownames(EM_newnames)[rowofall] <- "Stem"
    }
    if(clonenames_new_order[2*i-1,2] == "Normal"){
      clonenames_new_order[2*i-1,1] <- "Normal"
    }
    i <- i+1
  }
  
  #Adding an "ALL" cell to the clonenames_new_order matrix in the cases where we only have 2 subclones. Otherwise we will not be able to construct any phylogenetic trees.
  if(as.numeric(nrow(clonenames_new_order))/2 < 2){ #Changed it to 2 rather than 3.
    ALL <- matrix(0,2,ncol(clonenames_new_order))
    ALL[1,1] <- "ALL"
    ALL[1,2:(ncol(samples)+1)] <- samples
    ALL[2,1] <- "100"
    ALL[2,2:(ncol(samples)+1)] <- sampleTC
    clonenames_new_order <- rbind(ALL,clonenames_new_order)
    print("Warning message: Your dataset only contain two subclones. An ALL subclone has been added in order to be able to reconstruct a phylogenetic tree. This has the same events as Stem.")
    
  }
  
  #Creating the new event matrix with the new subclones.
  i <- 1
  EM_saved <- EM_newnames
  for(i in 1:nrow(EM_newnames)){
    therow <- which(clonenames_new_order == rownames(EM_newnames)[i], arr.ind = T)
    therow <- therow[1]
    rownames(EM_newnames)[i] <- clonenames_new_order[therow,1]
    i <- i+1
  }
  
  EM_test_newnames <- EM_newnames
  #root <- "Stem"
  #root <- "Normal"
  if(root != "Stem"){
    EM_test_newnames <- stem(EM_test_newnames,stem_co,root) #Adding the root to the event matrix.
    allrow <- match("Stem",rownames(EM_test_newnames))
    
    if(is.na(allrow) == TRUE){
      allrow <- match("ALL",rownames(EM_test_newnames))
    }
    
    thesumofstem <- sum(EM_test_newnames[allrow,])
    kvoten <- thesumofstem/(ncol(EM_test_newnames))
    print(thesumofstem)
    print(ncol(EM_test_newnames))
    if(kvoten < 0.5){
      print("Warning message: The stem is very long compared to the entire data set. Maybe you should root the tree in the stem instead of a Normal cell?")
      print("Warning message: If the normal cell is too different from the subclones, the plot viewer might show a graph where the tip labels drift off from the tips")
    }
  }
  #else{EM_test_newnames <- stem(EM_test_newnames,stem_co,root) } #Adding the root to the event matrix.
  output <- list()
  output[[1]] <- EM_test_newnames
  output[[2]] <- clonenames_new_order
  return(output)
}

#Creating a distribution-plot.
distribution <- function(overview){
  i <- 2
  empty <- 0
  for(i in 2:nrow(overview)){
    name <- overview[i,1]
    j <- 2
    for(j in 2:ncol(overview)){
      biopsy <- overview[1,j]
      value <- as.numeric(overview[i,j])
      if(value != 0){
        df_el <- t(replicate(value,c(biopsy,name)))
        if(empty == 0){
          df <- df_el
          empty <- 1
        }else{
          df <- rbind(df,df_el)
        }
      }
      j <- j+1
    }
    i <- i+1
  }
  
  df <- as.data.frame(df)
  # Plot
  p <- ggplot(df, aes(y=V2, x=V1,  fill=V2,height = stat(count))) +
    geom_density_ridges(alpha=0.8, stat="binline",bins=(ncol(overview)-1),scale=0.8)+theme_ridges()+
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 45))+
    scale_fill_viridis_d(direction = -1, guide = "none")+
    xlab("") +
    ylab("")
  p
  #ggsave(p,filename= "Distribution.png",width = w,height = h)
  return(p)
}

#Creating pies.
make_pie <- function(clonenames_new_order, root, samples, type, custom_col){
  if(root == "Normal"){# && "Normal" %in% clonenames_new_order[1,1] == FALSE){
    Normal <- matrix(0,2,ncol(clonenames_new_order)) #Adding the normal cell to the clonenames_new_order matrix.
    Normal[1,1] <- "Normal"
    Normal[1,2] <- "100"
    Normal[2,1] <- "Normal cells"
    Normal[2,2] <- "100"
    clonenames_new_order <- rbind(Normal,clonenames_new_order)}
  
  Subclones <- matrix(0,2,100)
  pies <- list() #Creating a list for all pie data.
  pie_images <- list()
  
  pie_empty <- matrix(0,length(samples),2) #Creating empty pies.
  pie_empty[,1] <- samples
  pie_empty[,2] <- "0"
  
  i <- 1 #The following loop will extract the size of the subclone in each sample.
  for(i in 1:(nrow(clonenames_new_order)/2)){ #Looping through the new subclones.
    j <- 1
    s <- 1
    for(j in 1:ncol(clonenames_new_order)){ #Looping through the samples in which the subclone exists.
      
      if(clonenames_new_order[(2*i-1),j] != "0"){ #We should not add all of the columns with zeros.
        
        if(j != 1){ #We're not in the first column. The data includes the subclones within the new subclone.
          Subclones[1,s] <- word(clonenames_new_order[(2*i-1),j],1) #The sample.
          Subclones[2,s] <- clonenames_new_order[(2*i),j] #The size of the subclone within that sample.
          s <- s+1
        }else{Subclones[1,s] <- clonenames_new_order[(2*i-1),j] #We're in the first position. This is the new subclone name.
        Subclones[2,s] <- clonenames_new_order[(2*i),j] #The mean size of the subclone.
        s <- s+1
        }
        
      }
      
      j <- j+1
    }
    
    Subclones <- Subclones[,Subclones[1,] != "0"] #Removing the rows with zeros.
    if(Subclones[1,1]!="Normal"){
      Subclones <- distinct(data.frame(t(Subclones))) #Adding the vector to a list after removing rows that are equal. This contains all of the pie data needed.
    }
    pies[[i]] <- Subclones #Adding the vector to a list after removing rows that are equal. This contains all of the pie data needed.
    Subclones <- matrix(0,2,100) #Resetting the matrix.
    
    i <- i+1
  }
  
  image_names <- matrix(0,1,(nrow(clonenames_new_order)/2)) #Creating a vector that is to be used in order to save all of the file names.
  
  unique_biopsies <- unique(datasegment[,2]) #Unique biopsies.
  if(unique_biopsies[1]=="ALL"){ #Removing the ALL.
    unique_biopsies <- unique_biopsies[2:length(unique_biopsies)]}
  unique_biopsies <- c(c(unique_biopsies),c("Normal","Stem"))
  
  #This part should be looped for each matrix in "pies".
  j <- 1
  
  #Custom_colors
  blue <- c("#6bb5d8","#6fb9e7","#4d8dc6","#2c6f9a","#205575")
  red <- c("#ed6d70","#ea5456","#e52421","#a71916")
  yellow <- c("#f6c400","#ed8606","#e55514")
  #grey <- c("#b9b8b8","#9d9d9c","#706f6f","#3c3c3b") In article.
  green <- c("#add3a2","#6abfa4","#497f7a","#2c574a")
  brown <- c("#ca9e67","#936037")
  purple <- c("#d4bae0","#c9a0dc","#ae87d0","#7851a9","#522d80","#500691","#330066")
  grey <- c("#b9b8b8","#9d9d9c","#8a8a8a","#706f6f","#595858","#3c3c3b","#212121")
  
  #custom_col <- t(as.matrix(c(blue[2],red[2],yellow[2],green[2],grey[6],grey[6])))
  
  # custom_col <- t(as.matrix(c(blue[1],blue[3],blue[5],
  #                             green[1],green[2],green[3],
  #                             yellow[1],yellow[2],yellow[3],
  #                             red[2],red[3],red[4],
  #                             grey[1],grey[2],grey[3],grey[5],grey[6],grey[7],
  #                             purple[1],purple[2],purple[3],purple[4],purple[5],purple[6],purple[7],
  #                             red[1],red[1])))
  
  #PDX1  
  custom_col <- t(as.matrix(c(grey[2],grey[4], #Parental.
                              yellow[2],yellow[3],red[3], #Control.
                              blue[1],blue[3],blue[4],blue[5], #Treated.
                              grey[5],grey[6])))
  
  # #PDX2
  # custom_col <- t(as.matrix(c(yellow[2],yellow[3],red[3], #Control. (Parental har bara ALL!)
  #                             blue[3],blue[4],blue[5], #Treated.
  #                             grey[5],grey[6])))
  
  
  # custom_col <- t(as.matrix(c(green[1],blue[1],blue[3],blue[5],
  #                             yellow[1],yellow[2],yellow[3],grey[6],grey[6])))
  
  # custom_col <- t(as.matrix(c(yellow[1],yellow[2],yellow[3],red[4],  #PDX3
  #                            grey[6],grey[6])))
  
  # custom_col <- t(as.matrix(c(green[1],green[3],green[4], #M8
  #                             blue[1],blue[3],blue[4],blue[5],
  #                             yellow[1],yellow[2],
  #                             red[1],red[3],
  #                             grey[6],grey[6])))
  
  # custom_col <- t(c(green[1],green[3],blue[1],blue[3],blue[5],
  #                   yellow[1],yellow[2],yellow[3],red[4],grey[6],grey[6]))
  
  for(j in 1:length(pies)){ #Looping though all of the matrices in pies.
    Subclone <- pies[j] #Extracting the matrix j representing the data for a particular subclone that is to be presented in the shape of pie charts.
    Subclone <- as.matrix(as.data.frame(Subclone[1])) #Transforming it to a data frame.
    
    #The subclone vector.
    i <- 2
    for(i in 2:nrow(Subclone)){ #Looping through the samples that the subclone are present in.
      if(i == 2){
        y <- as.vector(rep(Subclone[i,1],2)) #Extracting the first sample name and creating a vector where it appears two times.
        a <- as.vector(rep(Subclone[i,2],2)) #Extracting the TC.
        a[2] <- (100 - as.numeric(a[1])) #The first position will be the TC and the other one 100-TC and hence the other slice of the pie chart.
        
      }else{
        x <- as.vector(rep(Subclone[i,1],2))
        b <- as.vector(rep(Subclone[i,2],2))
        b[2] <- (100 -as.numeric(b[1]))
        y <- c(y,x) #Combining the sample vectors.
        a <- c(a,b)} #Combining the TC values. This gives a vector with all the TC values that is to be used when dividing the pie charts.
      
      i <- i+1
    }
    
    #Creating the pie colors for this particular biopsy.
    sp <- brewer.pal(11,"Spectral")
    if(type=="col"){
      colors_biopsies <- as.matrix(cbind(c("B1",sp[2],"white"),c("B2",sp[10],"white"),c("B3",sp[9],"white"),c("B4",sp[5],"white"),c("B5",sp[1],"white"),c("B6",sp[11],"white"),c("B7",sp[8],"white"),c("B8",sp[3],"white"),c("B9",sp[6],"white"),c("B10",sp[4],"white"),c("B11",sp[7],"white"),c("B12","#008080","white"),c("B13","#800000","white"),c("B14","#808080","white")))
      #colors_biopsies <- as.matrix(cbind(c("B1","indianred1","white"),c("B2","#619CFF","white"),c("B3","#00BA38","white"),c("B4","#00BFC4","white"),c("B5","indianred1","white"),c("B6","#5fc400","white"),c("B7","#F564E3","white"),c("B8","#000485","white")))
      colors_biopsies[1,1:length(unique_biopsies)]<- unique_biopsies
    }else if(type=="nocol"){
      colors_biopsies <- as.matrix(cbind(c("B1","indianred1","white"),c("B2","indianred1","white"),c("B3","indianred1","white"),c("B4","indianred1","white"),c("B5","indianred1","white"),c("B6","indianred1","white"),c("B7","indianred1","white"),c("B8","indianred1","white"),c("B9","indianred1","white"),c("B10","indianred1","white"),c("B11","indianred1","white"),c("B12","indianred1","white"),c("B13","indianred1","white"),c("B14","indianred1","white"),c("B15","indianred1","white")))
      colors_biopsies[1,1:length(unique_biopsies)]<- unique_biopsies
    }else if(type=="custom"){
      names <- t(as.matrix(c(paste( c("B"), 1:as.numeric(length(unique_biopsies)), sep=""))))
      white <- t(as.matrix(c(rep("white",length(unique_biopsies)))))
      colors_biopsies <- rbind(names,custom_col,white)
      colors_biopsies[1,1:length(unique_biopsies)]<- unique_biopsies
    }else{
      print("You have not chosen a correct color mode.")
    }
    
    c <- 1
    for(c in 1:(length(y)/2)){
      column <- match(y[c+(c-1)],colors_biopsies[1,])
      if(c != 1){
        color_matrix <- c(c(color_matrix),c(colors_biopsies[2:3,column]))
      }else{
        if(Subclone[1,1] != "Normal"){
          color_matrix <- colors_biopsies[2:3,column]
        }else{
          color_matrix <- colors_biopsies[2:3,1]
        }
      }
      c <- c+1
    }
    
    test <- data.frame(Subclone = y,
                       Names = c(rep(c("Sample 1","Sample 2"),length(y)/2)),
                       TC = as.numeric(a),colour = color_matrix) #Creating a data frame with the samples in which the subclone exist, sample names, the TC for each and the colors.
    test$Subclone <- factor(test$Subclone, levels = unique(Subclone[2:nrow(Subclone),1]))
    
    x <- ggplot(test, aes(x="", y = TC, group = Names, fill = colour)) +
      geom_bar(width = 10, stat = "identity")+
      geom_col(position = "fill")+scale_fill_identity()+facet_grid(.~Subclone)+
      coord_polar("y", start=0) +theme_void()+#theme(strip.text.x = element_text(size = 200))+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank(), legend.position="none")+
      theme(plot.margin=unit(c(0,0,0,0),units = "lines"))
    
    if(type=="nocol" && Subclone[1,1] != "Normal"){
      x <- x+theme(strip.text.x = element_text(size = 200))
    }
    
    pie_images[[j]] <- x
    names(pie_images)[j] <- Subclone[1,1]
    pie_images[[j]]$Subclone <- Subclone[1,1]
    
    #Changes the element_text. Standard is 200 for all images. When you have more samples you might need to change it.
    w <- 49
    #w <- 10*(as.numeric(nrow(Subclone))-1)
    s <- 10
    ggsave(x,filename=paste(Subclone[1,1],".png",sep=""),width = w,height = s) #Testade att strunta i vidgningen av bilden.
    image_names[j] <- paste(Subclone[1,1],".png",sep="")
    
    if(Subclone[1,1]=="Stem"){
      
      x <- ggplot(test, aes(x="", y = TC, group = Names, fill = colour)) +
        geom_col(position = "fill")+facet_grid(.~Subclone)+#scale_fill_manual(values=color_matrix)+
        coord_polar("y", start=0) +theme_void()+labs("Samples")+
        scale_fill_identity(guide="legend",labels=c(samples),breaks = colors_biopsies[2,1:length(samples)],name="Samples")+
        guides(labels = guide_legend(override.aes = list(shape = 15)))+
        theme(plot.margin=unit(c(0,0,0,0),units = "lines"))
      plot(x)
      
      if(type != "nocol"){
        legend <- cowplot::get_legend(x)
        ggsave(legend,filename="legend.pdf",width=4,height=10,units = "cm")}
      
    }
    
    j <- j+1
  }
  pieData <- list()
  pieData[[1]] <- image_names
  pieData[[2]] <- pie_images
  return(pieData)
}

#Adding the pies.
pie_it <- function(Tree,pieData, offset, size,col){
  image_names <- pieData[[1]]
  pie_images <- pieData[[2]]
  p <- Tree
  labels <- as.matrix(p$data$label)
  labels <- as.matrix(labels[!is.na(labels)])
  positions <- matrix(0,nrow(labels),1) #Empty matrix in which the positions are to be saved.
  pie <- matrix(0,nrow(labels),1)
  
  #Extracting the subclone that each image belongs to.
  i <- 1
  s <- 1
  for(i in 1:length(image_names)){ #Looping through the image names.
    thesubclone <- word(image_names[i],1,sep = ".png") #Extracting the subclone that each image belongs to.
    thesubclone_pos <- match(thesubclone,labels)
    pie_image <- match(thesubclone,names(pie_images))
    positions[s,1] <- thesubclone_pos
    pie[s,1] <- pie_image
    s <- s+1
    i <- i+1
  }
  
  d <- data.frame(node = positions,images = c(image_names),pie_add = c(pie))
  
  #img <- readPNG("legend.png")
  if(missing(offset)==TRUE){
    offset <- 1
  }
  
  if(missing(size)==TRUE){
    size <- 0.21
  }
  
  #img <- readPNG("legend.png")
  #img <- image_read("legend.png")
  
  new <- p %<+% d + geom_tiplab(aes(image=images), geom="image",offset = offset, size = size)
  # if(col=="yes"){
  #   get_png <- function(filename) {
  #     grid::rasterGrob(png::readPNG(filename), interpolate = TRUE)
  #   }
  #   l <- get_png("legend.png")
  # new <- new+annotation_custom(l,xmin = max(p$data$x), xmax = max(p$data$x)+1.5, ymin = 0, ymax = 3)
  # }else if(col == "no"){
  #   print("No legend.")
  # }
  return(new)
}
####

######################
#Files to be analyzed#
######################
setwd("~/") #Set the working directory to where your file is.
data <- load_matrix(filename="Segment_example.xlsx",sheetname ="Example_tumors")
x <- "Tumor1"

#Rule matrix. The first object is the mother that the second one the daughter it cannot have according to
#information we have from some source.
# rule <- matrix(0,1,2)
# rule[1,1] <- "17p13q12 LOSS (1+0)"
# rule[1,2] <- "17p13q12 GAIN (2+1)"
# colnames(rule) <- c("Mother","Daughter")

#####################################
#Generating event matrices and trees#
#####################################
event_co <- 100000 #Cutoff for the start and end position.
root <- "Stem"
x <- "Tumor1" #Declare the name of the tumor in the data set you want to analyze.
datasegment <- splitdata(data,name=x) #Extracting the beginning and end position of each sample in the segment file. #Declare which tumor you want to analyze. Specified by the first column in your data set.

#Creating the event matrix.
EM <- DEVOLUTION(datasegment,event_co,datatypes=c("All"), eps = 0.5) #Creating an event matrix based on the segment file chosen.

#The final event matrix
EM_dev <- subclones(EM,file_samples_subclones,root = "ALL") #The first element in this list is the new event matrix. The second one is used for making pie charts.
DB <- distribution(overview_stem) # If you want to plot an overview.
plot(DB)
ggsave(DB,filename= "Distribution.png",width = 10,height = 10)

#Visualizing the trees without pies and saving them
EM_phy <- phydatevent(EM_dev[[1]]) #Transforming the EM to phyDat format.
EM_mptree <- mp_tree(EM_phy,root) #Constructing the maximum parsimony tree.
EM_mltree <- ml_tree(EM_phy,root) #Constructing the maximum likelihood tree.
limitmp <- xlim(c(0, 10)) #Here you can determine the limits for the graph for mp.
limitml <- xlim(c(0,1.2)) #Here you can determine the limits for the graph for ml.
Treemp <- MP_treeplot(EM_mptree,limitmp) #Illustrating the maximum parsimony tree.
Treeml <- ML_treeplot(EM_mltree,limitml) #Illustrating the maximum likelihood tree.

s = 10
ggsave(Treemp,filename=paste("Tree","_mp",".png",sep=""),width = s,height = s)
ggsave(Treeml,filename=paste("Tree","_ml",".png",sep=""),width = s,height = s)

#Saving the event matrix and the clustering annotation of each event.
#It is basically an updated version of the input segment file with the clustering of events.
write.xlsx(as.data.frame(t(EM_dev[[1]])),"DEVOLUTION.xlsx",sheetName="Event matrix")
write.xlsx(Clustering,append = TRUE, "DEVOLUTION.xlsx",sheetName = "Clustering") #Saving the data set that has been used in order to make the EM. It includes information about the subclonal belonging.
write.xlsx(as.data.frame(t(EM)),append = TRUE,"DEVOLUTION.xlsx",sheetName="Event matrix samples")
write.xlsx(as.data.frame(overview_stem),append = TRUE,"DEVOLUTION.xlsx",sheetName="Overview")

#Creating pie charts and saving the final tree.
coltype <- "col" #Choose how you want your pies. nocol = Just red pie charts with a biopsy name above. col = colored pies. custom = create your own color scheme.
samples <- as.vector(unique(datasegment[datasegment[,2]!="ALL",2])) #Or just write it.
pieData <- make_pie(EM_dev[[2]],root,samples,type=coltype) #Creates the pie charts.
pietree <- pie_it(Treemp,pieData,offset=1,size=0.17,col="yes") #Adds pie charts to the tree. 0.17 is preferred. Lower to 0.09 for large phylogenies.
ggsave(pietree,filename=paste(x,"Test_tree",".png",sep=""),width = s,height = s)
