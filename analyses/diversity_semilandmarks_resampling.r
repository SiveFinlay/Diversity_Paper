#23/03/15

#Sensitivity analysis in response to reviewer's comments

#Concern that the number of semilandmarks is over-sampling and introduces error in measurement of shape variation

#steps:
  #1) Read in a clean up raw landmark data  (select just tenrecs and golden moles)
  
  #2) Procrustes superimposition of tenrecs and golden moles
  
  #3) Randomly resample the number of semilandmarks
    
  #4) Find the average Procrustes shape coordinates for each species
  #5) PCA of the shape coordinates
  #6) Select PC axes that account for 95% of the variation

  #Repeat steps 3-6 for 25%, 50% and 75% of the total number of semilandmarks, 100 samples of each %

  #Output: table of the number of PC axes that account for 95% of the variation after every analysis
  
library(geomorph)
library(vegan)
library(plotrix)

source("C:/Users/sfinlay/Desktop/Diversity/functions/Morpho_diversity_functions.r")

setwd("C:/Users/sfinlay/Desktop/Diversity/data/")

#######################
#Read in the data: 3 views of skulls
#################
#SkDors
    #1) Landmarks
      #land <- readland.tps(file="skdors/Skdors_16_12_13_10landmarks+4curves_edited.TPS")
    #2) Sliders
      #curves <- as.matrix(read.table("skdors/Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))
    #3) Taxonomy
      #taxa <- read.csv ("skdors/Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)
    #4) Specimens to remove
      #Null

#-----------------------------------------------------
#SkVent
  #1) Landmarks
    #land <- readland.tps(file="skvent/SkVent_30_10_13_13landmarks+1curve_edited.TPS")
  #2) Sliders
    #curves <- as.matrix(read.table(file="skvent/SkVent_1skull_13landmarks+1curve_sliders_edited.tps", header=TRUE))     #this is a tps file and the others are nts but it doesn't make a difference
  #3) Taxonomy
    #taxa <- read.csv("skvent/SkVent_30_10_13_imagelist+specimens.csv" , header=TRUE)
  #4) Specimens to remove
    #rem <- read.csv("skvent/SkVent_remove_spec.csv", header=T)

#--------------------------------------------------------
#SkLat
  #1) Landmarks
    land <- readland.tps(file="sklat/SkLat_08_11_13_9landmarks_2curves_edited.TPS")
  #2) Sliders
    curves <- as.matrix(read.table(file="sklat/SkLat_08_11_13_9landmarks_2curves_sliders_edited.NTS", header=TRUE))
  #3) Taxonomy
    taxa <- read.csv("sklat/SkLat_08_11_13_Specimens+images.csv", header=TRUE)
  #4) Specimens to remove
    rem <- read.csv("sklat/SkLat_remove_spec.csv", header=T)


#################################################
#CLEAN UP THE DATA
#################################################
#Combine the taxonomic information with the array of landmark data
  combine <- list(land=land, curves=curves, ID=taxa$ID,SpecID=taxa$SpecID, Order=taxa$Order_05,
                  Fam=taxa$Family_05, Genus=taxa$Genus_05, Species=taxa$Species_05, Binom=taxa$Binomial_05)

# Remove the _sp. specimens
 sp <- which(combine$Species=="sp.")

 combine <- remove.from.list(combine, sp)
 combine <- droplevels.from.list(combine)

#********************************************
#OPTION; depending on the data and the analysis
#************************************
#Remove the specimens listed in rem (sklat, and skvent data)
  #these are photographs of damaged specimens which could not be included in the landmark analyses
  #doesn't apply to the skdors data because rem is NULL

#find the ID numbers of specimens with missing data
  matching <- matching.id(rem$SpecID, combine$SpecID)
    combine <- remove.from.list(combine, matching)
    combine <- droplevels.from.list(combine)

##################################
#Select specimens that you want
#####################################
#Select the tenrecs and golden moles only  (full data sets include other small mammal species)
  tc.gm <- c(which(combine$Fam=="Chrysochloridae"), which(combine$Fam=="Tenrecidae"))

  mydata <- select.from.list(combine, tc.gm)
  mydata <- droplevels.from.list(mydata)
  
##############################################
#Procrustes superimposition 
##########################################
#General Procrustes Alignment of all of the scaled coordinates
  mydataGPA <- gpagen(mydata$land, curves=mydata$curves, ProcD=TRUE,)
  #ProcD=TRUE means that the coordinates are aligned by procrustes distance rather than bending energy
      # which means that RWA is equivalent to PCA (Zelditch 2012, page 150)

#List the coordinates with the taxonomic information
  Proc.co <- list(coords=mydataGPA$coords,csize=mydataGPA$Csize,ID=mydata$ID,SpecID=mydata$SpecID,
                  Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

#################################################
#Random re-sampling of Procrustes-superimposed semilandmarks
#################################################

#Five steps in a wrapper function
  #Slow function because it includes running a PCA
  #1) Select a random sample of rows in a matrix (random selection of the semi landmark rows
  #2) Combine landmark data with subsample of semilandmark data (~NB centroids issue)
  #3) Species averaging
      # group the arrays of coordinates according to species
      #average coordinate values for each species
  #4) PCA
  #5) Find the number of PC axes that account for 95% of the variation

#------------------------------------------------------------
#Select which analyses to run

#Skulls Dorsal (10 landmarks, 44 semilandmarks)
  #sk.landmarks <- 10
  #sk.semilandmarks <- 44
  
#Skulls Ventral (13 landmarks, 60 semilandmarks)
  #sk.landmarks <- 13
  #sk.semilandmarks <- 60
  
#Skulls Lateral (9 landmarks, 35 semilandmarks)
  sk.landmarks <- 9
  sk.semilandmarks <- 35

#---------------------------------------------------------------

#Run the function for the same percentage multiple times
samp90 <- NULL
  for (i in 1:100){
    samp90[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=9/10)
    }
    
samp80 <- NULL
  for (i in 1:100){
    samp80[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=4/5)
    }
    
samp75 <- NULL
  for (i in 1:100){
    samp75[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=3/4)
    }
    
samp70 <- NULL
  for (i in 1:100){
    samp70[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=7/10)
    }
    
samp65 <- NULL
  for (i in 1:100){
    samp65[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=65/100)
    }
    
samp60 <- NULL
  for (i in 1:100){
    samp60[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=3/5)
    }
    
samp50 <- NULL
  for (i in 1:100){
    samp50[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=1/2)
    }
    
samp40 <- NULL
  for (i in 1:100){
    samp40[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=2/5)
    }
    
samp30 <- NULL
  for (i in 1:100){
    samp30[i] <- semilandmark.subsample.PC95axes(landmarks=sk.landmarks, semilandmarks=sk.semilandmarks, percentage=3/10)
    }
#---------------------------------------------
#Summary output table       (plotting a summary is awkward because the number of axes are binary results)
  PC.table <- matrix(nrow=4, ncol=9)
  rownames(PC.table) <- c("7PCaxes","6PCaxes","5PCaxes","4PCaxes")
  colnames(PC.table) <- c("samp90", "samp80", "samp75","samp70","samp65","samp60","samp50","samp40","samp30")
  
  PC.table[,1] <- c(length(which(samp90==7)),length(which(samp90==6)), length(which(samp90==5)), length(which(samp90==4)))
  PC.table[,2] <- c(length(which(samp80==7)),length(which(samp80==6)), length(which(samp80==5)), length(which(samp80==4)))
  PC.table[,3] <- c(length(which(samp75==7)),length(which(samp75==6)), length(which(samp75==5)), length(which(samp75==4)))
  PC.table[,4] <- c(length(which(samp70==7)),length(which(samp70==6)), length(which(samp70==5)), length(which(samp70==4)))
  PC.table[,5] <- c(length(which(samp65==7)),length(which(samp65==6)), length(which(samp65==5)), length(which(samp65==4)))
  PC.table[,6] <- c(length(which(samp60==7)),length(which(samp60==6)), length(which(samp60==5)), length(which(samp60==4)))
  PC.table[,7] <- c(length(which(samp50==7)),length(which(samp50==6)), length(which(samp50==5)), length(which(samp50==4)))
  PC.table[,8] <- c(length(which(samp40==7)),length(which(samp40==6)), length(which(samp40==5)), length(which(samp40==4)))
  PC.table[,9] <- c(length(which(samp30==7)),length(which(samp30==6)), length(which(samp30==5)), length(which(samp30==4)))

##########################
#Save the summary table


setwd("C:/Users/sfinlay/Desktop/Diversity/output/semilandmark_resampling")

#Skdors
    #capture.output(PC.table, file= "skdors_PCresamp.tab.txt")
    
#Skvent
    #capture.output(PC.table, file= "skvent_PCresamp.tab.txt")

#Sklat
    #capture.output(PC.table, file= "sklat_PCresamp.tab.txt")