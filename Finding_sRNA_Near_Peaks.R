#LiviuCengher
#CheungLab December 2018

#For every sRNA, looks for a close sarA ChIP peak
#A) upstream of + sRNAs ~500bp: Positive strand start -500
#B) downstream of - sRNAs ~500bp: Negative strand end +500 bp

library(tidyr)
library(dplyr)
library(stringr)


###### Looking at sRNAs ######## 

## Distance from 5' of sRNA to peaks as PosN
PosN <- as.numeric(100) 

## Load Table
sRNA_8325 <-read.csv("sRNAs_NTCT8325.csv")
Pull_DT <- read.csv("ChIP Pulldown Results MinO_1.csv")

## Check Structure and view tables
str(sRNA_8325)
str(Pull_DT)
View(sRNA_8325)
View(Pull_DT)

## Turn integer into numeric (correct formatting)
Pull_DT$CenterOfPeak <-as.numeric(Pull_DT$CenterOfPeak)
sRNA_8325$Start <- as.numeric(sRNA_8325$Start)
sRNA_8325$End <- as.numeric(sRNA_8325$End)

## Seperate sRNAs into + and - Strand sRNAs   
sRNA_positive <- as.data.frame(subset(sRNA_8325, Strand =="+"))
View(sRNA_Output_Positive) #Look for positive control (teg88 with 1 peak before it, 5')
sRNA_positive[53,7]

sRNA_negative <- as.data.frame(subset(sRNA_8325, Strand =="-"))
View(sRNA_Output_negative)#Look For negative control (teg49 with 1 peak after it, 5')
sRNA_negative[95,7]

## For both Positive and negative sRNAs I pull out relevant info, enter into new data frame, followed by peak info

######## For Positive strand sRNAs ######## 

sRNA_Output_Positive <- apply(sRNA_positive, 1 , function(x){ #1 = by row

  PeakSearch <- (as.numeric(x[3])-PosN) #PeakSearch defines the start of the search for each sRNA (depending on posN)
  

  PeaksCounted <- sapply(Pull_DT$CenterOfPeak, function(y){
  if(PeakSearch < y & y < as.numeric(x[3])){y} # Looks at peaks upstream of the Peaksearch value and downstream of the start site
  })

c(x[1], x[2], x[3], x[5], x[7], x[8], x[9], paste(PeaksCounted, collapse= ",")) # Compiles the information I want and collapses all peak results into a single box

} )


#View(t_sRNA_Output_Positive) #UGH this is a mess. Some peak results are "Null" instead of "", and it is oriented the wrong way 

## Transpose for correct orientation
t_sRNA_Output_Positive <- as.data.frame(t(sRNA_Output_Positive)) # Changes table to something easier to work with
colnames(t_sRNA_Output_Positive)[8]<- "CenterOfPeak"

str(t_sRNA_Output_Positive)

t_sRNA_Output_Positive$CenterOfPeak <- str_replace_all(t_sRNA_Output_Positive$CenterOfPeak, fixed("NULL"), " ") #Removes Null values in in CenterOfPeak
t_sRNA_Output_Positive$CenterOfPeak <- str_replace_all(t_sRNA_Output_Positive$CenterOfPeak, fixed(" ,"), "")#Removes commas that had null to the left in in CenterOfPeak
t_sRNA_Output_Positive$CenterOfPeak <- str_replace_all(t_sRNA_Output_Positive$CenterOfPeak, fixed(", "), "")#Removes commas that had null to the right in in CenterOfPeak
t_sRNA_Output_Positive$CenterOfPeak <- str_replace_all(t_sRNA_Output_Positive$CenterOfPeak, fixed(" "), "")#Removes all spaces in CenterOfPeak

#View(t_sRNA_Output_Positive)
str(t_sRNA_Output_Positive)


Cleaned_sRNA_Output_Positive <- t_sRNA_Output_Positive[!t_sRNA_Output_Positive$CenterOfPeak == "",]
#View(Cleaned_sRNA_Output_Positive)

##Writes file, lists search area and number of hits found
write.csv(Cleaned_sRNA_Output_Positive, file = paste(
  "sRNA Output Positive ",as.character(PosN) , "bp (",
  nrow(Cleaned_sRNA_Output_Positive), ").csv", sep = "", collapse = NULL), row.names = FALSE)

######## For NEGATIVE strand sRNAs ######## 
sRNA_Output_Negative <- apply(sRNA_negative, 1 , function(x){ #1 = by row

  # Pull out whole row, enter into new data frame, followed by sRNA name/location/etc used
  ## Script: 
  PeakSearch <- (as.numeric(x[4]) + PosN) #PeakSearch defines the start of the search for each sRNA (depending on posN)
  
  PeaksCounted <- sapply(Pull_DT$CenterOfPeak, function(y){
    if(PeakSearch > y & y > as.numeric(x[4])){y}# Looks at peaks upstream of the end site value and downstream of the peaksearch
      })

  c(x[1], x[2], x[3], x[5], x[7], x[8], x[9], paste(PeaksCounted, collapse= " "))# Compiles the information I want and collapses all peak results into a single box
} )

## Transpose for correct orientation, Clean up strings
t_sRNA_Output_Negative <- as.data.frame(t(sRNA_Output_Negative))
colnames(t_sRNA_Output_Negative)[8]<- "CenterOfPeak"

View(test_t_sRNA_Output_Negative)
write.csv(test_t_sRNA_Output_Negative, "(200313) test_t_sRNA_Output_Negative.csv")
str(t_sRNA_Output_Negative)

t_sRNA_Output_Negative$CenterOfPeak <- str_replace_all(t_sRNA_Output_Negative$CenterOfPeak, fixed("NULL"), " ") #same as for positive
t_sRNA_Output_Negative$CenterOfPeak <- str_replace_all(t_sRNA_Output_Negative$CenterOfPeak, fixed(" "), "")#same as for positive

View(t_sRNA_Output_Negative)
str(t_sRNA_Output_Negative)

Cleaned_sRNA_Output_Negative <- t_sRNA_Output_Negative[!t_sRNA_Output_Negative$CenterOfPeak == "",]
str(Cleaned_sRNA_Output_Negative)
View(Cleaned_sRNA_Output_Negative)

##Writes file, lists search area and number of hits found
write.csv(Cleaned_sRNA_Output_Negative, file = paste(
  "sRNA Output Negative ", as.character(PosN), "bp (",
  nrow(Cleaned_sRNA_Output_Negative), ").csv", sep = "", collapse = NULL), row.names = FALSE)


######## Merging sRNA Positive and Negative strand list ######## 

# In file "Gene Overlap RNAseq ChIP.R"
