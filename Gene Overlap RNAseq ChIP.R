#LiviuCengher
#CheungLab December 2018

#Look for overlap between ChIP and RNAseq Data for sarA
#This should be an adjustable alignment of ChIP and RNAseq data


library(tidyr)
library(dplyr)
library(stringr)

#Defining what lengh I'm looking at to start
Distance_sRNAs<- 100 #PosN #Use PosN instead of a defined integer for integrated Program

#Read in Files
sarA_Seq_2h <- read.csv("sarA_KO_HG003_Felden_RNASeq_2h.csv", 
                        stringsAsFactors=FALSE)
sarA_Seq_2h$FPKM.WT<- as.numeric(sarA_Seq_2h$FPKM.WT)
sarA_Seq_2h$FPKM.SarA<- as.numeric(sarA_Seq_2h$FPKM.SarA)
str(sarA_Seq_2h)
View(sarA_Seq_2h)

sarA_Seq_4h30min <- read.csv("sarA_KO_HG003_Felden_RNASeq_4h30min.csv", 
                             stringsAsFactors=FALSE)
sarA_Seq_4h30min$FPKM.WT<- as.numeric(sarA_Seq_4h30min$FPKM.WT)
sarA_Seq_4h30min$FPKM.SarA<- as.numeric(sarA_Seq_4h30min$FPKM.SarA)
str(sarA_Seq_4h30min)
View(sarA_Seq_4h30min)

Pull_DT <-read.csv("ChIP Pulldown Results MinO_1.csv")# I'm not using this now, but could be useful in the future. Is already in the finding sRNAs program

colnames(sarA_Seq_2h)
colnames(sarA_Seq_4h30min)
#ChIPsarA_AllsRNA <- rbind(Cleaned_sRNA_Output_Positive, Cleaned_sRNA_Output_Negative) # Use this one for integrated program
RNAseq_sarA_AllExpression <- merge(sarA_Seq_2h, sarA_Seq_4h30min,
                                   by = "ID")
  ## I broke this somehow, need to add back: , .x= "2h", .y= "4.5h")
#.x columns are 2h timepoint values, .y columns are 4h30min timepoint values 
str(RNAseq_sarA_AllExpression)
View(RNAseq_sarA_AllExpression)


########## This block condenses the sRNAs that I'm interested in. When I add this to the sRNA finding program it would be easier#######
positive_strand_sRNA<- read.csv("sRNA Output Positive 100bp (40).csv")
positive_strand_sRNA$SRD.identifier<- as.character(positive_strand_sRNA$SRD.identifier)
str(positive_strand_sRNA)
negative_strand_sRNA<- read.csv("sRNA Output Negative 100bp (42).csv")
negative_strand_sRNA$SRD.identifier<- as.character(negative_strand_sRNA$SRD.identifier)
str(negative_strand_sRNA)

View(positive_strand_sRNA)
View(negative_strand_sRNA)

ChIPsarA_AllsRNA <- rbind(positive_strand_sRNA, negative_strand_sRNA)
#ChIPsarA_AllsRNA <- rbind(Cleaned_sRNA_Output_Positive, Cleaned_sRNA_Output_Negative) # Use this one for integrated program

str(ChIPsarA_AllsRNA)
View(ChIPsarA_AllsRNA)

#####Overlap ChiP and RNAseq #####
# look at overlap between ChIPsarA_AllsRNA and RNAseq_sarA_AllExpression
#for all rows with RNAseq_sarA_AllExpression get rid of anything past 4 characters (8 scharacters total)
RNAseq_sarA_sRNAshortened_Expression <- RNAseq_sarA_AllExpression
RNAseq_sarA_sRNAshortened_Expression$ID <- substr(RNAseq_sarA_AllExpression$ID, 1, 8)
str(RNAseq_sarA_sRNAshortened_Expression)
View(RNAseq_sarA_sRNAshortened_Expression)

ChIPandSeq_intersect_sRNA <- as.data.frame(intersect(ChIPsarA_AllsRNA$SRD.identifier, 
                                                RNAseq_sarA_sRNAshortened_Expression$ID))
colnames(ChIPandSeq_intersect_sRNA)[1] <- "ID"

str(ChIPandSeq_intersect_sRNA)
View(ChIPandSeq_intersect_sRNA)

#ChIPandSeq_intersect_sRNA <- as.character(ChIPandSeq_intersect_sRNA$ID)
#str(ChIPandSeq_intersect_sRNA)
#View(ChIPandSeq_intersect_sRNA)

#Annotating the overlapping genes with the RNA-seq data
ChIPandSeq_intersect_sRNA_anno <- left_join(ChIPandSeq_intersect_sRNA, 
                                            RNAseq_sarA_sRNAshortened_Expression,
                                            by  = 'ID')
str(ChIPandSeq_intersect_sRNA_anno)
View(ChIPandSeq_intersect_sRNA_anno)

write.csv(ChIPandSeq_intersect_sRNA_anno, 
  "All sRNAs that are in RNAseq and have sarA under 100bp in 5' by ChIP.csv", 
  row.names = FALSE)

#####2 h Overlap ChiP and RNAseq #####

ChIPandSeq_intersect_sRNA_anno_significant2h <-ChIPandSeq_intersect_sRNA_anno[ChIPandSeq_intersect_sRNA_anno$p.value.adj.x <0.05, ]


write.csv(ChIPandSeq_intersect_sRNA_anno_significant2h, file = paste(
  "2h sRNAs that significant in RNAseq and have sarA under 100bp in 5' by ChIP.csv"), row.names = FALSE)

View(ChIPandSeq_intersect_sRNA_anno_significant2h)

#####4 h Overlap ChiP and RNAseq #####

ChIPandSeq_intersect_sRNA_anno_significant4h <-ChIPandSeq_intersect_sRNA_anno[ChIPandSeq_intersect_sRNA_anno$p.value.adj.y <0.05, ]


write.csv(ChIPandSeq_intersect_sRNA_anno_significant4h, file = paste(
  "4h sRNAs that significant in RNAseq and have sarA under 100bp in 5' by ChIP.csv"), row.names = FALSE)

View(ChIPandSeq_intersect_sRNA_anno_significant4h)

