C:\Users\Liviu.LAPTOP-K0MK2U5Q\Desktop\sRNA Overlap ChIP and RNAseq (sarA_2)
Finding_sRNA_nearPeaks.R
	Input
		"ChIP Pulldown Results MinO_1.csv"
		"sRNAs_NTCT8325.csv"
	Varing Variables
		PosN <- as.numeric(100) #
	Output
		"sRNA Output Negative..."	#output list final
		"sRNA Output Positive..."	#output list final

Gene Overlap RNAseq ChIP.R
	Input
		"sarA_KO_HG003_Felden_RNASeq_2h.csv"		 #the two RNA-seq list were R bound together
		"sarA_KO_HG003_Felden_RNASeq_4h30min.csv"
		"ChIP Pulldown Results MinO_1.csv"

		"sRNA Output Positive 100bp (40).csv" # R bound together
		"sRNA Output Negative 100bp (42).csv"
	Output
***		"All sRNAs that are in RNAseq and have sarA under 100bp in 5' by ChIP.csv"		#looks at all sRNAs
*		"2h sRNAs that significant in RNAseq and have sarA under 100bp in 5' by ChIP.csv"	#looks at only the ones relevant in the 2h sRNAs
*		"4h sRNAs that significant in RNAseq and have sarA under 100bp in 5' by ChIP.csv"	#looks at only the ones relevant in the 4h sRNAs4


Summary
1) sRNA_nearPeaks
	Finds sRNAs that are near peaks
	This is already done by CLC workbench automatically for mRNAs
2) Gene Overlap RNAseq ChIP
	Only sRNAs
	Final output for sRNAs
