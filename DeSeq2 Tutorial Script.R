##Next Steps
## - Deseq2 package- learn + practice
##Bulk RNA seq data is better - Find data sets
##Correlate GWAS findings to Literature Review - note to self: Stay focused on XChrom
##Data processing: FASTQ > align to ref. genome (hg19) > SAM > BAM > VCF
##If we can download BAM files that is great
##FASTQ will have to be run through Galaxy 
BiocManager::install("DESeq2")
browseVignettes("DESeq2")

##DESeq2 Tutorial
#install.packages("htmltools")
#library(htmltools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

install.packages("htmltools")
library(htmltools)
source("https://bioconductor.org/biocLite.R")
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)

##Tutorial: https://bioconnector.github.io/workshops/r-rnaseq-airway.html#deseq2_analysis

##I had to go download the CSV files for this directly because the HTML code wasn't working.

##Tutorial From https://bioconnector.github.io/workshops/r-rnaseq-airway.html#deseq2_analysis
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)

##Using sample data first.
file.choose()
countData <- read_csv("C:\\Users\\jgrib\\Desktop\\Netgene CoLab\\airway_scaledcounts.csv")
metaData <-  read_csv("C:\\Users\\jgrib\\Desktop\\Netgene CoLab\\airway_metadata.csv")
head(countData)
head(metaData)

mycounts <- as.data.frame(mycounts)
metadata <- as.data.frame(metadata)
head(mycounts)
head(metadata)
class(mycounts)
class(metadata)

##This checks that the IDs for "mycounts" are aligned to the correct metadata samples (sorted by SRR)
names(mycounts)[-1]
metadata$id
names(mycounts)[-1]==metadata$id
all(names(mycounts)[-1]==metadata$id)

dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)
dds
dds <- DESeq(dds)

res <- results(dds, tidy=TRUE)
res <- as_tibble(res)
res

res %>%
  arrange(padj)

##For our project, developing a piping technique for isolation of X Chromosome genes may be necessary
##example: res %>%, filter(X Chrom)
##Will need a complete cataloge of X chromosome genes to do this most likely

plotCounts(dds, gene="ENSG00000103196", intgroup="dex")

plotCounts(dds, gene="ENSG00000103196", intgroup="dex", returnData=TRUE) %>% 
  ggplot(aes(dex, count)) + geom_boxplot(aes(fill=dex)) + scale_y_log10() + ggtitle("CRISPLD2")

##Cool. Note that these are Ensembl gene codes, not the standard gene codes used in Seurat.
##For example, ENSG00000002586 = CD99. See below: 

plotCounts(dds, gene="ENSG00000002586", intgroup="dex")
plotCounts(dds, gene="ENSG00000103196", intgroup="dex", returnData=TRUE) %>% 
  ggplot(aes(dex, count)) + geom_boxplot(aes(fill=dex)) + scale_y_log10() + ggtitle("CD99")

##There may be a way to take all these Ensembl codes for the X Chromosome and filter by them. 
##I converted CD99 to Ensembl by looking it up on GeneCards. There is probably an easier way. 
##From here, there are several visualization options. 

##Application: GSE179427
file.choose()
countDataMS <- read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\GSE179427_countmtx.csv", sep=";")
##Make sure your .csv is actually a .csv!
countDataMS <-read_delim(("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\GSE179427_countmtx.csv"),delim = ';')
head(countDataMS)
countDataMS
metaDataMS <-  read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\SraRunTable.txt")
head(metaDataMS)
metaDataMS

##The metadata is good, but I don't know how to align the count data with it? Note that this metadata includes scRNA seq and bulk RNA seq runs.
##I realized here that I messed up my metadata- the original file contained both bulk and
##scRNAseq data. Instead of coding it out I just selected by sequencing type on SRA
##And downloaded a bulk-only file. "C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\SraRunTableBulkOnlyGSE179427.txt"

##Now head(metadataMS) below wills show only bulk RNAseq samples. This is what I want. 

metaDataMS <-  read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\SraRunTableBulkOnlyGSE179427.txt")
head(metaDataMS)

##That should fix that. Taken from SRA run selector, downloaded Metadata sorted by RNA seq type. 
##Now, can I allign the SRR runs with the names that were given?
##Metadata really have all duplicate values, excepting number of bases. DO I CARE ABOUT NUMBER OF BASES? - (for pheno data)

##First, let's get these in the proper format for Deseq2
countdataMS <- as.data.frame(countDataMS)
metadataMS <- as.data.frame(metaDataMS)
head(countdataMS)
head(metadataMS)

##For Deseq2 to work, both need to be class = data.frame. Check with this: 
class(countdataMS)
class(metadataMS)

##This is where things are going to get tricky. Names aren't the same, and I need to allign them
names(countdataMS)[-1]
metadataMS$id
##This should have the same readout. It doesn't. ******
##This means that my colData does not match my countData. I knew this was going to happen.
##To fix, I need to find a way to combine the 10+ SRR runs into the sample names provided. For example,
##SD042_13 corresponds to one patient, which has 8 attached SRR Runs: 
##SRR15038127, SRR15038128, SRR15038129, SRR1503830, SRR1503831, SRR1503832, SRR1503833, SRR1503834
##From what I can tell, the metadata I need should be identical for all 8 of these samples.
##Is there a way to consolidate all of them as "SD042_13" in R?
##Alternatively, I could make another .csv manually from the metadata on SRA and the patient IDs
##That would be tedious, but not impossible (~68 patients). But I'm sure there's a coding way to do this.
names(countdataMS)[-1]==metadataMS$id
all(names(countdataMS)[-1]==metadataMS$id)

##Trying to manipulate 1 run. These are all the SRR metadata files for patient ID S08_047, male, MS lesion group
metadata1 <- metadataMS %>%
  filter(Run == c("SRR15038303","SRR15038304","SRR15038305","SRR15038307"))
metadata1

metadata1 <- metadataMS %>%
  filter(Run == c("SRR15038303"))
metadata1
##For some reason, this falls apart when I add more than 5 samples to the above list?

##Filtering Metadata
{metadata2 <- metadataMS %>%
  filter(Age == 49)
metadata2
##I don't understand why this isn't working. It just so happens all 8 runs should be included
##in this filter (I manually noted the patient's age was 42)
##Nevermind he was 49 that's why it didn't work. Works now. 
##This is a way of separating out samples, but I don't think its the best way (2 patients could have the same age)
##Now that I'm looking at it- sorting by sample makes a lot more sense (GSM5417528)
##I am going to look at this sample set manually.
}

##Metadata and run data from one patient only
{
file.choose()
metadataONE <-  read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\MetadataONEONLY.txt")
head(metadataONE)
metadataONE <- as.data.frame(metadataONE)
head(metadataONE)
##now I manually change the name from SRR15038303 to S08_047 on the notepad
metadataONE <-  read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\MetadataONEONLY.txt")
head(metadataONE)
metadataONE <- as.data.frame(metadataONE)
head(metadataONE)
##"Run" is now S08_047. 

countdataONE <- countdataMS %>%
  select(S08_047)
countdataONE
##This pulls out all of the reads for S08_047. But I need the Ensembl column as well
##Labeled as "...1"
countdataONE <- countdataMS %>%
  select(...1,S08_047)
countdataONE
##And there we go. I have metadata and count data for a single sample. Let's run it through the workflow. 

##Ensure both are data.frame
class(countdataONE)
class(metadataONE)

##Alligned names
names(countdataONE)[-1]
metadataONE$Run
names(countdataONE)[-1]==metadataONE$Run
all(names(countdataONE)[-1]==metadataONE$Run)
##Make sure you have this aligned to "Run" if you do it the way I did it manually, by reassigning the Run name to the Sample Id. If you put in "id" it won't find it. 
##True, True, all checks out. 

dds <- DESeqDataSetFromMatrix(countData=countdataONE, 
                              colData=metadataONE, 
                              design=~dex, 
                              tidy=TRUE)
dds
##We may have another problem here- Can't do differential analysis on a single patient (Duh)
##Need to add in another (control)
##Now we have 3 patients- SO42_13 (control, female), SO8_047 (lesion, male), and S96_234 (lesion, female)

}

##Metatada and Run data from 3 patients (SD042_13,S08_047,S96_234)
##I manually made the "metadataTHREE" file by replacing SRR runs with patient tags in text editor. 
metadataTHREE <-  read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\MetadataTHREEONLY.txt")
head(metadataTHREE)
metadataTHREE <- as.data.frame(metadataTHREE)
head(metadataTHREE)
metadataTHREE

countdataTHREE <- countdataMS %>%
  select(...1,SD042_13,S08_047,S96_234)
countdataTHREE

names(countdataTHREE)[-1]
metadataTHREE$Run
names(countdataTHREE)[-1]==metadataTHREE$Run
all(names(countdataTHREE)[-1]==metadataTHREE$Run)
##Readout: True True True True we're good to go

colnames(countdataTHREE) <- c("Gene", "SD042_13","S08_047","S96_234")

##I accidentally removed the names from my columns with this, but realized I could use it to add them back 

ddsTHREE <- DESeqDataSetFromMatrix(countData=countdataTHREE, 
                              colData=metadataTHREE, 
                              design=~Group, 
                              tidy=TRUE,
                              ignoreRank = FALSE)
ddsTHREE

##Rabbit hole with naming Lesion group, etc. Did not affect workflow (minimized)
{
##So we have the right stuff, but its in the wrong format. Error:
##Error in DESeqDataSet(se, design = design, ignoreRank) : all variables in design formula must be columns in colData
##I'm sure there's a way to mutate rows into columns...if that is the problem?
##Not sure what is going on here.
##I think I got it- "dex" was the column from the vingette. I don't want that one here.
##Specified "Group" as the design (control vs. lesion)

##I have a hunch that the results will be easier to understand if I use factoextra and code 
##control and lesion groups as 0 and 1, respectively. Let's do that. 

##Factoextra conversion of group names to integers: 
{
library(factoextra)
metadataTHREE$Group <- ifelse(metadataTHREE$Group == "Lesion", 1, 0)
metadataTHREE

##Nice. Now control group is 0 and Lesion group is 1. Running the same code as before
dds <- DESeqDataSetFromMatrix(countData=countdataTHREE, 
                              colData=metadataTHREE, 
                              design=~Group, 
                              tidy=TRUE,
                              ignoreRank = FALSE)
dds
 ##This gives a new readout- the design formula contains one or more numeric variables with integer values,
##specifying a model with increasing fold change for higher values.
##did you mean for this to be a factor? if so, first convert
##this variable to a factor using the factor() function
##I don't know exaclty what a factor function does. I'm rabit holeing this for now
##Because the first time around it auto-converted
}
##ddsTHREE is still from the OG code before this, sticking with it for now. 

}

##Deseq
ddsTHREE <- DESeq(ddsTHREE)

##Continuing the workflow
res <- results(ddsTHREE, tidy=TRUE)
res <- tbl_df(res)
res

res %>% 
  filter(padj<0.05) %>% 
  write_csv("metadataTHREEsig.csv")
##And it is saved to the X Chrom Project folder on my computer. 

##Now for plotting counts
plotCounts(ddsTHREE, gene="ENSG00000102245", intgroup="Group", returnData = TRUE)
plotCounts(ddsTHREE, gene="ENSG00000102245", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("CD40LG")


plotCounts(ddsTHREE, gene="ENSG00000147050", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("KDM6A")

##Plugged in the Ensembl codes for CD40LG and KDM6A - X Chromosome candidate genes. 

##I reran the code to get rid of the factoextra step, since that was throwing things off
##lesson learned- let Deseq make factors for itself, don't try to do that on your own
##This gives a pretty wimpy boxplot but has both control and lesion groups (only 3 patients here; 2 lesion, 1 control)

# Create the new column
res <- res %>% mutate(sig=padj<0.05)

# How many of each?
res %>% 
  group_by(sig) %>% 
  summarize(n=n())

plotMA(ddsTHREE)
plotMA(res)
##I'm not sure I'm doing the MA plot correctly. Should it be the Deseq object or the results df? 
res$sig <- res$padj < 0.05
plotMA(res[,c("baseMean","log2FoldChange","sig")], colLine = "tan2", colSig = "red", colNonSig = "midnightblue", xlim=c(200,10000))
##This fixes that. It needs to run on 3 columns, 2 are baseMean and log2FoldChange, the other is the logical "sig" (T/F/NA), which I created above as described by the padj column
##Below was from the markdown code from the tutorial (I apparently did this the hard way)
res %>% ggplot(aes(baseMean, log2FoldChange, col=sig)) + geom_point() + scale_x_log10() + ggtitle("MA plot")

##Volcano
res %>% ggplot(aes(log2FoldChange, -1*log10(pvalue), col=sig)) + geom_point() + ggtitle("Volcano plot")

##The next step in the workflow is transformation. 

vsdataTHREE <- vst(ddsTHREE, blind=FALSE)

plotPCA(vsdata, intgroup="Group")

##Lets try with log transformation
rld <- rlogTransformation(ddsTHREE)
plotPCA(rld, intgroup="Group")

sessionInfo()
##That concludes the tutorial pretty much. Now the challenge is getting all the metadata
##Lined up with all the patient samples

##For Heatmaps:
install.packages("NMF")
library(NMF)

aheatmap(assay(rld)[arrange(res, padj, pvalue)$row[1:50], ], 
         labRow=arrange(res, padj, pvalue)$symbol[1:50], scale="row", distfun="pearson", annCol=dplyr::select(metadataTHREE, Group, gender), 
         col=c("green3","gray90","gray90","red2"))
##It's a mess but at least the steps in the process are clear now. Remember, this is only 3 patients. 
##Kind of cool to see how there are differences between the control and lesion patients from a male/female split, even in limited view. 
##The real question: are there differences between male and female patients across SLE/MS? 
##Is there a good way to extract? I want the Ensembl codes so I can look them up.

##Lets get to the big dataset now. We have the sequencing reads, but need the metadata alligned. 

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(NMF)
library(factoextra)

##for presentation
metadataALL%>%
  select(Group) -> groups1
groups1 %>%
  filter(Group=="Control") -> groups1control
count(groups1control)

file.choose()
metadataALL <-  read_csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\MetaDataALL_GSE179427.txt")
head(metadataALL)
metadataALL <- as.data.frame(metadataALL)
head(metadataALL)
rownames(metadataALL)
metadataALL==
head(countDataMS)
colnames(countDataMS)

names.conv <- read.csv("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\Manual Pt_Tag.SRR.Expt.csv")
names.conv
##What I would like to do is use the above sheet, and convert all of the SRRs to patient IDs
##Or vice versa, it does not matter. It needs to be a logical function
##Where if "SRR1538123", change name to "SD036_12" and so on. 

##All SRRs
namesSRR <- names.conv %>%
  select(SRR) 
namesSRR

(namesSRR)

##All Patient IDs
namesID <- names.conv %>%
  select(Patient.ID)
namesID

metadataALL %>%
  mutate(Run=recode(Run, 
                      SRR15038127 = 
                      "SD042_13"))

##This works for one sample at a time. Is there a way I can do this for all of them at once?
##Also, as a sidenote, is there a way to "UNDO" something in R? Like if I wanted to Ctrl+Z the last step?

##Lets try this

metadataALL %>%
  mutate(Run=recode(Run,
                    (namesSRR = namesID)))

##Not working. I am definitely missing a step. 

########################START HERE

##6.28.22
##From James- Join operations w/dplyr - https://dplyr.tidyverse.org/reference/mutate-joins.html

library(readr)
library(dplyr)
file.choose()
setwd(dir = "C:\\Users\\jgrib\\Desktop\\X Chromosome Project")
names.conv <- read_csv("Manual Pt_Tag.SRR.Expt.csv")
metadataALL <- read_csv("MetaDataALL_GSE179427.txt")
names.conv
metadataALL

metaDataMS <- left_join(metadataALL,
                           names.conv,
                           by =c("Run"="SRR"),
                           suffix = c("_meta","_names"),
                           keep = FALSE,
                           na_matches = "never")
metaDataMS

##This works, now I have Patient.ID as a column on metadata
##We should be able to run the Deseq2 workflow with this metadata. I may need to reorganize
##So that "Patient.ID" is the first column. 


file.choose()
countDataMS <-read_delim(("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\GSE179427_countmtx.csv"),delim = ';')
head(countDataMS)
countDataMS
head(metaDataMS)
metaDataMS

countdataMS <- as.data.frame(countDataMS)
metadataMS <- as.data.frame(metaDataMS)
head(countdataMS)
head(metadataMS)
metadataMS$`Patient ID`
rename(metadataMS, "Patient.ID" = "Patient ID")

##Moving Patient.ID to first column
metadataMS <- metadataMS %>%
  relocate(Patient.ID)
metadataMS
##Now Patient.ID is the first column on the sheet (just easier for me to check this way)
metadataMS %>%
  arrange(Patient.ID)
##This arranges samples in alphabetical order on the metadata sheet.
####Actually does it? Check that. 
##Now I need to arrange the rows in alphabetical order on the count sheet

##trying a different way
metadataMS <- metadataMS[order(metadataMS$Patient.ID),]
metadataMS
##This worked, now using the same method to order metadata:

countdataMS <- countdataMS[ , order(names(countdataMS))]    # Apply order & names
countdataMS                                   
##This rearranged the sheet, but what arrange() and order(names) accomplished 
##is different. Alphabetical order is differently interpreted. 

##Dplyr method- didn't work as well as the df[,order(names)df] method aboves, saving for reference
{##lets try the dplyr method 

countdataMS <- countdataMS %>%
  select(sort(names(.)))
countdataMS
##Still the same issue. 
##I am sure there is a not super hard dplyr way around this. I will find it!
}

##Checking names
names(countdataMS)[-1]
metadataMS$Patient.ID
names(countdataMS)[-1]==metadataMS$Patient.ID
##Error- resolved. Needed to use the order function above to align dfs
{
##False all around, as expected. I need to rearrange the tibble.
##The issue is that the names are not in the same order on the metadata and countdata sheets
##Can I use the arrange() piping function to reorganize one of these tables?
##Can probably use arrange() and just sort by alphabetical order, this will let both datasets allign
##The sheets are both rearranged, but somehow sorted differently? Not quite alphabetical
}
#Now all readouts are true. 
all(names(countdataMS)[-1]==metadataMS$Patient.ID)

##We are all set. Let's do Deseq workflow and ultimately get a heatmap for all samples now. 

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(NMF)
library(factoextra)

ddsALL <- DESeqDataSetFromMatrix(countData=countdataMS, 
                                   colData=metadataMS, 
                                   design=~Group, 
                                   tidy=TRUE,
                                   ignoreRank = FALSE)
##Error here- cannot have "NA" in count matrix. Which sample is it? 
##I am going to have to get rid of the samples that have "N/A" counts 

countdataMS[is.na(countdataMS)]
##There is 1 NA count for some reason
countdataMS[is.na(countdataMS)] <- 0
##Converted the count that was "na" to "0"
countdataMS[is.na(countdataMS)]
##now it says character(0)- no more "na" counts

##Now the ddsALL line above works. 
##Deseq
ddsALL <- DESeq(ddsALL)

##Continuing the workflow
res <- results(ddsALL, tidy=TRUE)
res <- tbl_df(res)
res

res %>% 
  filter(padj<0.05) %>% 
  write_csv("metadataALLsig.csv")
##And it is saved to the X Chrom Project folder on my computer. 

##Now for plotting counts
plotCounts(ddsALL, gene="ENSG00000102245", intgroup="Group", returnData = TRUE)
plotCounts(ddsALL, gene="ENSG00000102245", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("CD40LG")


plotCounts(ddsALL, gene="ENSG00000147050", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("KDM6A")

##Plugged in the Ensembl codes for CD40LG and KDM6A - X Chromosome candidate genes. 

##I reran the code to get rid of the factoextra step, since that was throwing things off
##lesson learned- let Deseq make factors for itself, don't try to do that on your own
##This gives a pretty wimpy boxplot but has both control and lesion groups (only 3 patients here; 2 lesion, 1 control)

# Create the new column
res <- res %>% mutate(sig=padj<0.05)

# How many of each?
res %>% 
  group_by(sig) %>% 
  summarize(n=n())

plotMA(ddsALL)
plotMA(res)
##I'm not sure I'm doing the MA plot correctly. Should it be the Deseq object or the results df? 
res$sig <- res$padj < 0.05
plotMA(res[,c("baseMean","log2FoldChange","sig")], colLine = "tan2", colSig = "red", colNonSig = "midnightblue", xlim=c(200,10000))
##This fixes that. It needs to run on 3 columns, 2 are baseMean and log2FoldChange, the other is the logical "sig" (T/F/NA), which I created above as described by the padj column
##Below was from the markdown code from the tutorial (I apparently did this the hard way)
res %>% ggplot(aes(baseMean, log2FoldChange, col=sig)) + geom_point() + scale_x_log10() + ggtitle("MA plot")

##Volcano
res %>% ggplot(aes(log2FoldChange, -1*log10(pvalue), col=sig)) + geom_point() + ggtitle("Volcano plot- MS")

##The next step in the workflow is transformation. 

vsdataALL <- vst(ddsALL, blind=FALSE)

plotPCA(vsdataALL, intgroup="Group")

##Lets try with log transformation
rld <- rlogTransformation(ddsALL)
plotPCA(rld, intgroup="Group")
##When I run this, I get: 
##rlog() may take a long time with 50 or more samples,
##vst() is a much faster transformation
##Did vst in the previous step. rld is needed for the heatmap step I believe?  

sessionInfo()
##That concludes the tutorial pretty much. Now the challenge is getting all the metadata
##Lined up with all the patient samples

##For Heatmaps:
install.packages("NMF")
library(NMF)

aheatmap(assay(rld)[arrange(res, padj, pvalue)$row[1:50], ], 
         labRow=arrange(res, padj, pvalue)$symbol[1:50], scale="row", distfun="pearson", annCol=dplyr::select(metadataALL, Group, gender), 
         col=c("steelblue1","wheat1","wheat1","red2")) -> heatmap1
heatmap1

##Note that the graphics in R cannot display this. Using the below saves the output plot at a pdf in Documents
##can also save as a png by changing the first function to png(file = "name.png"), PDF had better resolution. 

pdf(file = "Heatmap1deseq.pdf")
aheatmap(assay(rld)[arrange(res, padj, pvalue)$row[1:50], ], 
         labRow=arrange(res, padj, pvalue)$symbol[1:50], scale="row", distfun="pearson", annCol=dplyr::select(metadataALL, Group, gender), 
         col=c("steelblue1","wheat1","wheat1","red2"))
dev.off()

##We have a heatmap. The next step- can we pull out just the X chromsome genes from the count data, and do this all over again? 

countdataMS
##My instinct is to use the join functions in dplyr that James described, matched to other data that lists X chrom genes.
##Also, we can filter with significance. From this workflow I already know significant p-value genes: 

sigpMS <- read.csv("metadataALLsig.csv")
sigpMS

##That said, I would also be interested to know about non-sig X chrom genes, so we can see if they ARE significant in SLE
##Can I use a reference genome (example: hg38) that is already annotated, and pull out all X chrom? 
file.choose()
##I found the following file on the ensembel database for hg38, ChromX: 
##"C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\Homo_sapiens.GRCh38.106.chromosome.X.gff3.gz"
##I have no idea what a .gff3 file is.
##Looked it up- this is like an in depth explorer of the genes with positions, 5' or 3' reads, introns/exons, etc.
##Probably not exactly what I'm looking for but useful for later. 

##the BioMart database might be the way to go here. 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl
filters = listFilters(ensembl)
filters[1:5,]

##This is promising. Chromosome is first filter function of many. 
genenamesMS <- countdataMS %>%
  select(columns = ...1)
genenamesMS
##These are all the gene names in the MS count matrix. 
##Next step is to get all the X chrom names out. Can I use biomaRt to get a list? 
##ensembl_gene_id  and chromosome_name


##Package:ensembldb
BiocManager::install("ensembldb")
library(ensembldb)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
edb
organism((edb))
##Homo sapiens
##Genome build is GRCh38

supportedFilters(edb)
XIDs <- transcriptsBy(edb, by = "gene", filter = SeqNameFilter(c("X")))
XIDs
##We have a list
##This package can do way more stuff as well- pulling out sequences, transcript, gene/exon models, etc. https://rdrr.io/bioc/ensembldb/f/vignettes/ensembldb.Rmd
##I don't know what a GRanges file is though. Need to make it a data frame so I understand? 
as_tibble(XIDs) -> XIDsdf
XIDsdf
##That worked
XIDsdf %>%
  dplyr::select(columns = gene_id, seqnames) -> XChromGenesdf
XChromGenesdf
##And we have a list we can work from 
##Now we use join functions to pull out the XChrom genes from the MS count data
##ISSUE- there are duplicate gene counts. Sorted by "group" on this df
##Can we eliminate the duplicates? I don't know if this will affect count data in deseq or not.

##6.29.22
library(dplyr)
XcountdataMS <- left_join(XChromGenesdf,
                        countdataMS,
                        by = c("columns" = "...1"),
                        suffix = c("_X","_Chrom"),
                        keep = FALSE,
                        na_matches = "never")
XcountdataMS

XcountdataMS %>%
  dplyr::filter(seqnames == "X")

##This made a tibble that is 6089 genes long- same number as on the XChromgenesdf tibble
##For it to work in Deseq2, I think I need to remove the seqnames column
##ISSUE- there are duplicates in the genes listed. Need to clean up the XIDsdf data frame first. 

##Consolidating the duplicate rows for some genes into 1 row because I only really care about the gene ID
XIDsdf <- XIDsdf %>%
  distinct(group_name, .keep_all = TRUE)
##This cleans it up so that its just one gene per row
XIDsdf

##Now cleaning up the conversion df
XIDsdf %>%
  dplyr::select(columns = gene_id, seqnames) -> XChromGenesdf
XChromGenesdf

##left join
XcountdataMS <- left_join(XChromGenesdf,
                          countdataMS,
                          by = c("columns" = "...1"), ##Match "columns" column to the ".001" column in countdata. These both correspond ton Ensembl gene IDs
                          suffix = c("_X","_Chrom"), ##I don't fully understand what this line does, but it doesn't seem to affect my output
                          keep = FALSE,
                          na_matches = "never")
XcountdataMS
##Now 2,389 distinct rows. To make sure:
XcountdataMS %>%
  distinct(columns, .keep_all = TRUE)
##Is the same, tibble 2399 x 70. 

##Removing the seqnames column so we can process with DEseq2
XcountdataMSFINAL <- dplyr::select(XcountdataMS, -c("seqnames"))
XcountdataMSFINAL
##Now the DEseq2 workflow.

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(NMF)
library(factoextra)

##First format as data frame for deseq2
countdataXMS <- as.data.frame(XcountdataMSFINAL)


##Remove "NA" counts
countdataXMS[is.na(countdataXMS)]
##There is 1 NA count for some reason
countdataXMS[is.na(countdataXMS)] <- 0
##Converted the count that was "na" to "0"
countdataXMS[is.na(countdataXMS)]
##now it says character(0)- no more "na" counts

ddsX <- DESeqDataSetFromMatrix(countData=countdataXMS, 
                                 colData=metadataMS, 
                                 design=~Group, 
                                 tidy=TRUE,
                                 ignoreRank = FALSE)
##Deseq
ddsX <- DESeq(ddsX)

##Continuing the workflow
resX <- results(ddsX, tidy=TRUE)
resX <- tbl_df(resX)
resX

resX %>% 
  dplyr::filter(padj<0.05) %>% 
  write_csv("XchromONLYsig.csv")
##And it is saved to the X Chrom Project folder on my computer. 

##Now for plotting counts
plotCounts(ddsX, gene="ENSG00000102245", intgroup="Group", returnData = TRUE)
plotCounts(ddsX, gene="ENSG00000102245", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("CD40LG")


plotCounts(ddsX, gene="ENSG00000147050", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("KDM6A")

##Plugged in the Ensembl codes for CD40LG and KDM6A - X Chromosome candidate genes. 
##These two charts look the same as they did when I examined the whole genome sample- and they should. They are both on the X chromsome. 
##Let's try a non-X chromosome gene. It should NOT work for this sample.
##Lactate Dehydrogenase A (LDHA)- Chromosome 11
plotCounts(ddsX, gene="ENSG00000134333", intgroup="Group", returnData=TRUE) %>% 
  ggplot(aes(Group, count)) + geom_boxplot(aes(fill="Samp")) + scale_y_log10() + ggtitle("LDHA")
##Error in counts(dds, normalized = normalized, replaced = replaced)[gene,  : 
##subscript out of bounds
##That is what I was hoping for! 


# Create the new column
resX <- resX %>% mutate(sig=padj<0.05)

# How many of each?
resX %>% 
  group_by(sig) %>% 
  summarize(n=n())

plotMA(ddsX)

resX$sig <- resX$padj < 0.05
plotMA(resX[,c("baseMean","log2FoldChange","sig")], colLine = "tan2", colSig = "red", colNonSig = "midnightblue", xlim=c(200,10000))
##This fixes that. It needs to run on 3 columns, 2 are baseMean and log2FoldChange, the other is the logical "sig" (T/F/NA), which I created above as described by the padj column
##Below was from the markdown code from the tutorial (I apparently did this the hard way)
resX %>% ggplot(aes(baseMean, log2FoldChange, col=sig)) + geom_point() + scale_x_log10() + ggtitle("MA plot")

##Volcano
resX %>% ggplot(aes(log2FoldChange, -1*log10(pvalue), col=sig)) + geom_point() + ggtitle("Volcano plot")

##The next step in the workflow is transformation. 

vsdataX <- vst(ddsX, blind=FALSE)
##This gave a warning message, but we can use plotPCA from rlogTransformation

plotPCA(vsdataX, intgroup="Group")

##Lets try with log transformation
rldX <- rlogTransformation(ddsX)
plotPCA(rldX, intgroup="Group")

##For Heatmaps:
library(NMF)

aheatmap(assay(rldX)[arrange(resX, padj, pvalue)$row[1:20], ], 
         labRow=arrange(resX, padj, pvalue)$symbol[1:20], scale="row", distfun="pearson", annCol=dplyr::select(metadataALL, Group, gender), 
         col=c("steelblue1","wheat1","wheat1","red2")) -> heatmapX
heatmapX
##Error: Error in assay(rldX)[arrange(res, padj, pvalue)$row[1:50], ] : 
##subscript out of bounds
##Fixed- had res instead of resX in the first line

jpeg(file = "HeatmapXdeseq1.jpeg")
aheatmap(assay(rldX)[arrange(resX, padj, pvalue)$row[1:20], ], 
         labRow=arrange(resX, padj, pvalue)$symbol[1:20], scale="row", distfun="pearson", annCol=dplyr::select(metadataALL, Group, gender), 
         col=c("steelblue1","wheat1","wheat1","red2"))
dev.off()

##And there we have it. A heatmap of the top 20 differentially expressed X chromosome genes in this dataset of MS patients.
sessionInfo()
