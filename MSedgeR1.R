##MS edgeR analysis
library(readr)
library(dplyr)
countDataMS <-read_delim(("C:\\Users\\jgrib\\Desktop\\X Chromosome Project\\GSE179427_countmtx.csv"),delim = ';')
countdataMS <- as.data.frame(countDataMS)
rename(countdataMS, "genes" = "...1") ##just renaming first column
##Metadata- from previous
{
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

metaDataMS %>%
  rename("Patient.ID" = "Patient ID") -> metadataMS
}
metadataMS
##Need count matrix, eliminate gene names and have raw counts:
countdataMS[-1] -> countsonlyMS
countsonlyMS
select(countsonlyMS,-c("S07_254","S09_301","S11_110","S12_122")) -> countsonlyMStrimmed ##now we have trimmed out the gender NA patients

##Here, we will make seperate data frames for males and females
##First- find all males
metadataMS %>%
  filter(gender == "male")-> MALE.metadata.MS

MALE.metadata.MS
c(NAWM="MS", Lesion="MS", Control = "control") -> conv1
MALE.metadata.MS$Group <- as.character(conv1[MALE.metadata.MS$Group])
##This makes it just "MS" and "control" groups. Let's try running again

##need patient names
MALE.metadata.MS$Patient.ID
list(MALE.metadata.MS$Patient.ID)

##then need to get the counts for just these patient IDs
countsonlyMStrimmed %>%
  select(MALE.metadata.MS$Patient.ID) -> male.counts.MS
colnames(male.counts.MS)
MALE.metadata.MS$Patient.ID
##Checking alignment
all(names(male.counts.MS)==MALE.metadata.MS$Patient.ID)
##All True
##na.counts
male.counts.MS[is.na(male.counts.MS)] <-0
##Single N.A count is now 0. 

##Now matrix 
MS.Matrix.Male <- as.matrix(male.counts.MS)


dge.MS.Male <- DGEList(counts=MS.Matrix.Male, genes=countdataMS[,1])

dim(dge.MS.Male)
dge.MS.Male$samples$lib.size <- colSums(dge.MS.Male$counts)
rownames(dge.MS.Male$counts) <- rownames(dge.MS.Male$genes) 
y3 <- calcNormFactors(dge.MS.Male)
y3$samples
barplot(y3$samples$lib.size*1e-6, names=1:28, ylab="Library size (millions)")

designMS <- model.matrix(~Group,data=MALE.metadata.MS)
y3 <- estimateDisp(y3, designMS, robust=TRUE)
y3$common.dispersion
plotBCV(y3)

fitmale <- glmFit(y3, designMS)
lrtmale <- glmLRT(fitmale)
topTags(lrtmale, n=1000) -> toptagsMSmale
toptagsMSmale

##When I do this, there is no P value or FDR for  the "~0+Group" design matrix and I don't know why
##But when I change to just "~Group", it works. I continue to not be 100% sure what I'm telling R to do here but I know I need Pvals
as.data.frame(toptagsMSmale) -> toptagsdfMSmale
toptagsdfMSmale$genes

o3 <- order(lrtmale$table$PValue)
cpm(y3)[o3[1:10],]
summary(decideTests(lrtmale))

plotMD(lrtmale)
abline(h=c(-1, 1), col="blue")

##need Venn diagram package to compare male to female
##also need to annotate gene names
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##Gives all possible annotation parameters
listAttributes(mart)
genes <- toptagsdfMSmale$genes
toptagsdfMSmale$id <- NA
GenesAnnotatedMaleMSTopTags <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                         "description", "chromosome_name","hgnc_symbol"),values=genes,mart= mart)
##To get Pvals, FDRs, logCPM scores etc on the same DF:
left_join(GenesAnnotatedMaleMSTopTags,toptagsdfMSmale, by =c("ensembl_gene_id"="genes")) -> FullAnnoTopTagsMaleMS
FullAnnoTopTagsMaleMS
write.table(FullAnnoTopTagsMaleMS, file = "top.tags.MS.Male.txt", sep = ",")

GenesAnnotatedMaleMSTopTags
top100XchromsMale <- FullAnnoTopTagsMaleMS %>%
  filter(chromosome_name == "X")
top100XchromsMale


##in the top 100 genes regulated differently for the NAWM group, 5 are on the X chromosome for males (I think this is what I can conclude?)
##of these, IGSF1 is probably the most interesting- immunoglobulin superfamily member 1.
##Also HIV tat factor 1- does someone in this dataset have HIV? 
##I would expect very different results for females. 
##I think I may want to run this again as just "control" vs "MS" so that I'm not throwing out the lesion group in consideration. This would entail just chaning the names "NAWM" and "Lesion" to "MS" in the metadata df. 



##This annotates the top 100 genes, but I would like to be a little more specific about Xchrom if I can. Should talk over what to expect at meeting gene wise. As these are Males.

##females
library(dplyr)
metadataMS %>%
  filter(gender =="female") -> FEMALE.metadata.MS
c(NAWM="MS", Lesion="MS", Control = "control") -> conv1
FEMALE.metadata.MS$Group <- as.character(conv1[FEMALE.metadata.MS$Group])
##patient names
FEMALE.metadata.MS$Patient.ID
countsonlyMStrimmed %>%
  dplyr::select(FEMALE.metadata.MS$Patient.ID) -> female.count.MS
##Checking allignment
all(names(female.count.MS)==FEMALE.metadata.MS$Patient.ID)
##All True

##Matrix
MS.Matrix.Female <- as.matrix(female.count.MS)


dge.MS.Female <- DGEList(counts=MS.Matrix.Female, genes=countdataMS[,1])

dim(dge.MS.Female)
dge.MS.Female$samples$lib.size <- colSums(dge.MS.Female$counts)
rownames(dge.MS.Female$counts) <- rownames(dge.MS.Female$genes) 
y4 <- calcNormFactors(dge.MS.Female)
y4$samples
barplot(y4$samples$lib.size*1e-6, names=1:36, ylab="Library size (millions)")

designMS.F <- model.matrix(~Group,data=FEMALE.metadata.MS)
y4 <- estimateDisp(y4, designMS.F, robust=TRUE)
y4$common.dispersion
plotBCV(y4)

fitfemale <- glmFit(y4, designMS.F)
lrtfemale <- glmLRT(fitfemale)
topTags(lrtfemale, n=1000) -> toptagsMSfemale
toptagsMSfemale


##When I do this, there is no P value or FDR for  the "~0+Group" design matrix and I don't know why
##But when I change to just "~Group", it works. I continue to not be 100% sure what I'm telling R to do here but I know I need Pvals
as.data.frame(toptagsMSfemale) -> toptagsdfMSfemale
toptagsdfMSfemale$genes

o4 <- order(lrtfemale$table$PValue)
cpm(y4)[o4[1:10],]
summary(decideTests(lrtfemale))

plotMD(lrtfemale)
abline(h=c(-1, 1), col="blue")

##need Venn diagram package to compare male to female
##also need to annotate gene names
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##Gives all possible annotation parameters
listAttributes(mart)
genes <- toptagsdfMSfemale$genes
toptagsdfMSfemale$id <- NA
GenesAnnotatedFemaleMSTopTags <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                               "description", "chromosome_name","hgnc_symbol"),values=genes,mart= mart)
##To get Pvals, FDRs, logCPM scores etc on the same DF:
left_join(GenesAnnotatedFemaleMSTopTags,toptagsdfMSfemale, by =c("ensembl_gene_id"="genes")) -> FullAnnoTopTagsFemaleMS
FullAnnoTopTagsFemaleMS
write.table(FullAnnoTopTagsFemaleMS, file = "top.tags.MS.female.txt", sep = ",")
read.csv("top.tags.MS.female.txt") -> MStagsF


##X Chrom isolation
GenesAnnotatedFemaleMSTopTags
top100XchromsFemale <- FullAnnoTopTagsFemaleMS %>%
  filter(chromosome_name == "X")
top100XchromsFemale


##At one point we may want to call all MS patients just "MS" rather than separate groups of MS patients

##Will need individual matrix for M and F as well. 

##Venn diagram
install.packages("GOplot")
library(GOplot)
GOVenn(toptagsdfMSfemale, toptagsdfMSmale, title = "MALE VS. FEMALE GENE EXPRESSION: MULTIPLE SCLEROSIS (TOP 1000 TAGS)", label = c("Female", "Male"))

##Male and Female together
{
MS.Matrix <- as.matrix(countsonlyMStrimmed)
##Error: NA counts not allowed
##I remember from before- there was one single N/A count in this dataset. 
countdataMS[is.na(countdataMS)] <- 0

{##What I would like to do at this point is make sure that the counts align with metadata (trimmed)
names(countsonlyMStrimmed)
metadataMStrimmed$Patient.ID
all(names(countsonlyMStrimmed)==metadataMStrimmed$Patient.ID)
##True. So the dfs align.I don't know if this is how edgeR lines up the design matrix, but this step doesn't take much time and is a carryover from the DEseq2 workflow so I did it just in case.
} ##Aligning metadata names to count df names

dge.MS <- DGEList(counts=MS.Matrix, genes=countdataMS[,1])

##Then rerun counts only and as.matrix above. After this dge.MS works. 
##NOTE- only 1 count is N/A- IMPUTE OR NO?
dim(dge.MS)
dge.MS$samples$lib.size <- colSums(dge.MS$counts)
rownames(dge.MS$counts) <- rownames(dge.MS$genes) 
y2 <- calcNormFactors(dge.MS)
y2$samples
barplot(y2$samples$lib.size*1e-6, names=1:64, ylab="Library size (millions)")
##A few of these are kind of small

plotMDS(y2) ##this step takes a while, and I'm not really sure how worth it it is? 

nrow(designMS)
rownames(designMS)
metadataMS[is.na(metadataMS$gender),] %>%
  select(Patient.ID) -> NAgenderpts ##these are the four patients that have "NA" gender 

##error: there are patients who have "NA" gender
metadataMS %>%
  filter(is.na(metadataMS$gender)) ##this returns the 4 patients as well. I want to cut them out of the data frame
subset(metadataMS, gender == "male" | gender == "female") -> metadataMStrimmed
##cool that worked! This is saying "make a subset of the data frame, using all results where gender is male OR (|) gender is female (excludes the NA values)
##getting rid of gender = NA patients from the metadataMS df
##but now I realize- its the y2 matrix that has 4 extra rows not the metadata sheet. 
##this recapitulates the same problem as below- how do I know what rows in y2 correspond to which patients? And are they lined up correctly?
##I could go all the way back to the source, and eliminate the 4 patient IDs from the count df as well and limit to 64. Then rerun all the code. 
##I can do this now, but it seems not the best way. I need to get better at understanding Matrices. 
ncol(y2)
colnames(y2)
##I think this happened because on the metadata there were a few rows that had "NA" for gender.
##I can either redesign the matrix, or I can eliminate those 4 patients from the counts (nrow 64, ncol 68 above, so 4 patients unaccounted for)
##On eye test, this seems true. For example, rowname "9" is missing in the design matrix. 
##I'd like to talk about this design matrix step as it is crucial and I don't understand it well.
##For example: how are the counts and the patient names lining up? I dont understand how
##Is it just lining up the metadata sheet with the counts? because I don't know if they are in the same order
##Also, blocking factors, etc.- what am I really telling R to do here? 

##Design Matrix
designMS <- model.matrix(~0+gender+Group,data=metadataMStrimmed)
y2 <- estimateDisp(y2, designMS, robust=TRUE)
y2$common.dispersion
plotBCV(y2)

##Dif expression - glmFit
fit2 <- glmFit(y2, designMS)
lrt2 <- glmLRT(fit2)
topTags(lrt2, n=20) -> toptagsMS
toptagsMS
as.data.frame(toptags) -> toptagsdfMS
toptagsdfMS$genes
toptagsdfMS

o2 <- order(lrt2$table$PValue)
cpm(y2)[o2[1:10],]
summary(decideTests(lrt2))

plotMD(lrt2)
abline(h=c(-1, 1), col="blue")

##The analysis here is done with the label NAWM, which is one of the two MS dx groups (compared to control). I'm not sure it this means we're only examining this group?
}
