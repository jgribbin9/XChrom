##MS edgeR analysis
library(readr)
library(dplyr)
library(edgeR)
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
##This makes it just "MS" and "control" groups. Let's try running it

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

dim(dge.MS.Male) ## 36974, 28
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
topTags(lrtmale, n=100000) -> toptagsMSmale


##When I do this, there is no P value or FDR for  the "~0+Group" design matrix and I don't know why
##But when I change to just "~Group", it works. I continue to not be 100% sure what I'm telling R to do here but I know I need Pvals
as.data.frame(toptagsMSmale) -> toptagsdfMSmale
toptagsdfMSmale %>%
  dplyr::filter(PValue < 0.05) -> toptagsdfMSmale
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
topXchromsMale <- FullAnnoTopTagsMaleMS %>%
  filter(chromosome_name == "X")
topXchromsMale ##All X chromosomal genes, comparing Male patients to male controls

##females

metadataMS %>%
  filter(gender =="female") -> FEMALE.metadata.MS
c(NAWM="MS", Lesion="MS", Control = "control") -> conv1 ##Here, I reclassified all lesion types as MS, vs. controls who do not have MS
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

dim(dge.MS.Female) ##36974, 36
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
topTags(lrtfemale, n=100000) -> toptagsMSfemale
toptagsMSfemale


##When I do this, there is no P value or FDR for  the "~0+Group" design matrix and I don't know why
##But when I change to just "~Group", it works. I continue to not be 100% sure what I'm telling R to do here but I know I need Pvals
as.data.frame(toptagsMSfemale) -> toptagsdfMSfemale
toptagsdfMSfemale %>%
  dplyr::filter(PValue < 0.05) -> toptagsdfMSfemale
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
topXchromsFemale <- FullAnnoTopTagsFemaleMS %>%
  filter(chromosome_name == "X")
topXchromsFemale
write.csv(topXchromsFemale, )

##Testing
topXchromsFemale%>%filter(hgnc_symbol == "TSIX")
topXchromsMale%>%filter(hgnc_symbol == "TSIX")

##For the GOVenn function to work, I need df's that have the following columns: 
##genes      logFC      logCPM       LR       PValue          FDR id
topXchromsFemale %>%
  dplyr::select(ensembl_gene_id,logFC,logCPM,LR,PValue,FDR,id) ->XF2


topXchromsMale%>%
  dplyr::select(ensembl_gene_id,logFC,logCPM,LR,PValue,FDR,id) ->XM2


##I cross-referenced these results vs. what was on my poster last summer- they are consistent.


##Venn diagram
install.packages("GOplot")
library(GOplot)
GOVenn(toptagsdfMSfemale, toptagsdfMSmale, title = "MALE VS. FEMALE GENE EXPRESSION: MULTIPLE SCLEROSIS (ALL SIGNIFICANT GENES)", label = c("Female", "Male")) -> VennALL
GOVenn(XF2, XM2, title = "MALE VS. FEMALE GENE EXPRESSION: MULTIPLE SCLEROSIS (X CHROMOSOME ONLY)", label = c("Female", "Male"), plot=F) -> VennX ##See above for generating the "XF2" and "XM2" dfs that are needed for this to work


##Fixing Graphics for Venn Diagram
##Needs to be colorblind accessible
colors <- c("#CC79A7", "#56B4E9") ##needs a character vector to change colors
colorsinside <- c("#009E73", "#F0E442", "#D55E00")

GOVenn(XF2, XM2, title = "MALE VS. FEMALE GENE EXPRESSION: MULTIPLE SCLEROSIS (X CHROMOSOME ONLY)", label = c("Female", "Male"), plot=T, circle.col = colors, lfc.col = colorsinside) -> VennXCB
GOVenn(toptagsdfMSfemale, toptagsdfMSmale, title = "MALE VS. FEMALE GENE EXPRESSION: MULTIPLE SCLEROSIS (ALL SIGNIFICANT GENES)", label = c("Female", "Male"), plot=T, circle.col = colors, lfc.col = colorsinside) -> VennALLCB
##these were pulled from a stack overflow post about colorblind pallets- not sure how offical that is? https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

##Can I pull out the specific genes that are shown in this Venn diagram? In particular, the ones in the middle that are contra-regulated are of high interest. 
VennX$table -> VennXgenelist

VennX$table$A_only -> VennXfemaleGenesOnly ##DOUBLE CHECK THIS TO MAKE SURE THIS IS FEMALE!!! (I think it is becaue "A" is first and "female" is listed first)
as.data.frame(VennXfemaleGenesOnly) -> VXFG

VennX$table$B_only -> VennXmaleGenesOnly
as.data.frame(VennXmaleGenesOnly) -> VXMG

VennX$table$AB -> VennF.M_overlap
as.data.frame(VennF.M_overlap) -> VXOL

#############6.28.23
##We need to annotate the above to know what we are dealing with from a GO perspective
VXFG
rownames(VXFG) -> VXFG.genes
VXFG %>% mutate(new_col=VXFG.genes) -> VXFG
rownames(VXFG)
VXFG$new_col
all(rownames(VXFG) == VXFG$new_col)

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##Gives all possible annotation parameters
listAttributes(mart)
genes <- VXFG$new_col
VXFG.Anno <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                                 "description", "chromosome_name","hgnc_symbol"),values=genes,mart= mart)

VXMG
rownames(VXMG) -> VXMG.genes
VXMG %>% mutate(new_col=VXMG.genes) -> VXMG
rownames(VXMG)
VXMG$new_col
all(rownames(VXMG) == VXMG$new_col)

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##Gives all possible annotation parameters
listAttributes(mart)
genes <- VXMG$new_col
VXMG.Anno <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                             "description", "chromosome_name","hgnc_symbol"),values=genes,mart= mart)

VXOL
rownames(VXOL) -> VXOL.genes
VXOL %>% mutate(new_col=VXOL.genes) -> VXOL
rownames(VXOL)
VXOL$new_col
all(rownames(VXOL) == VXOL$new_col)

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##Gives all possible annotation parameters
listAttributes(mart)
genes <- VXOL$new_col
VXOL.Anno <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                             "description", "chromosome_name","hgnc_symbol"),values=genes,mart= mart)
##Same issue as before- adds EVERY GO term to each gene. Can (should?) I narrow this down to one or two?

##Exporting for IPA
write.csv(topXchromsFemale, "TOP.X.CHROMS.FEMALE.MS.csv")
write.csv(topXchromsMale, "TOP.X.CHROMS.MALE.MS.csv")

VXFG.Anno
VXMG.Anno
VXOL.Anno

left_join(VXFG.Anno, VXFG, by = c("ensembl_gene_id" = "new_col")) ->VXFG.Anno.Final
left_join(VXMG.Anno, VXMG, by = c("ensembl_gene_id" = "new_col")) ->VXMG.Anno.Final
left_join(VXOL.Anno, VXOL, by = c("ensembl_gene_id" = "new_col")) ->VXOL.Anno.Final

write.csv(VXFG.Anno.Final, "Female_Exclusive_XGenes.csv")
write.csv(VXMG.Anno.Final, "Male_Exclusive_XGenes.csv")
write.csv(VXOL.Anno.Final, "Common_MaleFemale_XGenes.csv")

##Tables

VXFG.Anno.Final
VXMG.Anno.Final
VXOL.Anno.Final

VXOL.Anno.Final[order(-VXOL.Anno.Final$logFC_A),] -> VXOL.Anno.Final
colnames(VXOL.Anno.Final) <- c("Ensembl GeneID", "Description", "Chromsome", "HGNC Symbol", "logFC-Female", "logFC-Male", "Trend")

library(kableExtra)

VXFG.Anno.Final%>%
  kbl(caption = "Unique Female Genes- X Chromosome")%>%
  kable_classic(full_width = F, html_font = "Cambria") -> FemaleExclusiveGenesChart

VXMG.Anno.Final%>%
  kbl(caption = "Unique Male Genes- X Chromosome")%>%
  kable_classic(full_width = F, html_font = "Cambria") -> MaleExclusiveGenesChart

VXOL.Anno.Final%>%
  kbl(caption = "Male/Female Overlap Genes- X Chromosome")%>%
  kable_classic(full_width = F, html_font = "Cambria") -> MaleFemaleOverlapGenesChart



##LAST STEP: USE LEFT_JOIN TO JOIN VOXL.ANNO AND OTHER ANNOS TO THEIR OG DATA FRAME (Example: left_join VOXL.ANNO, VOXL, by ensembel id [just need to refresh how to do this])
##THis will annotate all the genes with HGNC symbols and ensembel IDs. 



############

 ##Male and Female together
{
MS.Matrix <- as.matrix(countsonlyMStrimmed)
##Error: NA counts not allowed
##I remember from before- there was one single N/A count in this dataset. 
countdataMS[is.na(countdataMS)] <- 0
any(is.na(countdataMS))

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

##Let's try chromoMap
install.packages("chromoMap")
library(chromoMap)

##For this, I need a data frame with genes, labeled with chromosome start and end loci
##I think the bioMart anno has the ability to include this?

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##Gives all possible annotation parameters
listAttributes(mart)
genes <- toptagsdfMSfemale$genes
toptagsdfMSfemale$id <- NA
GenesAnnotatedFemaleMSTopTags <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                                 "description", "chromosome_name","hgnc_symbol", "start_position", "end_position"),values=genes,mart= mart)
##I am not sure if I should use "gene start" or "transcript start", both are avaliable with bioMart. Will try with gene start for now
##This is the same step as before, just now has start and end on the df. Adding to the rest

##To get Pvals, FDRs, logCPM scores etc on the same DF:
left_join(GenesAnnotatedFemaleMSTopTags,toptagsdfMSfemale, by =c("ensembl_gene_id"="genes")) -> FullAnnoTopTagsFemaleMSCHROMO
FullAnnoTopTagsFemaleMSCHROMO
write.table(FullAnnoTopTagsFemaleMS, file = "top.tags.MS.female.CHROMO.txt", sep = ",")
read.csv("top.tags.MS.female.CHROMO.txt") -> MStagsFCHROMO


##X Chrom isolation

top100XchromsFemaleCHROMO <- FullAnnoTopTagsFemaleMSCHROMO %>%
  filter(chromosome_name == "X")
top100XchromsFemaleCHROMO

##This df now has start_position and end_position for all female X chromosomal genes on it. This is a prerequisit for chromoMap
##Format. Columns must be organized 1.) Chromosome name, 2.) Chromosome start, 3.) chromosome end, 4.) centromere start (optional, excluded here)
top100XchromsFemaleCHROMO%>%
  dplyr::select(c(chromosome_name, start_position, end_position, hgnc_symbol, )) ->FemaleChromo1
##this has the info. 
##Need to update: chromosome start and chromosome end should specify the beginning and end of the entire X chromosome
##Start and end are the "element start" and "element end" parameters
##According to https://www.ncbi.nlm.nih.gov/grc/human/data, the X chromosome is 156,040,895 BP long
##so I need a data frame that says X chrome, 1, and 156,040,895
cbind(top100XchromsFemaleCHROMO, Start=c(1), End=c(156040895)) -> top100XchromsFemaleCHROMO
##easy enough. now to reformat
top100XchromsFemaleCHROMO%>%
  dplyr::select(c(chromosome_name, Start, End, hgnc_symbol, chromosome_name, start_position, end_position, logFC)) ->FemaleChromo1

chr_file_X <- top100XchromsFemaleCHROMO%>%
  dplyr::select(c(chromosome_name, Start, End))
dim(chr_file_X)
##I think I only need the first row of this df
chr_file_X_short <- chr_file_X[1,]

anno_file_X <- top100XchromsFemaleCHROMO%>%
  dplyr::select(c(hgnc_symbol, chromosome_name, start_position, end_position, logFC))
dim(anno_file_X)
##The vingette says I need to drop the headers

chromoMap(list(chr_file_X_short),list(anno_file_X)) -> p1
##seems to work this way

colorlist <- list(c("blue", "white", "red"))

chromoMap(list(chr_file_X_short),list(anno_file_X),
          data_based_color_map = T, data_colors = colorlist,
          data_type = "numeric") -> p2


chromoMap(list(chr_file_X_short),list(anno_file_X),
          data_based_color_map = T, data_colors = colorlist,
          data_type = "numeric", chr.2D.plot = T) -> p3

##Let's get all X chrom genes for males and for females


##IL13RA1 jumps out as a potential gene of interest. It is also right next to NKRF, which is increased in expression (this should be ANTI inflammatory)
top100XchromsFemale %>%
  dplyr::filter(hgnc_symbol == "IL13RA1", "NKRF", "TLR7", "TLR8")
##Other biologic things to point out: TSIX expression appears to be down while XIST expression is slightly up
##TLR7 and TLR8 expression are both up. These genes are neighbors on the X chromosome. TLR7 is a known escape gene


##Trial GO annotation
FemTopTagsGOTRIAL <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                                 "description", "chromosome_name","hgnc_symbol", "go_id"),values=genes,mart= mart)
##why does it look like that? I think it added all associated GO ids with each gene? and that's why the dataset expanded so much?
listAttributes(mart)

FemTopTagsGOTRIAL <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                                     "description", "chromosome_name","hgnc_symbol", "name_1006", "definition_1006", "go_linkage_type"),values=genes,mart= mart)
##GO segregation trial
VXFG.Anno %>% filter(name_1006 == "immune system process")
VXFG.Anno %>% filter(name_1006 == "innate immune response")
VXFG.Anno %>% filter(name_1006 == "immune response")

##3.23.24 goanna
goana(lrtmale, geneid=rownames(lrtmale)) ->goana.Male
goana(lrtfemale, geneid=rownames(lrtfemale)) -> goana.Female
topGO(goana.Male) -> topGOMale
write.csv(topGOMale, "topGOMale.csv")
topGO(goana.Female) -> topGOFemale
write.csv(topGOFemale, "topGOFemale.csv")

topGOMale%>%
  kbl(caption = "Gene Ontology- Male MS vs. Control")%>%
  kable_classic(full_width = F, html_font = "Cambria") -> TopGoMaleTable

topGOFemale%>%
  kbl(caption = "Gene Ontology- Female MS vs. Control")%>%
  kable_classic(full_width = F, html_font = "Cambria") -> TopGoFemaleTable
