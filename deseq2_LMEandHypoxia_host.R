library("DESeq2") 
library("arrayQualityMetrics") 
library("pheatmap") 
library("VennDiagram") 
library("tidyverse")
library("vegan") 
library("ape")
library(plyr)
library(Rmisc)
library(Venn)
library(gridExtra)
library(ggprism)
library(ggpubr) 
library(dplyr) 
library(cowplot)
library(base)


# Load count data from LME---------------------------------------------------------------
countdata <- read.delim("LME_all_sample_count_ofav.tsv")
countdata <- countdata %>% filter(!grepl("__",gene_id)) %>%
  column_to_rownames("gene_id")

head(countdata) 
length(countdata[,1])
# 32574 isogroups mapped

# split coral and symbiont ---------------
head(countdata)
tail(countdata)

# make conditions table --------
conds <- data.frame(samID = names(countdata))

conds <- conds %>%
  mutate(species = case_when(
    str_detect(samID, "FR") ~ "frank",
    str_detect(samID, "OF") ~ "fav"
  )) %>%
  mutate(lesion = case_when(
    str_detect(samID, "U") ~ "U",
    str_detect(samID, "AH") ~ "AH",
    str_detect(samID, "AL") ~ "AL"
  )) %>%
   mutate(bank = case_when(
    str_detect(samID, "W") ~ "west",
    TRUE ~ "east"
  )) %>%
  mutate(genet = paste(species, bank, as.numeric(str_extract(samID, "\\d+")), sep="_"))

# check sample order
table(conds$samID == names(counts))

#----------Total counts?
totalCounts = colSums(counts)
head(totalCounts)
totalCounts
min(totalCounts) #92271
mean(totalCounts) #592082.5
max(totalCounts)  #1085198

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conds,
  design = ~ species + lesion + bank)

# Set base mean minimum ---------------------------------------------------
means <- apply(counts,1,mean)
table(means>3)
# FALSE  TRUE 
# 21209 11365  

means3 <- names(means[means>3])
head(means3)
length(means3)
#11365

coral.countFilt <- counts[row.names(counts) %in% means3,]
head(coral.countFilt)

totalCountsFilt <- colSums(coral.countFilt)
totalCountsFilt

min(totalCountsFilt) #91181
max(totalCountsFilt) #1074174
mean(totalCountsFilt) #584525.4

# check sample order
table(names(coral.countFilt) == as.vector(conds$sam))
# TRUE 76

# Reconstruct data object (filtered) ------------------------------
ddsFilt <- DESeqDataSetFromMatrix(
  countData = coral.countFilt,
  colData = conds,
  design = ~ species + lesion + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("species","lesion","bank"), force=T)

# outliers that fail 2 tests
out2 <- c("FRAAH2", "OFAH3", "OFAH4", "OFAH5")

# remove out2
dim(coral.countFilt)
coral.count.out <- coral.countFilt %>% select(-one_of(c(out2)))
dim(coral.count.out)
# 11365    72

dim(conds)
head(conds)
conds.out <- conds %>% filter(!samID %in% c(out2))
dim(conds.out)
# 72 5

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = coral.count.out,
  colData = conds.out,
  design = ~ species + lesion + bank)

# remake vsd
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)


#---Results
# log2 fold change (MLE): species FR vs OF
resSpecies = results(deds, independentFiltering = F, contrast=c("species","frank","fav"))
resSpecies
table(resSpecies$padj<0.05)
# outliers removed: 1415

# log2 fold change (MLE): bank west vs east 
resBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
resBank
table(resBank$padj<0.05)
# outliers removed: 41

# log2 fold change (MLE): lesion AL vs AH 
resALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
resALvAH
table(resALvAH$padj<0.05)
# outliers removed: 181 

#log2 fold change (MLE): lesion AL vs U 
resALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
resALvU
table(resALvU$padj<0.05)
# outliers removed: 499

#log2 fold change (MLE): lesion AH vs U 
resAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
resAHvU
table(resAHvU$padj<0.05)
# outliers removed: 0

# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# outliers removed
# FALSE  TRUE 
# 11121   244 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE outliers removed
# 11245   120 

table(abs(resALvAH$log2FoldChange)>1.5)
# FALSE  TRUE outliers removed
# 11152   213

table(abs(resAHvU$log2FoldChange)>1.5)
# FALSE  TRUE outliers removed
# 11261   104

table(abs(resALvU$log2FoldChange)>1.5)
# FALSE  TRUE outliers removed
#  11069   296 

# differential expression as percentage of # of reads ----

isogroupsMapped <- nrow(counts)

# how many DEGs as percentage of transcriptome size?
round((table(resBank$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.13 
round((table(resSpecies$padj<0.05)[2]/isogroupsMapped*100),2)
# 4.34
round((table(resALvAH$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.56 
round((table(resALvU$padj<0.05)[2]/isogroupsMapped*100),2)
# 1.53 
round((table(resAHvU$padj<0.05)[2]/isogroupsMapped*100),2)
# NA (no DEGs) 

# Write results for GO/KOG analysis -------------------------------------------

# by LFC
LFCspecies <- data.frame(cbind("gene"=row.names(resSpecies),"LFC"=resSpecies$log2FoldChange))
write.table(LFCspecies,quote=F,row.names=F,file="GO_Species_LFC.csv",sep=",")

LFCbank <- data.frame(cbind("gene"=row.names(resBank),"LFC"=resBank$log2FoldChange))
write.table(LFCbank,quote=F,row.names=F,file="GO_Bank_LFC.csv",sep=",")

LFCAHvU <- data.frame(cbind("gene"=row.names(resAHvU),"LFC"=resAHvU$log2FoldChange))
write.table(LFCAHvU,quote=F,row.names=F,file="GO_AHvU_LFC.csv",sep=",")

LFCALvAH <- data.frame(cbind("gene"=row.names(resALvAH),"LFC"=resALvAH$log2FoldChange))
write.table(LFCALvAH,quote=F,row.names=F,file="GO_ALvAH_LFC.csv",sep=",")

LFCALvU <- data.frame(cbind("gene"=row.names(resALvU),"LFC"=resALvU$log2FoldChange))
write.table(LFCALvU,quote=F,row.names=F,file="GO_ALvU_LFC.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(resSpecies),"logP"=round(-log(resSpecies$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resSpecies$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4497 5000 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Species_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resBank),"logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4905 4592  
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Bank_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resAHvU),"logP"=round(-log(resAHvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resAHvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4458 5039 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_AHvU_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resALvAH),"logP"=round(-log(resALvAH$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resALvAH$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 5365 4132   
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_ALvAH_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resALvU),"logP"=round(-log(resALvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resALvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 5103 4394 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_ALvU_logP.csv",sep=",")

# Subset data by species ------------------------------------------------------

OFconds= conds.out[conds.out$species== "fav",]
dim(OFconds) #33 5
FRconds= conds.out[conds.out$species== "frank",]
dim(FRconds) #39  5

OFcoral.count.out <- coral.count.out[, OFconds$sam]
all((OFconds$sam) == colnames(OFcoral.count.out))
dim(OFcoral.count.out) #11365    33

FRcoral.count.out <- coral.count.out[, FRconds$sam]
all((FRconds$sam) == colnames(FRcoral.count.out))
dim(FRcoral.count.out) #11365    39


# Run DESeq for Orbicella faveolata ------------------------------

OFddsFiltOut <- DESeqDataSetFromMatrix(
  countData = OFcoral.count.out,
  colData = OFconds,
  design = ~ lesion + bank)

# remake vsd
OFvsd <- varianceStabilizingTransformation(OFddsFiltOut, blind=TRUE)

deds <- DESeq(OFddsFiltOut)

# log2 fold change (MLE): bank west vs east 
OFresBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
OFresBank
table(OFresBank$padj<0.05)
# 298

# log2 fold change (MLE): lesion AL vs AH 
OFresALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
OFresALvAH
table(OFresALvAH$padj<0.05)
# 64

#log2 fold change (MLE): lesion AL vs U 
OFresALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
OFresALvU
table(OFresALvU$padj<0.05)
# 160

#log2 fold change (MLE): lesion AH vs U 
OFresAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
OFresAHvU
table(OFresAHvU$padj<0.05)
#0

vals= cbind(OFresBank$pvalue, OFresBank$padj, OFresALvAH$pvalue, OFresALvAH$padj, OFresALvU$pvalue, OFresALvU$padj, OFresAHvU$pvalue, OFresAHvU$padj)
head(vals)
colnames(vals)=c("pval.b", "padj.b","pval.alah", "padj.alah","pval.alu", "padj.alu","pval.ahu", "padj.ahu")
length(vals[,1])
table(complete.cases(vals))
OFvsdpvals=cbind(assay(OFvsd), vals)
head(OFvsdpvals)
dim(OFvsdpvals)
# 11365    41
table(complete.cases(OFvsdpvals))

OFvsdpvals =as.data.frame(OFvsdpvals)
transcript <- rownames(OFvsdpvals)
OFvsdpvals $transcript= transcript #make a new column that is transcript so it can be merged with the genenames file
dim(OFvsdpvals)

# Load gene name annotations
gg <- read.delim("OfavGenome_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs
OFvsdpvals =merge(OFvsdpvals, gg, by="transcript", all=T)
dim(OFvsdpvals)
head(OFvsdpvals)
OFvsdpvals = OFvsdpvals[!is.na(OFvsdpvals$OFAH1),] #removes all the rows with no GE data
dim(OFvsdpvals)

write.csv(OFvsdpvals, "fgbdieoff_bm3_wald_VSDandPVALS_OFhost_geneNames_MS.csv", quote=F)

#output for GO and KOG
#AL v. AH. for Orbicella faveolata
OF_LFCALvAH <- data.frame(cbind("gene"=row.names(OFresALvAH),"LFC"= OFresALvAH $log2FoldChange))
write.table(OF_LFCALvAH,quote=F,row.names=F,file="GO_ALvAH_LFC_OF.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(OFresALvAH),"logP"=round(-log(OFresALvAH$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[OFresALvAH$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4759 4738 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_OF_ALAH_logP.csv",sep=",")

#AL v. U. for Orbicella faveolata
OF_LFCALvU <- data.frame(cbind("gene"=row.names(OFresALvU),"LFC"= OFresALvU $log2FoldChange))
write.table(OF_LFCALvU,quote=F,row.names=F,file="GO_ALvU_LFC_OF.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(OFresALvU),"logP"=round(-log(OFresALvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[OFresALvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4961 4536 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_OF_ALU_logP.csv",sep=",")

#AH v. U. for Orbicella faveolata
OF_LFCAHvU <- data.frame(cbind("gene"=row.names(OFresAHvU),"LFC"= OFresAHvU $log2FoldChange))
write.table(OF_LFCAHvU,quote=F,row.names=F,file="GO_AHvU_LFC_OF.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(OFresAHvU),"logP"=round(-log(OFresAHvU $pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[OFresAHvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
#4912 4585  
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_OF_AHU_logP.csv",sep=",")

# Run DESeq for Orbicella franksi ------------------------------

FRddsFiltOut <- DESeqDataSetFromMatrix(
  countData = FRcoral.count.out,
  colData = FRconds,
  design = ~ lesion + bank)

# remake vsd
FRvsd <- varianceStabilizingTransformation(FRddsFiltOut, blind=TRUE)

deds <- DESeq(FRddsFiltOut)

# log2 fold change (MLE): bank west vs east 
FRresBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
FRresBank
table(FRresBank$padj<0.05)
# 2

# log2 fold change (MLE): lesion AL vs AH 
FRresALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
FRresALvAH
table(FRresALvAH$padj<0.05)
# 17

#log2 fold change (MLE): lesion AL vs U 
FRresALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
FRresALvU
table(FRresALvU$padj<0.05)
# 199

#log2 fold change (MLE): lesion AH vs U 
FRresAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
FRresAHvU
table(FRresAHvU$padj<0.05)
# 1

vals= cbind(FRresBank$pvalue, FRresBank$padj, FRresALvAH$pvalue, FRresALvAH$padj, FRresALvU$pvalue, FRresALvU$padj, FRresAHvU$pvalue, FRresAHvU$padj)
head(vals)
colnames(vals)=c("pval.b", "padj.b","pval.alah", "padj.alah","pval.alu", "padj.alu","pval.ahu", "padj.ahu")
length(vals[,1])
table(complete.cases(vals))
FRvsdpvals=cbind(assay(FRvsd), vals)
head(FRvsdpvals)
dim(FRvsdpvals)
# 11365    47
table(complete.cases(FRvsdpvals))

FRvsdpvals =as.data.frame(FRvsdpvals)
transcript <- rownames(FRvsdpvals)
FRvsdpvals $transcript= transcript #make a new column that is transcript so it can be merged with the genenames file
dim(FRvsdpvals)

# Load gene name annotations
gg <- read.delim("OfavGenome_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs
FRvsdpvals =merge(FRvsdpvals, gg, by="transcript", all=T)
dim(FRvsdpvals)
head(FRvsdpvals)
FRvsdpvals = FRvsdpvals[!is.na(FRvsdpvals$FRAAH10),] #removes all the rows with no GE data
dim(FRvsdpvals)

write.csv(FRvsdpvals, "fgbdieoff_bm3_wald_VSDandPVALS_FRhost_geneNames.csv", quote=F)

#output for GO and KOG
#AL v. AH. for Orbicella faveolata
FR_LFCALvAH <- data.frame(cbind("gene"=row.names(FRresALvAH),"LFC"= FRresALvAH $log2FoldChange))
write.table(FR_LFCALvAH,quote=F,row.names=F,file="GO_ALvAH_LFC_FR.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(FRresALvAH),"logP"=round(-log(FRresALvAH $pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[FRresALvAH $log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 5250 4247 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_FR_ALAH_logP.csv",sep=",")

#AL v. U. for Orbicella faveolata
FR_LFCALvU <- data.frame(cbind("gene"=row.names(FRresALvU),"LFC"= FRresALvU $log2FoldChange))
write.table(FR_LFCALvU,quote=F,row.names=F,file="GO_ALvU_LFC_FR.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(FRresALvU),"logP"=round(-log(FRresALvU $pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[FRresALvU $log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4949 4548 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_FR_ALU_logP.csv",sep=",")

#AH v. U. for Orbicella faveolata
FR_LFCAHvU <- data.frame(cbind("gene"=row.names(FRresAHvU),"LFC"= FRresAHvU $log2FoldChange))
write.table(FR_LFCAHvU,quote=F,row.names=F,file="GO_AHvU_LFC_FR.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(FRresAHvU),"logP"=round(-log(FRresAHvU $pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[FRresAHvU $log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
#4504 4993  
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_FR_AHU_logP.csv",sep=",")


# Explore with plots ------------------------------------------------------

#Sample distance heatmap
pheatmap(cor(assay(vsd)),border_color=NA, main="SampleHeatmap")

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds, main="Dispersion Plot Response")

#MA plot
quartz()
par(mfrow=c(2,3))
plotMA(resSpecies, ylim = c(-1, 1), main="MA Plot Species") 
plotMA(resBank, ylim = c(-1, 1), main="MA Plot Bank") 
plotMA(resALvAH, ylim = c(-1, 1), main="MA Plot ALvAH") 
plotMA(resALvU, ylim = c(-1, 1), main="MA Plot ALvU") 
plotMA(resAHvU, ylim = c(-1, 1), main="MA Plot AHvU") 


###--------------Get pvals for all comparisons from larger dataset

head(resSpecies)
valsSpecies = cbind(resSpecies$pvalue, resSpecies$padj)
head(valsSpecies)
colnames(valsSpecies)=c("pval.s", "padj.s")
length(valsSpecies[,1])
table(complete.cases(valsSpecies))

head(resBank)
valsBank = cbind(resBank$pvalue, resBank$padj)
head(valsBank)
colnames(valsBank)=c("pval.b", "padj.b")
length(valsBank[,1])
table(complete.cases(valsBank))

head(resALvAH)
valsALvAH = cbind(resALvAH$pvalue, resALvAH$padj)
head(valsALvAH)
colnames(valsALvAH)=c("pval.alah", "padj.alah")
length(valsALvAH[,1])
table(complete.cases(valsALvAH))

head(resALvU)
valsALvU = cbind(resALvU$pvalue, resALvU$padj)
head(valsALvU)
colnames(valsALvU)=c("pval.alu", "padj.alu")
length(valsALvU[,1])
table(complete.cases(valsALvU))

head(resAHvU)
valsAHvU = cbind(resAHvU$pvalue, resAHvU$padj)
head(valsAHvU)
colnames(valsAHvU)=c("pval.ahu", "padj.ahu")
length(valsAHvU[,1])
table(complete.cases(valsAHvU))

#Make rlogdata and pvals table
vsdpvals=cbind(assay(vsd),valsSpecies, valsBank, valsALvAH, valsALvU, valsAHvU)
head(vsdpvals)
dim(vsdpvals)
# 11365    82
table(complete.cases(vsdpvals))
vsdpvals=as.data.frame(vsdpvals)
transcript <- rownames(vsdpvals)
vsdpvals$transcript= transcript #make a new column that is transcript so it can be merged with the genenames file
dim(vsdpvals)

# Load gene name annotations
gg <- read.delim("OfavGenome_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs
vsdpvals=merge(vsdpvals, gg, by="transcript", all=T)
dim(vsdpvals)
head(vsdpvals)
vsdpvals = vsdpvals[!is.na(vsdpvals$FRAAH10),] #removes all the rows with no GE data
dim(vsdpvals)

write.csv(vsdpvals, "fgbdieoff_bm3_wald_VSDandPVALS_geneNames.csv", quote=F)

################################################## Ex situ deoxygenation experiment
# Load data ---------------------------------------------------------------
files <- list.files(path="remapping_May2023/hypoxia_host_counts/host_htseq_counts/")

counts <- read.delim("hypoxia_counts_ofav.tsv")
counts <- counts %>% filter(!grepl("__",gene_id)) %>%
  column_to_rownames("gene_id")
# 48522 isogroups mapped

# make conditions table
conds <- data.frame(
  "sample" = substr(files,1,3),
  "treat" = substr(files,1,1),
  "genet" = substr(files,2,2),
  "rep" = substr(files,3,3)
)

# make sure sample names match up
table(conds$sample == names(counts))
# yep

# Set base mean minimum ---------------------------------------------------
means <- apply(counts,1,mean)
table(means>3)
# FALSE  TRUE 
# 18398 14176  

means3 <- names(means[means>3])
head(means3)
length(means3)
#14176

coral.countFilt <- counts[row.names(counts) %in% means3,]
head(coral.countFilt)

totalCountsFilt <- colSums(coral.countFilt)
totalCountsFilt

min(totalCountsFilt) # 1514583
max(totalCountsFilt) # 4047750
mean(totalCountsFilt) # 2752642

# check sample order
table(names(coral.countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = coral.countFilt,
  colData = conds,
  design = ~ treat + genet)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("treat", "genet"), force=T)

# outliers that fail 2 tests
out2 <- c("H5B")

# remove out2
dim(coral.countFilt) #14176    12
coral.count.out <- coral.countFilt %>% select(-one_of(c(out2)))
dim(coral.count.out) #14176    11

dim(conds)
head(conds)
conds.out <- conds %>% filter(!sample %in% c(out2))
dim(conds.out)

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = coral.count.out,
  colData = conds.out,
  design = ~ treat + genet)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)

#---Results
# log2 fold change (MLE): cond H vs C
res <- results(deds, independentFiltering = F, contrast=c("treat","H","C"))
res

# how many genes pass multiplicity-corrected 0.05 FDR cutoff?
table(res$padj < 0.05)
# FALSE  TRUE 
# 12267  1907

# Extract log-fold changes --------
table(abs(res$log2FoldChange)>1.5)
# FALSE  TRUE 
#  12744  1430 

isogroupsMapped <- nrow(counts)

# how many DEGs as percentage of transcriptome size?
round((table(res$padj<0.05)[2]/isogroupsMapped*100),2) #5.85

# new vsd ------
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# Explore with plots ------------------------------------------------------

#Sample distance heatmap
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.9)(100)
pheatmap(cor(assay(vsd)),border_color=NA, color = heat.colors, main="SampleHeatmap")


###--------------Get pvals
head(res)
vals <- cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval", "padj")
length(vals[,1])
table(complete.cases(vals))

#Make rlogdata and pvals table
vsdpvals <- as.data.frame(cbind(assay(vsd),vals))
head(vsdpvals)
dim(vsdpvals)
# 14176    13
table(complete.cases(vsdpvals))
vsdpvals=as.data.frame(vsdpvals)
transcript <- rownames(vsdpvals)
vsdpvals$transcript= transcript #make a new column that is transcript so it can be merged with the genenames file
dim(vsdpvals)

# Load gene name annotations
gg <- read.delim("OfavGenome_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs
vsdpvals=merge(vsdpvals, gg, by="transcript", all=T)
dim(vsdpvals)
head(vsdpvals)
vsdpvals = vsdpvals[!is.na(vsdpvals$C1B),] #removes all the rows with no GE data
dim(vsdpvals)
 
write.csv(vsdpvals, "hypoxia_bm3_wald_VSDandPVALS_geneNames.csv", quote=F)

# Write results for GO/KOG analysis -------------------------------------------

# by LFC

LFCtreat <- data.frame(cbind("gene"=row.names(res),"LFC"=res$log2FoldChange))
head(LFCtreat)
write.table(LFCtreat,quote=F,row.names=F,file="GO_hypoxia_LFC.csv",sep=",")

# by -log p-value
logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[res$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
#  7399 6237  
logs$logP <- logs$logP*sign
logs$gene <- gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_hypoxia_host_logP.csv",sep=",")


# PCoA ---------

# variance stabilized expression data
exp <- data.frame(assay(vsd))
head(exp)

# condition data
table(conds.out$sam == names(exp))
head(conds.out)

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA 
set.seed(1)
adonisRes <- adonis2(t(exp)~treat+genet,data=conds.out,method="manhattan")
adonisRes
#adonis2(formula = t(exp) ~ treat + genet, data = conds.out, method = "manhattan")
#         Df  SumOfSqs      R2      F Pr(>F)    
#treat     1  91747342 0.22684 3.9334  0.001 ***
#genet     2 149428849 0.36946 3.2031  0.001 ***
#Residual  7 163277543 0.40370                  
#Total    10 404453735 1.00000                  


# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors

# plotting PCoA
margin <- 1

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA
pdf(file="hypoxia_pca_host.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
     main = "Host Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
ordiellipse(scores,conds.out$treat,label=F, draw="polygon",
             col = c("orange", "purple4"), alpha = 0.5, border = F)
ordispider(scores,conds.out$treat,label=F)
points(scores[grep("C1",conds.out$sam),xaxis],
       scores[grep("C1",conds.out$sam),yaxis], col="orange", pch=15) +
  points(scores[grep("C2",conds.out$sam),xaxis],
         scores[grep("C2",conds.out$sam),yaxis], col="orange", pch=16) +
  points(scores[grep("C5",conds.out$sam),xaxis],
         scores[grep("C5",conds.out$sam),yaxis], col="orange", pch=17) +
  points(scores[grep("H1",conds.out$sam),xaxis],
         scores[grep("H1",conds.out$sam),yaxis], col="purple4", pch=15) +
  points(scores[grep("H2",conds.out$sam),xaxis],
         scores[grep("H2",conds.out$sam),yaxis], col="purple4", pch=16) +
  points(scores[grep("H5",conds.out$sam),xaxis],
         scores[grep("H5",conds.out$sam),yaxis], col="purple4", pch=17)

# legend of sites 
legend("bottomleft", 
       c("Control","Hypoxia"),
       pch=c(19,19), 
       col=c("orange","purple4"), cex=1.5, bty = "n")

legend("bottomright", 
       c("Genet 1","Genet 2", "Genet 5"),
       pch=c(15,16,17), 
       col=c("black","black","black"), cex=1.5, bty = "n")

#insert p value 
legend("topleft", inset=.02, 
       paste("Condition p=",adonisRes$`Pr(>F)`[1], sep=" "), cex=1.5, bty='n') 
legend("topright", inset=.02, 
       paste("Genet p=",adonisRes$`Pr(>F)`[2], sep=" "), cex=1.5, bty='n') 
dev.off()

###################### make venn diagram with the two species responses
OFvsdpvals <- read.csv("fgbdieoff_bm3_wald_VSDandPVALS_OFhost_geneNames.csv")
row.names(OFvsdpvals) <- OFvsdpvals$transcript
OFvsdpvals$X <- NULL
OFvsdpvals <- as.data.frame(OFvsdpvals)
head(OFvsdpvals)

FRvsdpvals <- read.csv("fgbdieoff_bm3_wald_VSDandPVALS_FRhost_geneNames.csv")
row.names(FRvsdpvals) <- FRvsdpvals$transcript
FRvsdpvals $X <- NULL
FRvsdpvals <- as.data.frame(FRvsdpvals)
head(FRvsdpvals)

hypovsdpvals <- read.csv("hypoxia_bm3_wald_VSDandPVALS_geneNames.csv")
row.names(hypovsdpvals) <- hypovsdpvals$transcript
hypovsdpvals $X <- NULL
hypovsdpvals <- as.data.frame(hypovsdpvals)
head(hypovsdpvals)

OF <- row.names(OFvsdpvals[OFvsdpvals $padj.alu<0.05 & !is.na(OFvsdpvals $padj.alu) | OFvsdpvals$padj.alah<0.05 & !is.na(OFvsdpvals $padj.alah),])
FR <- row.names(FRvsdpvals[FRvsdpvals $padj.alu<0.05 & !is.na(FRvsdpvals $padj.alu) | FRvsdpvals $padj.alah<0.05 & !is.na(FRvsdpvals $padj.alah),])
hyo <- row.names(hypovsdpvals[hypovsdpvals$padj<0.05 & !is.na(hypovsdpvals $padj),])

candidates <- list("O.faveolata LME"= OF, "O.franksi LME"= FR, "O.faveolata Hypoxia"= hyo)
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("red","orange","blue"),
  alpha = 0.3,
  label.col = rep("black", 7),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred","darkorange","darkblue"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.08),
  cat.pos = 1
);
quartz()
grid.draw(prettyvenn)

overlap <- calculate.overlap(candidates)

#$a5
#[1] "XM_020748780.1" "XM_020749191.1" "XM_020751503.1" "XM_020752970.1" "XM_020758673.1" "XM_020761566.1" "XM_020762325.1" "XM_020765809.1"

############ plot the intersection genes 
vsdpvals <- read.csv("MarieReanalysis2022/HostGE_REANALYSIS/REANALYSIS_fgbdieoff_bm3_wald_VSDandPVALS_geneNames_MS.csv") #none of these genes have annotations (but check with actual genome reference)
vsdpvals$X <- NULL
vsdpvals <- as.data.frame(vsdpvals)
head(vsdpvals)
vsdpvals <- vsdpvals[ -c(84) ]

Cands1= vsdpvals[grepl("XM_020748780.1",vsdpvals$transcript), ]
Cands2= vsdpvals[grepl("XM_020749191.1",vsdpvals$transcript), ]
Cands3= vsdpvals[grepl("XM_020751503.1",vsdpvals$transcript), ]
Cands4= vsdpvals[grepl("XM_020752970.1",vsdpvals$transcript), ]
Cands5= vsdpvals[grepl("XM_020758673.1",vsdpvals$transcript), ]
Cands6= vsdpvals[grepl("XM_020761566.1",vsdpvals$transcript), ]
Cands7= vsdpvals[grepl("XM_020762325.1",vsdpvals$transcript), ]
Cands8= vsdpvals[grepl("XM_020765809.1",vsdpvals$transcript), ]

Cands=rbind(Cands1, Cands2, Cands3, Cands4, Cands5, Cands6, Cands7, Cands8)
names(Cands)[names(Cands) == 'transcript'] <- 'transcr'

exp <- Cands %>% select(-grep("p",colnames(Cands))) # get rid of columns with "p"
head(exp)

geneL= exp %>%
	gather("sample", "expr", -transcr)
geneL$species = str_sub(geneL$sample, 1,2)
geneL$lesion <- ifelse( grepl("AH", geneL$sample), "AH", 
                  ifelse( grepl("AL", geneL$sample), "AL", "U"))
geneL$bank <- ifelse( grepl("W",geneL$sample), "west", "east")

geneL$lesion_f = factor(geneL$lesion, levels=c("U","AH","AL"))

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("transcr","species","lesion_f"))

pd <- position_dodge(0.1)
LMEgenes=ggplot(summ,aes(x= lesion_f,y=expr, color=species))+
	geom_point(aes(shape= species),size=3,position=pd)+
	geom_line(aes(group= species,linetype= species),position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
	scale_shape_manual(values=c(16,17))+
	scale_color_manual(values=c("darkorange","darkred"))+
	theme_minimal()+
	facet_wrap(~ transcr,scales="free_y", ncol=8)+
	theme(legend.text=element_text(size=10)) +
	theme(legend.key = element_blank())+
	theme(legend.direction = 'horizontal', legend.position = 'top')


### do the hypoxia ones
vsdpvals <- read.csv("hypoxia_bm3_wald_VSDandPVALS_geneNames.csv")
vsdpvals$X <- NULL
vsdpvals <- as.data.frame(vsdpvals)
head(vsdpvals)

vsdpvals <- vsdpvals[ -c(15) ]

Cands1= vsdpvals[grepl("XM_020748780.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata dimethylglycine dehydrogenase, mitochondrial-like (LOC110043332), no abbreviated annotation
Cands2= vsdpvals[grepl("XM_020749191.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata prostatic spermine-binding protein-like (LOC110043704), no abbreviated annotation
Cands3= vsdpvals[grepl("XM_020751503.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata uncharacterized LOC110045874 (LOC110045874), no abbreviated annotation
Cands4= vsdpvals[grepl("XM_020752970.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata uncharacterized LOC110047223 (LOC110047223), no abbreviated annotation
Cands5= vsdpvals[grepl("XM_020758673.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata galaxin-like (LOC110052545), no abbreviated annotation
Cands6= vsdpvals[grepl("XM_020761566.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata glutamic acid-rich protein-like (LOC110055197), no abbreviated annotation
Cands7= vsdpvals[grepl("XM_020762325.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata pancreatic secretory granule membrane major glycoprotein GP2-like (LOC110055901), no abbreviated annotation
Cands8= vsdpvals[grepl("XM_020765809.1",vsdpvals$transcript), ] #PREDICTED: Orbicella faveolata calumenin-like (LOC110059122), , no abbreviated annotation

Cands=rbind(Cands1, Cands2, Cands3, Cands4, Cands5, Cands6, Cands7, Cands8)
names(Cands)[names(Cands) == 'transcript'] <- 'transcr'

exp <- Cands %>% select(-grep("p",colnames(Cands))) # get rid of columns with "p"
head(exp)

geneL= exp %>%
	gather("sample", "expr", -transcr)
geneL$treat = str_sub(geneL$sample, 1,1)

geneL$treat_f = factor(geneL$treat, levels=c("C","H"))

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("transcr","treat_f"))
summ$paired=(c(rep("1",2),rep("2",2),rep("3",2),rep("4",2),rep("5",2),rep("6",2),rep("7",2),rep("8",2)))


pd <- position_dodge(0.1)
HPXgenes=ggplot(summ,aes(x= treat_f,y=expr))+
	geom_point(shape=17, color="darkblue",fill="darkblue",size=3,position=pd)+
	geom_line(aes(group = paired),color="darkblue", position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd, color="darkblue")+
	theme_minimal()+
	facet_wrap(~ transcr,scales="free_y", ncol=8)+
	theme(legend.text=element_text(size=10)) +
	theme(legend.key = element_blank())+
	theme(legend.direction = 'horizontal', legend.position = 'top')

fig3 <- plot_grid(LMEgenes, HPXgenes, ncol = 1,rel_heights = c(1.2, 1))

####################################################################################### PCoA ---------
dfs=load("resultsLME.Rdata")

# variance stabilized expression data
exp <- data.frame(assay(vsd))
head(exp)

# condition data
table(conds.out$samID == names(exp))
head(conds.out)

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)

adonisRes <- adonis2(div.dd.veg~species+bank+lesion,
                     data=conds.out,
                     sim.function="vegdist",
                     method="manhattan",by="margin")
adonisRes
#adonis2(formula = div.dd.veg ~ species + bank + lesion, data = conds.out, method = "manhattan", by = "margin", sim.function = "vegdist")
#         Df SumOfSqs      R2      F Pr(>F)    
#species   1    83.88 0.05633 4.3544  0.001 ***
#bank      1    27.42 0.01841 1.4233  0.028 *  
#lesion    2    65.70 0.04412 1.7052  0.001 ***
#Residual 67  1290.62 0.86677                  
#Total    71  1489.00 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors
head(scores[,1])

# plotting PCoA----
margin <- 1

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA plot
quartz()
pdf(file="lme_pca_host.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Host Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
ordiellipse(scores,conds.out$lesion,label=F, draw="polygon",
             col = c("orange", "red", "blue"), alpha = 0.5, border = F)
ordispider(scores, conds.out $species)
# plot affected fav and frank
points(scores[conds.out $lesion=="AL" & conds.out $species=="fav",xaxis],
       scores[conds.out $lesion=="AL" & conds.out $species=="fav",yaxis],
       col="red", pch=17, cex=1) +
points(scores[conds.out $lesion=="AL" & conds.out $species=="frank",xaxis],
         scores[conds.out $lesion=="AL" & conds.out $species=="frank",yaxis],
         col="red", pch=16, cex=1) +
# plot affected-healthy fav and frank  
points(scores[conds.out $lesion=="AH" & conds.out $species=="fav",xaxis],
         scores[conds.out $lesion=="AH" & conds.out $species=="fav",yaxis],
         col="orange", pch=17, cex=1) +
points(scores[conds.out $lesion=="AH" & conds.out $species=="frank",xaxis],
         scores[conds.out $lesion=="AH" & conds.out $species=="frank",yaxis],
         col="orange", pch=16, cex=1) +
# plot affected-healthy fav and frank  
points(scores[conds.out $lesion=="U" & conds.out $species=="fav",xaxis],
         scores[conds.out $lesion=="U" & conds.out $species=="fav",yaxis],
         col="blue", pch=17, cex=1) +
points(scores[conds.out $lesion=="U" & conds.out $species=="frank",xaxis],
         scores[conds.out $lesion=="U" & conds.out $species=="frank",yaxis],
         col="blue", pch=16, cex=1)
# legend of sites 
legend("topright", 
       c("Affected-Lesion", "Affected-Healthy", "Unaffected"),
       pch=c(8,8), 
       col=c("red","orange","blue"), cex=1.5, bty = "n")
legend("bottomright", 
       c(expression(italic("O. faveolata")),expression(italic("O. franksi"))),
       pch=c(17,16), 
       col=c("black","black"), cex=1.5, bty = "n")

#insert p value 
legend("topleft",
       paste("Host Species p=",adonisRes$`Pr(>F)`[1], sep=" "), 
       cex=1.5, bty='n')  

legend("bottomleft", 
       paste("Lesion p = ",adonisRes$`Pr(>F)`[3], sep=" "), 
       cex=1.5, bty='n')

legend("bottomleft", 
       paste("Colony p = ",adonisRes$`Pr(>F)`[4], sep=" "), 
       cex=1.5, bty='n')  
dev.off()


#### Plot LME PCA separated by species
#Start with Ofav
OF=load("REANALYSIS_resultsLME_OF.Rdata")

# variance stabilized expression data
exp <- data.frame(assay(OFvsd))
head(exp)

# condition data
table(OFconds$sam == names(exp))
head(OFconds)

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis2(div.dd.veg~bank+lesion+genet,
                     data= OFconds,
                     sim.function="vegdist",
                     method="manhattan",by="margin")
adonisRes
#adonis2(formula = div.dd.veg ~ bank + lesion + genet, data = OFconds, method = "manhattan", by = "margin", sim.function = "vegdist")
#         Df SumOfSqs      R2      F Pr(>F)   
#bank      0     0.00 0.00000    Inf          
#lesion    2    40.51 0.07484 1.4006  0.008 **
#genet    16   270.98 0.50066 1.1711  0.018 * 
#Residual 13   188.00 0.34734                 
#Total    32   541.26 1.00000                 


# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors
head(scores[,1])

# plotting PCoA----
margin <- 0.75

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA plot
quartz()
pdf(file="lme_pca_host_OF.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Host Gene Expression in Orbicella faveolata",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
ordiellipse(scores, OFconds $lesion,label=F, draw="polygon",
             col = c("orange", "red", "blue"), alpha = 0.5, border = F)
ordispider(scores, OFconds$lesion)
# plot affected fav 
points(scores[OFconds$lesion=="AL",xaxis],
		scores[OFconds $lesion=="AL",yaxis],
       col="red", pch=17, cex=1) +
# plot affected-healthy fav 
points(scores[OFconds$lesion=="AH",xaxis],
		scores[OFconds $lesion=="AH",yaxis],
         col="orange", pch= 17, cex=1) +
# plot affected-healthy fav 
points(scores[OFconds $lesion=="U",xaxis],
		scores[OFconds $lesion=="U",yaxis],
         col="blue", pch= 17, cex=1)
# legend of sites 
legend("topleft", 
       c("Affected-Lesion", "Affected-Healthy", "Unaffected"),
       pch=c(8,8), 
       col=c("red","orange","blue"), cex=1.5, bty = "n")

#insert p value 
legend("bottomleft", 
       paste("Lesion p = ",adonisRes$`Pr(>F)`[2], sep=" "), 
       cex=1.5, bty='n')
legend("bottomright", 
       paste("Colony p = ",adonisRes$`Pr(>F)`[3], sep=" "), 
       cex=1.5, bty='n')  
dev.off()


##### PCoA for Orbicella franksi
FR=load("resultsLME_FR.Rdata")

exp <- data.frame(assay(FRvsd))
head(exp)

# condition data
table(FRconds$sam == names(exp))
head(FRconds)

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis2(div.dd.veg~bank+lesion+genet,
                     data= FRconds,
                     sim.function="vegdist",
                     method="manhattan",by="margin")
adonisRes
#adonis2(formula = div.dd.veg ~ bank + lesion + genet, data = FRconds, method = "manhattan", by = "margin", sim.function = "vegdist")
#         Df SumOfSqs      R2      F Pr(>F)  
#bank      0     0.00 0.00000    Inf         
#lesion    2    54.05 0.06602 1.3703  0.021 *
#genet    18   397.68 0.48572 1.1202  0.085 .
#Residual 17   335.28 0.40951                
#Total    38   818.73 1.00000                
                     
# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors
head(scores[,1])

# plotting PCoA----
margin <- 0.75

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA plot
quartz()
pdf(file="lme_pca_host_FR.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Host Gene Expression in Orbicella franksi",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
ordiellipse(scores, FRconds$lesion,label=F, draw="polygon",
             col = c("orange", "red", "blue"), alpha = 0.5, border = F)
ordispider(scores, FRconds $lesion)
# plot affected fav 
points(scores[FRconds $lesion=="AL",xaxis],
		scores[FRconds $lesion=="AL",yaxis],
       col="red", pch=16, cex=1) +
# plot affected-healthy fav 
points(scores[FRconds $lesion=="AH",xaxis],
		scores[FRconds $lesion=="AH",yaxis],
         col="orange", pch=16, cex=1) +
# plot affected-healthy fav 
points(scores[FRconds $lesion=="U",xaxis],
		scores[FRconds $lesion=="U",yaxis],
         col="blue", pch=16, cex=1)
# legend of sites 
legend("topleft", 
       c("Affected-Lesion", "Affected-Healthy", "Unaffected"),
       pch=c(8,8), 
       col=c("red","orange","blue"), cex=1.5, bty = "n")

#insert p value 
legend("bottomleft", 
       paste("Lesion p = ",adonisRes$`Pr(>F)`[2], sep=" "), 
       cex=1.5, bty='n')
legend("bottomright", 
       paste("Colony p = ",adonisRes$`Pr(>F)`[3], sep=" "), 
       cex=1.5, bty='n')  
dev.off()


################ plot hypoxia genes from Alderice 2020
vsdpvals <- read.csv("fgbdieoff_bm3_wald_VSDandPVALS_geneNames.csv")
vsdpvals$X <- NULL
vsdpvals <- as.data.frame(vsdpvals)
head(vsdpvals)

vsdpvals <- vsdpvals[ -c(84) ]

#HIFs
#Cands1= vsdpvals[grepl("XM_020760929.1",vsdpvals$transcript), ] #not present in dataset (must be lowly expressed)
#Cands2= vsdpvals[grepl("XM_020760974.1",vsdpvals$transcript), ] #not present in dataset (must be lowly expressed)

Cands3= vsdpvals[grepl("XM_020774709.1",vsdpvals$transcript), ] #padj.alah 0.2859135, HIF1a, padj.alu 0.08474852

#HSP90s
#Cands4= vsdpvals[grepl("XM_020751073.1",vsdpvals$transcript), ]#not present in dataset (must be lowly expressed)
#Cands5= vsdpvals[grepl("XM_020751076.1",vsdpvals$transcript), ]#not present in dataset (must be lowly expressed)
#Cands6= vsdpvals[grepl("XM_020751075.1",vsdpvals$transcript), ]#not present in dataset (must be lowly expressed)

#plot one of the top DEGs for the two species. 
Cands3= vsdpvals[grepl("XM_020775010.1",vsdpvals$transcript), ] #padj.alah, BHMT
Cands3= vsdpvals[grepl("XM_020774546.1",vsdpvals$transcript), ] #padj.alah, PRDX6
Cands3= vsdpvals[grepl("XM_020765263.1",vsdpvals$transcript), ] #padj.alah, CHL1
Cands3= vsdpvals[grepl("XM_020748813.1",vsdpvals$transcript), ] #padj.alah, COMP

names(Cands3)[names(Cands3) == 'transcript'] <- 'transcr'

exp <- Cands3 %>% select(-grep("p",colnames(Cands3))) # get rid of columns with "p"
head(exp)

geneL= exp %>%
	gather("sample", "expr", -transcr)
geneL$species = str_sub(geneL$sample, 1,2)
geneL$lesion <- ifelse( grepl("AH", geneL$sample), "AH", 
                  ifelse( grepl("AL", geneL$sample), "AL", "U"))
geneL$bank <- ifelse( grepl("W",geneL$sample), "west", "east")

geneL$lesion_f = factor(geneL$lesion, levels=c("U","AH","AL"))

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("transcr","species","lesion_f"))

pd <- position_dodge(0.1)
LMEgenes=ggplot(summ,aes(x= lesion_f,y=expr, color=species))+
	geom_point(aes(shape=species),size=3,position=pd)+
	geom_line(aes(group= species,linetype= species),position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
	scale_shape_manual(values=c(16,17))+
	scale_color_manual(values=c("darkorange","darkred"))+
	theme_minimal()+
	theme(legend.text=element_text(size=10)) +
	theme(legend.key = element_blank())+
	theme(legend.direction = 'horizontal', legend.position = 'top')+
	annotate("text",  x=Inf, y = Inf, label = "padj ahal = 0.29, padj alu = 0.084", vjust=1, hjust=1)

################ hypoxia dataset
vsdpvals <- read.csv("hypoxia_bm3_wald_VSDandPVALS_geneNames.csv")
vsdpvals$X <- NULL
vsdpvals <- as.data.frame(vsdpvals)
head(vsdpvals)

#HIFs
#Cands1= vsdpvals[grepl("XM_020760929.1",vsdpvals$transcript), ] #not present in dataset (must be lowly expressed)
Cands2= vsdpvals[grepl("XM_020760974.1",vsdpvals$transcript), ] #padj 0.2673318, HIF1AN
Cands3= vsdpvals[grepl("XM_020774709.1",vsdpvals$transcript), ] #padj 1.975848e-08, HIF1A

#HSP90s
#Cands4= vsdpvals[grepl("XM_020751073.1",vsdpvals$transcript), ]#not present in dataset (must be lowly expressed)
#Cands5= vsdpvals[grepl("XM_020751076.1",vsdpvals$transcript), ]#not present in dataset (must be lowly expressed)
#Cands6= vsdpvals[grepl("XM_020751075.1",vsdpvals$transcript), ]#not present in dataset (must be lowly expressed)

#Cands=rbind(Cands2,Cands3)
Cands <- Cands3[ -c(15) ]
names(Cands)[names(Cands) == 'transcript'] <- 'transcr'

exp <- Cands %>% select(-grep("p",colnames(Cands))) # get rid of columns with "p"
head(exp)

geneL= exp %>%
	gather("sample", "expr", -transcr)
geneL$treat = str_sub(geneL$sample, 1,1)

geneL$treat_f = factor(geneL$treat, levels=c("C","H"))

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("transcr","treat_f"))
#summ$paired=(c(rep("1",2),rep("2",2)))
summ$paired=(c(rep("1",2)))

pd <- position_dodge(0.1)
HPXgenes=ggplot(summ,aes(x= treat_f,y=expr))+
	geom_point(shape= 17, color="purple4",fill="purple4",size=3,position=pd)+
	geom_line(aes(group = paired),color="purple4", position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd, color="purple4")+
	#scale_shape_manual(values=c(17))+
	theme_minimal()+
	annotate("text",  x=Inf, y = Inf, label = "padj = 3.9e-12", vjust=1, hjust=1)

plot_grid(LMEgenes, HPXgenes, ncol = 1,rel_heights = c(1.2, 1))



