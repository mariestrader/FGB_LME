library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("pheatmap") # for sample heatmap
library("VennDiagram") # for Venn diagram of #s of DEGs
library("tidyverse") # for wrangling, plotting
library("vegan") # for PCoA
library("ape") # for PCoA

# Load count data from LME---------------------------------------------------------------
counts <- read.delim("LME_FGB_OfavB1.txt")

counts <- counts %>%
  select(-X.1) %>%
  column_to_rownames("X") %>%
  rename_all(~str_replace(., "_.*", ""))

names(counts)
head(counts) 
length(counts[,1])
# 39010 isogroups mapped

# make conditions table --------
conds <- data.frame(samID = names(counts))

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

table(conds$samID == names(counts))


# Symbiont ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts = colSums(counts)
head(totalCounts)
totalCounts
min(totalCounts) #37272
mean(totalCounts) #290424.7
max(totalCounts)  #758279

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conds,
  design = ~ species + lesion + bank)

# Set base mean minimum ---------------------------------------------------
means <- apply(counts,1,mean)
table(means>3)
# FALSE  TRUE 
# 26296 12714 

means3 <- names(means[means>3])
head(means3)
length(means3)
# 12714

sym.countFilt <- counts[row.names(counts) %in% means3,]
head(sym.countFilt)

totalCountsFilt <- colSums(sym.countFilt)
totalCountsFilt

min(totalCountsFilt) #35802
max(totalCountsFilt) #717757
mean(totalCountsFilt) #274134.9

# check sample order
table(names(sym.countFilt) == as.vector(conds$sam))
# TRUE 76

# Reconstruct data object (filtered) ------------------------------
ddsFilt <- DESeqDataSetFromMatrix(
  countData = sym.countFilt,
  colData = conds,
  design = ~ species + lesion + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("species","lesion","bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out2 <- c("FRAAL4", "FRAUE6", "OFAH4")
#out1 <- c("OFUW8")

# remove out2
dim(sym.countFilt)
sym.count.out <- sym.countFilt %>% select(-one_of(c(out2)))
dim(sym.count.out)
# 12714    73

dim(conds)
head(conds)
conds.out <- conds %>% filter(!samID %in% c(out2))
dim(conds.out)
# 73 5

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = sym.count.out,
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
# 49

# log2 fold change (MLE): bank west vs east 
resBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
resBank
table(resBank$padj<0.05)
# 21

# log2 fold change (MLE): lesion AL vs AH 
resALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
resALvAH
table(resALvAH$padj<0.05)
# 26

#log2 fold change (MLE): lesion AL vs U 
resALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
resALvU
table(resALvU$padj<0.05)
# 36

#log2 fold change (MLE): lesion AH vs U 
resAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
resAHvU
table(resAHvU$padj<0.05)
# 0

# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# 30

table(abs(resBank$log2FoldChange)>1.5)
# 43

table(abs(resALvAH$log2FoldChange)>1.5)
# 48

table(abs(resAHvU$log2FoldChange)>1.5)
# 15

table(abs(resALvU$log2FoldChange)>1.5)
# 59

# differential expression as percentage of # of reads ----

isogroupsMapped <- nrow(counts)

# how many DEGs as percentage of transcriptome size?
round((table(resBank$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.05
round((table(resSpecies$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.13
round((table(resALvAH$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.07
round((table(resALvU$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.09
round((table(resAHvU$padj<0.05)[2]/isogroupsMapped*100),2)
# NA (no DEGs) 

# Load and save ------
save(vsd, conds.out, sym.count.out, ddsFiltOut, deds, 
     resSpecies, resBank, resAHvU, resALvAH, resALvU, file="resultsLME_sym_MS.Rdata")
load("resultsLME_sym_MS.Rdata") 

###--------------Get pvals
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
# 12714    83
table(complete.cases(vsdpvals))
rownames(vsdpvals)<-sub("sym","", rownames(vsdpvals))
vsdpvals=as.data.frame(vsdpvals)
transcript <- rownames(vsdpvals)
vsdpvals$transcript= transcript #make a new column that is transcript so it can be merged with the genenames file

# Load gene name annotations
gg <- read.delim("Bmin_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs
vsdpvals=merge(vsdpvals, gg, by="transcript", all=T)
dim(vsdpvals) #51199    85
head(vsdpvals)
vsdpvals = vsdpvals[!is.na(vsdpvals$FRAAH10),] #removes all the rows with no GE data
dim(vsdpvals) #12714    85

write.csv(vsdpvals, "fgbdieoff_bm3_wald_VSDandPVALS_SymgeneNames.csv", quote=F)

# Write results for GO/KOG analysis -------------------------------------------
load("resultsLME_sym_MS.Rdata")

# by LFC
rownames(resSpecies)<-sub("sym","", rownames(resSpecies))
LFCspecies <- data.frame(cbind("gene"=row.names(resSpecies),"LFC"=resSpecies$log2FoldChange))
write.table(LFCspecies,quote=F,row.names=F,file="GO_Species_LFC.csv",sep=",")

rownames(resBank)<-sub("sym","", rownames(resBank))
LFCbank <- data.frame(cbind("gene"=row.names(resBank),"LFC"=resBank$log2FoldChange))
write.table(LFCbank,quote=F,row.names=F,file="GO_Bank_LFC.csv",sep=",")

rownames(resAHvU)<-sub("sym","", rownames(resAHvU))
LFCAHvU <- data.frame(cbind("gene"=row.names(resAHvU),"LFC"=resAHvU$log2FoldChange))
write.table(LFCAHvU,quote=F,row.names=F,file="GO_AHvU_LFC.csv",sep=",")

rownames(resALvAH)<-sub("sym","", rownames(resALvAH))
LFCALvAH <- data.frame(cbind("gene"=row.names(resALvAH),"LFC"=resALvAH$log2FoldChange))
write.table(LFCALvAH,quote=F,row.names=F,file="GO_ALvAH_LFC.csv",sep=",")

rownames(resALvU)<-sub("sym","", rownames(resALvU))
LFCALvU <- data.frame(cbind("gene"=row.names(resALvU),"LFC"=resALvU$log2FoldChange))
write.table(LFCALvU,quote=F,row.names=F,file="GO_ALvU_LFC.csv",sep=",")

# by -log p-value
logs=data.frame(cbind("gene"=row.names(resSpecies),"logP"=round(-log(resSpecies$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resSpecies$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 2476 2460 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Species_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resBank),"logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 2633 2303   
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Bank_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resAHvU),"logP"=round(-log(resAHvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resAHvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 2480 2456
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_AHvU_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resALvAH),"logP"=round(-log(resALvAH$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resALvAH$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 2527 2409   
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_ALvAH_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resALvU),"logP"=round(-log(resALvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resALvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 2564 2372 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_ALvU_logP.csv",sep=",")

################################### PCoA ---------
load("resultsLME_sym_MS.Rdata")

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
adonisRes <- adonis2(div.dd.veg~species+bank+lesion,
                     data=conds.out,
                     sim.function="vegdist",
                     method="manhattan",by="margin")
adonisRes
#adonis2(formula = div.dd.veg ~ species + bank + lesion, data = conds.out, method = "manhattan", by = "margin", sim.function = "vegdist")
#         Df SumOfSqs      R2      F Pr(>F)    
#species   1    28.00 0.02278 1.6716  0.001 ***
#bank      1    19.77 0.01609 1.1805  0.062 .  
#lesion    2    39.57 0.03219 1.1812  0.021 *  
#Residual 68  1138.85 0.92663                  
#Total    72  1229.02 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors
head(scores[,1])

# plotting PCoA----
margin <- 0.5

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA plot
#quartz()
pdf(file="lme_pca_sym_genet.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Breviolum minutum Gene Expression",
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
       col="red", pch=2, cex=1) +
points(scores[conds.out $lesion=="AL" & conds.out $species=="frank",xaxis],
         scores[conds.out $lesion=="AL" & conds.out $species=="frank",yaxis],
         col="red", pch=1, cex=1) +
# plot affected-healthy fav and frank  
points(scores[conds.out $lesion=="AH" & conds.out $species=="fav",xaxis],
         scores[conds.out $lesion=="AH" & conds.out $species=="fav",yaxis],
         col="orange", pch=2, cex=1) +
points(scores[conds.out $lesion=="AH" & conds.out $species=="frank",xaxis],
         scores[conds.out $lesion=="AH" & conds.out $species=="frank",yaxis],
         col="orange", pch=1, cex=1) +
# plot affected-healthy fav and frank  
points(scores[conds.out $lesion=="U" & conds.out $species=="fav",xaxis],
         scores[conds.out $lesion=="U" & conds.out $species=="fav",yaxis],
         col="blue", pch=2, cex=1) +
points(scores[conds.out $lesion=="U" & conds.out $species=="frank",xaxis],
         scores[conds.out $lesion=="U" & conds.out $species=="frank",yaxis],
         col="blue", pch=1, cex=1)
# legend of sites 
legend("topright", 
       c("Affected-Lesion", "Affected-Healthy", "Unaffected"),
       pch=c(8,8), 
       col=c("red","orange","blue"), cex=1.5, bty = "n")
legend("bottomright", 
       c(expression(italic("O. faveolata")),expression(italic("O. franksi"))),
       pch=c(2,1), 
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

# Subset data by species ------------------------------------------------------
load("resultsLME_sym_MS.Rdata")

OFconds= conds.out[conds.out$species== "fav",]
dim(OFconds) #35 5
FRconds= conds.out[conds.out$species== "frank",]
dim(FRconds) #38  5

OFsym.count.out <- sym.count.out[, OFconds$sam]
all((OFconds$sam) == colnames(OFsym.count.out))
dim(OFsym.count.out) #12714    35

FRsym.count.out <- sym.count.out[, FRconds$sam]
all((FRconds$sam) == colnames(FRsym.count.out))
dim(FRsym.count.out) # 12714    38


# Run DESeq for OF ------------------------------

OFddsFiltOut <- DESeqDataSetFromMatrix(
  countData = OFsym.count.out,
  colData = OFconds,
  design = ~ lesion + bank)

# remake vsd
OFvsd <- varianceStabilizingTransformation(OFddsFiltOut, blind=TRUE)

deds <- DESeq(OFddsFiltOut)

# log2 fold change (MLE): bank west vs east 
OFresBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
OFresBank
table(OFresBank$padj<0.05)
# 11

# log2 fold change (MLE): lesion AL vs AH 
OFresALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
OFresALvAH
table(OFresALvAH$padj<0.05)
# 0

#log2 fold change (MLE): lesion AL vs U 
OFresALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
OFresALvU
table(OFresALvU$padj<0.05)
# 12

#log2 fold change (MLE): lesion AH vs U 
OFresAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
OFresAHvU
table(OFresAHvU$padj<0.05)
# 2

save(OFvsd, OFconds, OFsym.count.out, OFddsFiltOut, deds, 
     OFresBank, OFresALvAH, OFresALvU, OFresAHvU, file="resultsLME_OF_syms.Rdata")

vals= cbind(OFresBank$pvalue, OFresBank$padj, OFresALvAH$pvalue, OFresALvAH$padj, OFresALvU$pvalue, OFresALvU$padj, OFresAHvU$pvalue, OFresAHvU$padj)
head(vals)
colnames(vals)=c("pval.b", "padj.b","pval.alah", "padj.alah","pval.alu", "padj.alu","pval.ahu", "padj.ahu")
length(vals[,1])
table(complete.cases(vals))
OFvsdpvals=cbind(assay(OFvsd), vals)
head(OFvsdpvals)
dim(OFvsdpvals)
# 12714    43
table(complete.cases(OFvsdpvals))

#remove sym from rownames
rownames(OFvsdpvals)<-sub("sym","", rownames(OFvsdpvals))
OFvsdpvals =as.data.frame(OFvsdpvals)
transcript <- rownames(OFvsdpvals)
OFvsdpvals $transcript= transcript #make a new column that is transcript so it can be merged with the genenames file
dim(OFvsdpvals)

# Load gene name annotations
gg <- read.delim("Bmin_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs

OFvsdpvals =merge(OFvsdpvals, gg, by="transcript", all=T)
dim(OFvsdpvals)
head(OFvsdpvals)
OFvsdpvals = OFvsdpvals[!is.na(OFvsdpvals$OFAH10),] #removes all the rows with no GE data
dim(vsdpvals) #12714    83
head(OFvsdpvals)

write.csv(OFvsdpvals, "fgbdieoff_bm3_wald_VSDandPVALS_OFsym_geneNames.csv", quote=F)

##########################################output for GO and KOG
#AL v. AH. for Orbicella faveolata
rownames(OFresALvAH)<-sub("sym","", rownames(OFresALvAH))
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
rownames(OFresALvU)<-sub("sym","", rownames(OFresALvU))
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
rownames(OFresAHvU)<-sub("sym","", rownames(OFresAHvU))
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

# Run DESeq for FR ------------------------------

FRddsFiltOut <- DESeqDataSetFromMatrix(
  countData = FRsym.count.out,
  colData = FRconds,
  design = ~ lesion + bank)

# remake vsd
FRvsd <- varianceStabilizingTransformation(FRddsFiltOut, blind=TRUE)

deds <- DESeq(FRddsFiltOut)

# log2 fold change (MLE): bank west vs east 
FRresBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
FRresBank
table(FRresBank$padj<0.05)
# 5

# log2 fold change (MLE): lesion AL vs AH 
FRresALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
FRresALvAH
table(FRresALvAH$padj<0.05)
# 11

#log2 fold change (MLE): lesion AL vs U 
FRresALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
FRresALvU
table(FRresALvU$padj<0.05)
# 10

#log2 fold change (MLE): lesion AH vs U 
FRresAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
FRresAHvU
table(FRresAHvU$padj<0.05)
# 0

save(FRvsd, FRconds, FRsym.count.out, FRddsFiltOut, deds, 
     FRresBank, FRresALvAH, FRresALvU, FRresAHvU, file="resultsLME_FR_syms.Rdata")

vals= cbind(FRresBank$pvalue, FRresBank$padj, FRresALvAH$pvalue, FRresALvAH$padj, FRresALvU$pvalue, FRresALvU$padj, FRresAHvU$pvalue, FRresAHvU$padj)
head(vals)
colnames(vals)=c("pval.b", "padj.b","pval.alah", "padj.alah","pval.alu", "padj.alu","pval.ahu", "padj.ahu")
length(vals[,1])
table(complete.cases(vals))
FRvsdpvals=cbind(assay(FRvsd), vals)
head(FRvsdpvals)
dim(FRvsdpvals)
# 12714    46
table(complete.cases(FRvsdpvals))

#remove sym from rownames
rownames(FRvsdpvals)<-sub("sym","", rownames(FRvsdpvals))
FRvsdpvals =as.data.frame(FRvsdpvals)
transcript <- rownames(FRvsdpvals)
FRvsdpvals $transcript= transcript #make a new column that is transcript so it can be merged with the genenames file
dim(FRvsdpvals)

# Load gene name annotations
gg <- read.delim("Bmin_iso2gene.tab", header=F,na.strings=c("","NA"))
head(gg)
dim(gg)
colnames(gg)=c("transcript","gene")

#combine the vsdpvals with the gene names file, keep all instances in both dfs
FRvsdpvals =merge(FRvsdpvals, gg, by="transcript", all=T)
dim(FRvsdpvals)
head(FRvsdpvals)
FRvsdpvals = FRvsdpvals[!is.na(FRvsdpvals$FRAAH10),] #removes all the rows with no GE data
dim(FRvsdpvals)
head(FRvsdpvals)

write.csv(FRvsdpvals, "fgbdieoff_bm3_wald_VSDandPVALS_FRsym_geneNames.csv", quote=F)

#output for GO and KOG
#AL v. AH. for Orbicella faveolata
rownames(FRresALvAH)<-sub("sym","", rownames(FRresALvAH))
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
rownames(FRresALvU)<-sub("sym","", rownames(FRresALvU))
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
rownames(FRresAHvU)<-sub("sym","", rownames(FRresAHvU))
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


####################################### Plot LME PCoA separated by species
#Start with Ofav
OF=load("resultsLME_OF_syms.Rdata")

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
#lesion    2    32.71 0.06432 1.0785  0.089 .
#genet    16   230.93 0.45406 0.9517  0.662  
#Residual 15   227.47 0.44726                
#Total    34   508.58 1.00000                
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1        

# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors
head(scores[,1])

# plotting PCoA----
margin <- 0.5

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA plot
quartz()
pdf(file="lme_pca_syms_OF_genet.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Breviolum minutum Gene Expression in Orbicella faveolata",
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
       col="red", pch=2, cex=1) +
# plot affected-healthy fav 
points(scores[OFconds$lesion=="AH",xaxis],
		scores[OFconds $lesion=="AH",yaxis],
         col="orange", pch=2, cex=1) +
# plot affected-healthy fav 
points(scores[OFconds $lesion=="U",xaxis],
		scores[OFconds $lesion=="U",yaxis],
         col="blue", pch=2, cex=1)
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

##### now for Ofrank
FR=load("resultsLME_FR_syms.Rdata")

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
#bank      0     0.00 0.00000    NaN       
#lesion    2    43.06 0.06397 1.0546  0.168
#genet    18   280.55 0.41682 0.7635  0.999
#Residual 16   326.63 0.48528              
#Total    37   673.07 1.00000              
              
                 
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
pdf(file="lme_pca_sym_FR_genet.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Breviolum minutum Gene Expression in Orbicella franksi",
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
       col="red", pch=1, cex=1) +
# plot affected-healthy fav 
points(scores[FRconds $lesion=="AH",xaxis],
		scores[FRconds $lesion=="AH",yaxis],
         col="orange", pch=1, cex=1) +
# plot affected-healthy fav 
points(scores[FRconds $lesion=="U",xaxis],
		scores[FRconds $lesion=="U",yaxis],
         col="blue", pch=1, cex=1)
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

################################################## Hypoxia Experiment
# Load data ---------------------------------------------------------------
#Get list of all files in directory
files <- list.files(path="/sym_htseqcounts/")

counts <- read.delim("hypoxia_counts_symD.tsv")

counts <- counts %>% filter(!grepl("_",gene_id)) %>%
  column_to_rownames("gene_id")
dim(counts) #33775    12

conds <- data.frame(
  "sample" = substr(files,1,3),
  "treat" = substr(files,1,1),
  "genet" = substr(files,2,2),
  "rep" = substr(files,3,3)
)

# make sure sample names match up
table(conds$sample == names(counts))
# yep


# SYM ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts = colSums(counts)
min(totalCounts) # 267438
mean(totalCounts) # 1016345
max(totalCounts)  # 1457192

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conds,
  design = ~ treat + genet)

# Set base mean minimum ---------------------------------------------------
means <- apply(counts,1,mean)
table(means>3)
# FALSE  TRUE 
# 17553 16222   

means3 <- names(means[means>3])
head(means3)
length(means3)
#12714

sym.countFilt <- counts[row.names(counts) %in% means3,]
head(sym.countFilt)

totalCountsFilt <- colSums(sym.countFilt)
totalCountsFilt

min(totalCountsFilt) # 265799
max(totalCountsFilt) # 1445031
mean(totalCountsFilt) # 1008498

# check sample order
table(names(sym.countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = sym.countFilt,
  colData = conds,
  design = ~ treat + genet)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("treat", "genet"), force=T)

# outliers that fail 2 tests H5B
out2 <- c("H5B")

# remove out2
dim(sym.countFilt)
sym.count.out <- sym.countFilt %>% select(-one_of(c(out2)))
dim(sym.count.out)
# 16222    11

dim(conds)
head(conds)
conds.out <- conds %>% filter(!sample %in% c(out2))
dim(conds.out)

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = sym.count.out,
  colData = conds.out,
  design = ~ treat + genet)

# remake vsd
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFilt)

#---Results
# log2 fold change (MLE): cond H vs C
res <- results(deds, independentFiltering = F, contrast=c("treat","H","C"))
res

# how many genes pass multiplicity-corrected 0.05 FDR cutoff?
table(res$padj < 0.05)
# FALSE 
#16222 

# Extract log-fold changes --------
table(abs(res$log2FoldChange)>1.5)
# FALSE  TRUE 
# 15909   313 
#write.csv( as.data.frame(res), file="results.csv" ) 

# new vsd ------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)

# Save/Load Data ----------------------------------------------------------
save(conds, sym.count.out, ddsFilt, deds, 
     res,vsd, file="results_sym_hypoxia.Rdata")
load("results_sym_hypoxia.Rdata")

###--------------Get pvals
head(res)
vals <- cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval", "padj")
length(vals[,1])
table(complete.cases(vals))

#Make vsd and pvals table
vsdpvals <- as.data.frame(cbind(assay(vsd),vals))
head(vsdpvals)
dim(vsdpvals)
# 16222    14
table(complete.cases(vsdpvals))
  
write.csv(vsdpvals, "hypoxia_bm3_wald_VSDandPVALS_sym.csv", quote=F)

# Write results for GO/KOG analysis -------------------------------------------
load("results_sym_hypoxia.Rdata")

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
#  763 841  
logs$logP <- logs$logP*sign
logs$gene <- gsub("sym", "", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_hypoxia_host_logP.csv",sep=",")

# PCoA ---------
load("results_sym_hypoxia.Rdata")

# variance stabilized expression data
exp <- data.frame(assay(vsd))
head(exp)

# condition data
table(conds$sam == names(exp))
head(conds)

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA 
set.seed(1)
adonisRes <- adonis2(t(exp)~treat+genet,data=conds,method="manhattan")
adonisRes
#adonis2(formula = t(exp) ~ treat + genet, data = conds, method = "manhattan")
#         Df SumOfSqs      R2      F Pr(>F)    
#treat     1  8850529 0.08859 1.0199  0.349    
#genet     2 21631191 0.21651 1.2463  0.001 ***
#Residual  8 69424901 0.69490                  
#Total    11 99906621 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
           

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
pdf(file="hypoxia_pca_sym.pdf", height = 7, width = 7)
plot(scores[,xaxis], scores[,yaxis],type="n", 
     main = "Durisdinium trenchii Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
ordiellipse(scores,conds$treat,label=F, draw="polygon",
             col = c("orange", "purple4"), alpha = 0.5, border = F)
ordispider(scores,conds$treat,label=F)
points(scores[grep("C1",conds$sam),xaxis],
       scores[grep("C1",conds$sam),yaxis], col="orange", pch=0) +
  points(scores[grep("C2",conds$sam),xaxis],
         scores[grep("C2",conds$sam),yaxis], col="orange", pch= 0) +
  points(scores[grep("C5",conds$sam),xaxis],
         scores[grep("C5",conds$sam),yaxis], col="orange", pch=1) +
  points(scores[grep("H1",conds$sam),xaxis],
         scores[grep("H1",conds$sam),yaxis], col="purple4", pch=1) +
  points(scores[grep("H2",conds$sam),xaxis],
         scores[grep("H2",conds$sam),yaxis], col="purple4", pch=2) +
  points(scores[grep("H5",conds$sam),xaxis],
         scores[grep("H5",conds$sam),yaxis], col="purple4", pch=2)

# legend of sites 
legend("bottomleft", 
       c("Control","Hypoxia"),
       pch=c(2,2), 
       col=c("orange","purple4"), cex=1.5, bty = "n")

legend("bottomright", 
       c("Genet 1","Genet 2", "Genet 5"),
       pch=c(0,1,2), 
       col=c("black","black","black"), cex=1.5, bty = "n")

#insert p value 
legend("topleft", inset=.02, 
       paste("Condition p=",adonisRes$`Pr(>F)`[1], sep=" "), cex=1.5, bty='n') 
legend("topright", inset=.02, 
       paste("Genet p=",adonisRes$`Pr(>F)`[2], sep=" "), cex=1.5, bty='n') 
dev.off()

