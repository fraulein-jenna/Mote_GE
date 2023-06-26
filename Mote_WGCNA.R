#### WGCNA GE Analysis of Mote OA x Temperature Stress Experiment
#November 2021
#Based on Wyatt's April 21 WGCNA script

#### data input and cleaning ####
setwd("~/Desktop/USC/Research Projects/Mote GE Analysis")
options(stringsAsFactors=FALSE)
library(tidyverse)
library(readxl)

#Reading in the table of counts per isogroup by sample - you must open R in the same folder as the input file or change the working directory
# host reads
countsHost=read.table("Cleaned New Ref/AllCountsHost.txt",header=TRUE,row.names=1)%>%
  dplyr::rename("G44_R2_T10"="G44_R2_T210")%>%
  select(everything(), -"G31_R3_T10")%>%
  select(sort(names(.)))

# figure out which samples didn't map to any isogroups
#countsHost_totals<-countsHost%>%
 # summarise_all(sum)
#iso.zero<-apply(countsHost_totals,2,function(x) all(x==0))
#iso.zero 
  
head(countsHost)
length(countsHost[,1]) # 38271 genes old ref, #34842 new ref
names(countsHost)

#create CSV file of GE sample names to compare to key sample names to see where mismatches are
#GEnames <- names(countsHost)
#write.csv(GEnames, file = "GEnames.csv", row.names = F)

#symb reads
countsSymbs=read.table("New Ref Mapping/AllCountsSymbs.txt",header=TRUE,row.names=1)
head(countsSymbs)
length(countsSymbs[,1]) # 39971 genes
names(countsSymbs)<-names(countsHost) # rename the symbiont count columns to match sample names from host counts
names(countsSymbs)
countsSymbs<-countsSymbs%>%select(everything(), -"G31_R3_T10")%>%select(sort(names(.)))

#creating a table with info from experiment - treatments etc.
key = read_xlsx("Datasheets/USC_NSF_2016_Acerv_Processing.xlsx", sheet = "Sheet3")%>%
  select(-"Sample ID")%>%
  unite(Sample, Genotype, "Row/Tank", sep ="_", remove=F)%>%
  unite(Treatment, Temperature, pH, sep = " ", remove =F)

#rename treatments to heat, pH, control, combined
key$Treatment <- sub("Heat Am", "Heat", key$Treatment)
key$Treatment <- sub("Heat Low", "Combined", key$Treatment)
key$Treatment <- sub("Control Am", "Control", key$Treatment)
key$Treatment <- sub("Control Low", "pH", key$Treatment)

key <- key %>% unite(meta, Genotype, Treatment, "Row/Tank", sep = "-", remove = F)%>%
  arrange(Sample)%>%
  filter(Sample %in% names(countsHost))

#create CSV of sample names in the key to compare to the GE sample names to see where mismatches are
sampleNames <- key$Sample
write.csv(sampleNames, file = "sampleNames.csv", row.names = F)

#create a pairing of sample names from counts file with metadata from experimental design
idx<-which(key$Sample %in% colnames(countsHost))
key$Sample[idx]==names(countsHost)
sampleID<-key$meta[idx]

#symbs
idx<-which(key$Sample %in% colnames(countsSymbs))
key$Sample[idx]==names(countsSymbs)
sampleID<-key$meta[idx]


#replacing GE sample names with full exp. metadata info
colnames(countsHost)<-sampleID

#symbs
colnames(countsSymbs)<-sampleID

#before creating the matrix, remove the samples with b, 21, etc.
#we only need unique names for sample names, for the matrix we just need genotype info
key$Genotype <- sub("G3.b", "G3", key$Genotype)
key$Genotype <- sub("G44.20", "G44", key$Genotype)
key$Genotype <- sub("G62.21", "G62", key$Genotype)
key$Genotype <- sub("G44.b", "G44", key$Genotype)
key$Genotype <- sub("G34.21", "G34", key$Genotype)

genotype<-key$Genotype
treatment<-key$Treatment
tank<-key$"Row/Tank"

conditions=data.frame(cbind(genotype,treatment))

####creating DESeq2 dataset - to log transform counts data for WGCNA based pipeline ####
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2) 

#Remove isogroups with low counts from dataset - count less than 10 in more than 90% of samples
countsHost$low = apply(countsHost[,1:190],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each isogroup (host) 
colnames(countsHost)

counts<-countsHost[-which(countsHost$low>171),] #171 is 90% of 190 samples - get rid of count less than 10 in more than 90% of samples

nrow(counts) #16958 genes pass filter => high expression isogroups ## new ref: 16125 genes pass filter

#symbs
countsSymbs$low = apply(countsSymbs[,1:190],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each isogroup (host) 
colnames(countsSymbs)

counts<-countsSymbs[-which(countsSymbs$low>171),] #171 is 90% of 190 samples - get rid of count less than 10 in more than 90% of samples

nrow(counts) #16329 genes pass filter => high expression isogroups

#Create model matrix of genotype*treatment, determine which combinations are zero values 
ml <- model.matrix(~genotype*treatment,conditions)
colnames(ml)
all.zero<-apply(ml,2,function(x) all(x==0))
all.zero # none of them were zero

#remove the zero values from the matrix if you have them 
idx2<-which(all.zero)
ml <- ml[,-idx2]

#now transform counts data using DESeq2, makes structure of data and conditions in a readable format for Deseq2
ddsCOUNTS<-DESeqDataSetFromMatrix(countData=counts[,1:190],colData=conditions,design = ml)
#model matrix isn't full rank because of samples like 44b, 63-21, etc. not being fully replicated across treatments
#so before creating matrix, renamed these to be just 44 etc. because we only need unique names for the samples.

# #rlog transform
rlogCOUNTS<-rlog(ddsCOUNTS,blind=TRUE) #use blind=TRUE to not account for experimental design
head(assay(rlogCOUNTS))
##Table of rlog transformed values for each gene for each sample

#make rlogCOUNTS a dataframe to be easily worked with 
dat=as.data.frame(assay(rlogCOUNTS))
colnames(dat)<-names(counts[1:190])
dat <- dat%>%select(sort(names(.)))
boxplot(dat)
#how do expression plots look overall? most genes are from 7-10 for most samples. 
#Host: about 14 have very high outliers (>20), for Symbs only 2 samples have high outliers.
# new host ref: about 5 have high outliers

##### WGCNA Installation
#BiocManager::install("WGCNA", force = TRUE)

#### Data input, cleaning and pre-processing
library(WGCNA)
disableWGCNAThreads()
options(stringsAsFactors=F)
head(dat)
dim(dat)

#create CSV file of dat sample names to compare to clinA sample names to see where mismatches are
datnames <- names(dat)
write.csv(datnames, file = "datnames.csv", row.names = F)

#transpose expression data
datExpr0=as.data.frame(t(dat[,1:190]))

#check for genes with too many missing values
gsg=goodSamplesGenes(datExpr0, verbose=3)
gsg$allOK #TRUE

##cluster samples to see if there are any obvious outliers
##making experimental trait df for colonies
clinA<-read.csv("Datasheets/Phenotype_metrics_measured.csv",header=T)%>%
  select(-Temp, -pH)%>%
  dplyr::rename("Genotype"="Genotype..")

#rename treatments to heat, pH
clinA$Treatment <- sub("High Temperature", "Heat", clinA$Treatment)
clinA$Treatment <- sub("High pCO2", "pH", clinA$Treatment)

clinA <- clinA %>% unite(ID, Genotype, Treatment, Tank, sep = "-", remove = F)%>%
  arrange(ID)%>%
  filter(ID %in% names(dat))
       

#create CSV file of clinA sample names to compare to dat sample names to see where mismatches are
clinAnames <- clinA$ID
write.csv(clinAnames, file = "clinAnames.csv", row.names = F)

#removing colonies that arent present in RNA data
idx<-which(clinA$ID %in% names(dat))
clinA$ID[idx]==names(dat) #just double check that the names are the same
samID<-clinA$ID[idx]

alltraits<-clinA[idx,]

#we only need unique names for sample names, for the matrix we just need genotype info
alltraits$Genotype <- sub("G3.b", "G3", alltraits$Genotype)
alltraits$Genotype <- sub("G44.20", "G44", alltraits$Genotype)
alltraits$Genotype <- sub("G62.21", "G62", alltraits$Genotype)
alltraits$Genotype <- sub("G44.b", "G44", alltraits$Genotype)
alltraits$Genotype <- sub("G34.21", "G34", alltraits$Genotype)

#Coding in the colony data
Test<-alltraits
names(Test)
for (i in 1:4){ # a loop to store all of the treatments as individual columns denoted 1/0 Y/N
  goi<-unique(Test$Treatment)[i]
  tmp<-data.frame(ifelse(Test$Treatment==goi,1,0))
  colnames(tmp)<-goi
  Test<-cbind(Test,tmp)
}

for (i in 1:12){ # a loop to store all of the genos as individual columns denoted 1/0 Y/N
  goi<-unique(Test$Genotype)[i]
  tmp<-data.frame(ifelse(Test$Genotype==goi,1,0))
  colnames(tmp)<-goi
  Test<-cbind(Test,tmp)
}

Test$GxT<-paste(Test$Genotype,"x",Test$Treatment) #adding in all of the possible interactions terms
for (i in 1:48){
  goi<-unique(Test$GxT)[i]
  tmp<-data.frame(ifelse(Test$GxT==goi,1,0))
  colnames(tmp)<-goi
  Test<-cbind(Test,tmp)
}

head(Test)

Test2=Test[,c(1,5:32,34:81)]
Test=Test[c(1,5:32)]

head(Test)

allTraits<-Test
allInteract<-Test2

# Form a data frame analogous to expression data that will hold the clinical traits
Samples = rownames(datExpr0);
rownames(allTraits) = allTraits$ID;
traitRows = match(Samples, allTraits$ID);
datTraits = allTraits[traitRows, -1];

rownames(allInteract) = allInteract$ID
datInteract=allInteract[,-1]

table(rownames(datTraits)==rownames(datExpr0)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order

#sample dendrogram and trait heat map showing outliers 
library(flashClust)
A=adjacency(t(datExpr0),type="signed")                 #SELECT SIGNED OR UNSIGNED HERE
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

# Remove outlying samples from expression and trait data 
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]
#file
save(datExpr, datTraits,datInteract, file="MoteSymbSamplesAndTraitsAndInteractionssigned_Feb23.RData")

################Moving on!  Network construction and module detection (tutorial PDF Step 2)
#####################################################################

library(WGCNA)
options(stringsAsFactors = FALSE)
library(flashClust)
disableWGCNAThreads()

#host
lnames = load(file="MoteSamplesAndTraitsAndInteractionssigned_Jun22.RData") 
lnames

#Figure out proper SFT 
# # Choose a set of soft-thresholding powers
powers = c(seq(1,18,by=1)); #may need to adjust these power values to hone in on proper sft value
# # Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType="signed", verbose = 2) 
#want smallest value, to plateau closest to 0.9 (but still under)

# # Plot the results:
quartz()
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of .9
abline(h=0.9,col="red")
# # Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

########################## going with *14* for host, *8* for symbs

##################### ran this on HPC for host, see RNAseq_Mote.txt for details #########
###Calculate the adjacencies with the soft thresholding power
softPower=12
adjacency=adjacency(datExpr, power=softPower,type="signed") #must change method type here too!!

#translate the adjacency into topological overlap matrix (TOM) and calculate the corresponding dissimilarity: 
#(to minimize effects of noise and spurious associations)
TOM=TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

save(adjacency, TOM, dissTOM, file="MoteSamplesAndTraits_symbs_TOM_signed_Jun22.RData")

######################open this file instead if adjacency and TOM run on supercomputer, otherwise skip to geneTree:
lnames = load(file="MoteSamplesAndTraits_TOM_signed_Feb23.RData")
#######################################

#####Use hierarchical clustering to produce a tree of genes (based on what exactly? adjacency and similarity/topology)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(3,6)
quartz()
par(mar=c(1,1,1,1))
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
###each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

###Cutting up the tree into highly co-expressed gene modules (by branches)
minModuleSize=35 #we only want large modules, this is considered relatively high
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods) #lists the modules and how many genes are in each one

### HOST
#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
#3245 2006 1276 1158 1034  780  748  659  574  524  447  432  266  257  252  249  233  233  222  196 
#21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40 
#182  174  169  156  154  152  150  107   91   90   88   88   88   87   85   70   62   56   52   49 
#41 
#47 

## NEW HOST
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
#2 3112 1859 1255 1045  828  808  712  639  516  403  345  329  328  305  291  280  277  263  201  189  168 
#22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41 
#154  153  150  145  142  134  128  102   94   91   75   75   74   72   71   71   66   61   58   54 

### SYMBS
#  1    2    3    4    5    6    7    8    9 
#2130 1306  981  732  403  297  252  165   65 

##Plot module assignments under the gene tree
dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors sft=14, min Mod size=35")

#Merge modules whose expression profiles are very similar or choose not to merge
#calculate eigengenes
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 14)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigengenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
#plot
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#start with a 95% similarity merge to get an initial sense of module trait correlations
MEDissThres= 0 # 
MEDissThres= 0.10
MEDissThres=0.15
MEDissThres=0.2

##the lower the height cut, the more the modules have to be similar (more conservative)
#plot the cut line into tree
abline(h=MEDissThres, col="red")
abline(h=MEDissThres, col="blue")
abline(h=MEDissThres, col="green")

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
#merge module colors and find new eigengenes on the new merged modules
mergedColors= merge$colors
mergedMEs= merge$newMEs
length(unique(mergedColors)) #shows how many modules there are now, 

##plot new module colors on gene tree under previous colors to see how they change
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors= mergedColors
#create numerical lables corresponding to the colors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs #a dataframe of module Eigengenes for each module for each sample

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "MoteNetworkSymb_rlog_signed_unmerged_sft14_Feb23.RData")


#################Relating modules to traits and finding important genes
#######################################################################
#can start here if you restarted R session

library(WGCNA)
library(stringr)
library(flashClust)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Load the expression and trait data saved in the first part
lnames = load(file = "MoteSymbSamplesAndTraitsAndInteractionssigned_Feb23.RData") 
#The variable lnames contains the names of loaded variables

# Load network data saved in the second part. Change this based on the level of merge you want to work with
lnames = load(file = "MoteNetworkSymb_rlog_signed_unmerged_sft14_Feb23.RData")


datTraits$Treatment<-as.factor(sapply(str_split(rownames(datTraits),"_"),"[[",1))

#######################Replot module dendrogram
MEList= moduleEigengenes(datExpr, moduleColors, softPower=14)$eigengenes

MEs<- MEList

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
abline(h=.25,col="blue")
abline(h=.21,col="green")

#now for module trait heatmap
#correlate eigengenes with external clinical traits to look for most significant associations
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datTraits1=datTraits[,c(1:16)]
head(datTraits1)

datTraits2=datTraits[,c(13:16)]
head(datTraits2)

# rename the traits to look nice for the heatmap, and reorder them so they are grouped by symb etc
datTraits2=datTraits2%>%
  relocate("Control", .before = "Combined")


# Recalculate MEs with color labels (MEs=module eigengenes)
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=14)$eigengenes 
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes

# for overall heatmap
moduleTraitCor2 = cor(MEs, datTraits1, use = "p"); #p=pearsons
moduleTraitPvalue = corPvalueStudent(moduleTraitCor2, nSamples);
Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#represent module trait correlations as a heatmap
# module-trait correlations
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor2,2)
textMatrix = ps;

#displays only significant p values
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor2)

# Will display correlations and their p-values
#textMatrix = paste(signif(moduleTraitCor2, 2),  "\n(",
#signif(moduleTraitPvalue, 1), ")", sep = "")

# Display the correlation values within a heatmap plot - all modules and all traits
par(mar = c(4, 5.5, 1, 1));
quartz()
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits1),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab = 0.7,
               zlim = c(-0.7,0.7))

# module size barplot
labelShift=300 # increase to move module size labels to the right
quartz()
par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,4500),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

#selecting the 15 largest modules
MEs2 = MEs%>%select("MEgreen","MEpurple","MEblue","MEmagenta","MEbrown","MEyellow", "MEred",
                    "MEblack","MEpink","MEturquoise", "MEgreenyellow", "MEtan", "MEsalmon",
                    "MEcyan", "MEmidnightblue")
# for pretty fig 
moduleTraitCor = cor(MEs2, datTraits2, use = "p"); #p=pearsons
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs2)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor2=cor(MEs2,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor2, nSamples)

### creating a clustered complex heatmap - code from Chille et al. 2021 BMC Genomics
### https://github.com/echille/Mcapitata_Developmental_Gene_Expression_Timeseries/blob/v1.0.0/2a-WGCNA/Developmental_WGCNA.Rmd

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)
#Make list of heatmap row colors
htmap.colors <- names(MEs2)
htmap.colors <- gsub("ME", "", htmap.colors)

mod.sizes<-as.data.frame(mct)%>%
  filter(Freq>290)%>%
  arrange(moduleColors=c("turquoise","blue","yellow","green","red",
                                   "brown","black","pink","magenta","purple",
                         "greenyellow", "tan", "salmon", "cyan", "midnightblue"))

rownames(moduleTraitCor) = gsub("ME", "", rownames(moduleTraitCor))
#rownames(moduleTraitCor) = paste(rownames(moduleTraitCor), sep= " ", mod.sizes$Freq)

quartz()
ht=Heatmap(moduleTraitCor, name = "Corr.", 
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           #right_annotation = size.annot,
           #width = unit(4, "in"), height = unit(8.5, "in"), 
           cluster_columns = FALSE,
           #cluster_rows = METree, row_split = 6, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
             }
             else {
               grid.text(sprintf("-"), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 10),
           column_names_rot = 45,
           row_names_gp = gpar(fontsize = 10, border = FALSE))
draw(ht)

map.grob =  grid.grabExpr(draw(ht)) 

#dev.off()

#Gene relationship to trait and important modules:
#############Treatment
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Heat); 
#change to your trait name of interest: Combined, Control, Heat, pH, etc
names(weight) = "Heat"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")); 
#finds pearson correlations of GE and module eigengenes
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
#making dataframe of pvalues for gene module membership values in each module

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p")); 
# data frame of correlations between expression and trait of interest using pearson correlations
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)); 
#p-values for correlations between expression and trait of interest
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
#geneTraitSignificance is a dataframe of correlation/covariance of GE data and stage value
#GSPvalue is a list of correlation of gTS pvalue

#modules of interest for Combined
moduleCols=c( "yellow", "purple")

#modules of interest for Heat
moduleCols=c("midnightblue", "green")

#plot scatter plots of gene significance vs module membership for all of these modules of interest
#add correlation and p-value, use this to look at how strong the modules are.
quartz()
par(mfrow=c(1,1))
par(mar = c(2, 2, 2, 2));
par(bg = "grey");
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem", module),
                     ylab = "Gene Sig for Heat",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1.2, col = module)
} 


###Visualization of Gene Networks
############################################

#heatmaps of module expression with bar plot of eigengene
#use this to look at the different samples and make sure the modules make sense
#we want blocks by treatment, no single sample driving the differences.

#start with modules that looked strongest: "green","purple","cyan", "tan"

#sort ME by treatment
MEs$sort=rownames(MEs)
MEs.sorted <- MEs%>%
  separate(sort, into=c("garbage", "treatment", "garbage2"), sep= "-")%>%
  dplyr::select(-garbage, -garbage2)%>%
  arrange(treatment)%>%
  dplyr::select(-treatment)

which.module="green" #pick module of interest

ME=MEs.sorted[, paste("ME",which.module, sep="")]
genes=datExpr[,moduleColors==which.module ]

#sort genes by treatment
genes$sort=rownames(genes)
genes.sorted <- genes%>%
  separate(sort, into=c("garbage", "treatment", "garbage2"), sep= "-")%>%
  dplyr::select(-garbage, -garbage2)%>%
  arrange(treatment)%>%
  dplyr::select(-treatment)

quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes.sorted) ),nrgcols=30,rlabels=F, clabels=rownames(genes.sorted), rcols=which.module,)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

##########To output ME by sample

which.module="yellow" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
meyellow<-ME

which.module="purple" 
ME=MEs[, paste("ME",which.module, sep="")]
mepurple<-ME

which.module="midnightblue"  
ME=MEs[, paste("ME",which.module, sep="")]
memidnightblue<-ME

meout<-data.frame(cbind(rownames(datExpr),meyellow, memidnightblue, mepurple))

write.csv(meout,"MEbySample_modsofinterest_cleanref.csv",quote=F,row.names=F) 

# get MEs from mods of interest from previously created ME dataframe 
# sort by treatment for boxplots
MEs$sort=rownames(MEs)
MEs.sorted <- MEs%>%
  separate(sort, into=c("genotype", "treatment", "garbage"), sep= "-")%>%
  dplyr::select(-garbage)%>%
  arrange(treatment)%>%
  select(MEred, MEmagenta, MEcyan, treatment, genotype)%>%
  mutate_at("treatment", as.factor)%>%
  mutate_at("genotype", as.factor)%>%
  mutate(genotype = recode(genotype, "G3.b" = "G3", "G34.21" = "G34", "G44.20" = "G44", "G44.b" = "G44",
                          "G62.21" = "G62", "G63" = "G62"))

medians <- MEs.sorted%>%
  group_by(treatment)%>%
  dplyr::summarise(red = median(MEred, na.rm=TRUE), magenta = median(MEmagenta, na.rm=TRUE),
            cyan= median(MEcyan, na.rm=TRUE), .groups="keep")

treatment.order <- c("Control","Combined","Heat","pH")

cyan <- ggplot(MEs.sorted, aes(x = factor(treatment, level = treatment.order), y = MEcyan))+
  geom_boxplot(fill= "#00FFFF")+
  theme_classic()+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  geom_hline(yintercept =0.06751502, linetype = "dashed", color = "grey")+
  geom_segment(x= "Heat", xend = "Heat", y = 0.06751502, yend =-0.05848323,color = "#D22B2B",
              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  geom_segment(x= "pH", xend = "pH", y = 0.06751502, yend =0.03228030,color = "#6082B6",
               arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  geom_segment(x= "Combined", xend = "Combined", y =0.06751502, yend = -0.05848323,color = "#D22B2B",
             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
  geom_segment(x= "Combined", xend = "Combined", y = -0.05848323, yend = (-0.05848323-(0.06751502-0.03228030)),color = "#6082B6",
               arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))

cyan

magenta <- ggplot(MEs.sorted, aes(x = factor(treatment, level = treatment.order), y = MEmagenta))+
  geom_boxplot(fill= "#FF00FF")+
  theme_classic()+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  scale_y_continuous(breaks=c(-0.1,0, 0.1))+
  geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
  geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
               arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
               arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
  geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
               arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
  geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
               arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
magenta

red <- ggplot(MEs.sorted, aes(x = factor(treatment, level = treatment.order), y = MEred))+
  geom_boxplot(fill= "#FF0000")+
  theme_classic()+
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(breaks=c(-0.1,0, 0.1))+
  geom_hline(yintercept = -0.03785620, linetype = "dashed", color = "grey")+
  geom_segment(x= "Heat", xend = "Heat", y = -0.03785620, yend = 0.06115896,color = "grey",
             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
  #geom_segment(x= "pH", xend = "pH", y = -0.07594524, yend = -0.04762537,color = "#6082B6",
  #             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
red

library(cowplot)
library(grid)
library(gridExtra)

boxes <- plot_grid(cyan, red, magenta, ncol = 3, align = "v", axis="b",rel_heights = c(1, 1,1,1), labels = "AUTO")
boxes

#create common x and y labels
y.grob <- textGrob("Eigengene Expression", 
                   gp=gpar(fontsize=12), rot=90)

x.grob <- textGrob("Treatment", 
                   gp=gpar(fontsize=12))

#add to plot
allplots <-grid.arrange(arrangeGrob(boxes, left = y.grob, bottom = x.grob))

#combine into figure with the heatmap
fig1 <- plot_grid(map.grob, allplots, rel_widths = c(1.4,1), labels = c("A", ""))
fig1


#############Network Analysis to functional annotation and gene ontology###############
# based on Maria's larval WGCNA script

##read in table with the annotated genes for the isogroups we mapped our reads to
# this is from the refs folder on the HPC - saved to personal computer using scp
annot=read.table("Host GO New Ref/Acer_iso2go.tab",header=FALSE,sep="\t") #not all isogroups are annotated 
probes = names(datExpr)
probes2annot = match(probes,annot$V1) 
sum(is.na(probes2annot)) # 6155 NAs (host), 3653 NAs (symb), 6155 new ref host
datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExpr),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
#datOutput=datOutput[order((datOutput$MM.black),decreasing=T),]
#write.table(datOutput,"AnnotatedNetworkAnalysisResultsSymb_rlog_signed_sft8_merge90_Jun22.csv",row.names=F,sep=",")


#Creating files for GO analysis of modules
##############Execute this entire block of code through creating VSD files 
#Host: green, purple, tan, cyan

## can also to rank-based GO with 0 for absent gene and then kME for the genes present in the module.
#below: creating output files for that.

### need to change color name here to module of interest.
red=datOutput%>% 
  select("ProbeID","moduleColors","MM.red")%>%
  mutate(kME = case_when(moduleColors == "red" ~ MM.red))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")

write.csv(red, file = "unmerged_hostGO_MM_red_numeric.csv",row.names = F, quote=F)  

yellow=datOutput%>% 
  select("ProbeID","moduleColors","MM.yellow")%>%
  mutate(kME = case_when(moduleColors == "yellow" ~ MM.yellow))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")

write.csv(yellow, file = "unmerged_hostGO_MM_yellow_numeric.csv",row.names = F, quote=F)

midnightblue=datOutput%>% 
  select("ProbeID","moduleColors","MM.midnightblue")%>%
  mutate(kME = case_when(moduleColors == "midnightblue" ~ MM.midnightblue))%>%
  mutate(kME= ifelse(is.na(kME), 0, kME))%>%
  select("ProbeID", "kME")

write.csv(midnightblue, file = "unmerged_hostGO_MM_midnightblue_numeric.csv",row.names = F, quote=F)

## above files go into the GO MWU analysis

#looking at most significant GO terms for each module - using GO MWU output
red_terms <- read.csv("Host GO New Ref/MWU_BP_unmerged_hostGO_MM_red_numeric.csv",header=T, sep= " ")
purple_terms <- read.csv("Host GO New Ref/MWU_BP_unmerged_hostGO_MM_purple_numeric.csv",header=T, sep= " ")
midnightblue_terms <- read.csv("Host GO New Ref/MWU_BP_unmerged_hostGO_MM_midnightblue_numeric.csv",header=T, sep= " ")

#creating table with top 10 most significant GO terms for each module
#install.packages("flextable")
library(flextable)

red_terms_paper<-red_terms%>%
  arrange(pval)%>%
  select(pval, nseqs, term, name)

red_terms_top10<-red_terms_paper[c(1:10),]
red_terms_top10$module="red"

purple_terms_paper<-purple_terms%>%
  arrange(pval)%>%
  select(pval, nseqs, term, name)

purple_terms_top10<-purple_terms_paper[c(1:10),]
purple_terms_top10$module="purple"

midnightblue_terms_paper<-midnightblue_terms%>%
  arrange(pval)%>%
  select(pval, nseqs, term, name)

midnightblue_terms_top10<-midnightblue_terms_paper[c(1:10),]
midnightblue_terms_top10$module="midnightblue"

top10<-yellow_terms_top10%>%
  full_join(purple_terms_top10)%>%
  full_join(midnightblue_terms_top10)

theme_table <- function(x) {
  border_remove(x) %>% 
    valign(valign = "center", part = "all") %>% 
    align(align = "center", part = "all") %>% 
    fontsize(part = "all", size = 10) %>% 
    bold(part = "header", bold = TRUE) %>%
    bold(part = "body", j = 1, bold = TRUE) %>% 
    color(color = "black", part = "header") %>% 
    bg(part = "header", bg = "transparent")
}


ft <- as_grouped_data(red_terms_top10, groups = c("module"), expand_single = TRUE) %>% 
  as_flextable(hide_grouplabel = FALSE, col_keys = c("pval", "nseqs", "name")) %>% 
  set_header_labels(pval = "P-value", nseqs = "No. Sequences", name = "GO Term") %>%
  set_formatter(pval = function(x){
    formatC(x, format = "e", digits = 1)
  }) %>%
  theme_table() %>% 
  align(i = ~!is.na(name), align = "left", part = "body") %>% 
  bg(i = ~ module %in% "yellow", bg = "#FFFF8F") %>%
  bg(i = ~ module %in% "purple", bg = "#C3B1E1") %>%
  bg(i = ~ module %in% "midnightblue", bg = "#6396D9") %>%
  hline(i = rep(c(TRUE, TRUE, TRUE, TRUE), length = nrow_part(.)))
ft

width(ft, j = 3, width = 10, unit = "in")

save_as_image(width(ft, j = 3, width = 4, unit = "in"), "Final Figures/GO_table_red.pdf",  webshot = "webshot")

## smaller tables for figure 1
yellow_terms_top5 <- yellow_terms_top10[1:5,]
  
yt <- as_grouped_data(yellow_terms_top5, groups = c("module"), expand_single = TRUE) %>% 
  as_flextable(hide_grouplabel = TRUE, col_keys = c("pval", "name")) %>% 
  set_header_labels(pval = "P-value", name = "GO Term") %>%
  set_formatter(pval = function(x){
    formatC(x, format = "e", digits = 3)
  }) %>%
  theme_table() %>% 
  align(i = ~!is.na(name), align = "left", part = "body") %>% 
  bg(i = ~ module %in% "yellow", bg = "#AFE1AF") %>%
  hline(i = rep(c(TRUE, TRUE, TRUE, TRUE), length = nrow_part(.)))
yt

yellow_GO <-gen_grob(yt)

purple_terms_top5 <- purple_terms_top10[1:5,]

pt <- as_grouped_data(purple_terms_top5, groups = c("module"), expand_single = TRUE) %>% 
  as_flextable(hide_grouplabel = TRUE, col_keys = c("pval", "name")) %>% 
  set_header_labels(pval = "P-value", name = "GO Term") %>%
  set_formatter(pval = function(x){
    formatC(x, format = "e", digits = 3)
  }) %>%
  theme_table() %>% 
  align(i = ~!is.na(name), align = "left", part = "body") %>% 
  bg(i = ~ module %in% "purple", bg = "#E0B0FF") %>%
  hline(i = rep(c(TRUE, TRUE, TRUE, TRUE), length = nrow_part(.)))
pt

purple_GO <-gen_grob(pt)

midnightblue_terms_top5 <- midnightblue_terms_top10[1:5,]

mt <- as_grouped_data(midnightblue_terms_top5, groups = c("module"), expand_single = TRUE) %>% 
  as_flextable(hide_grouplabel = TRUE, col_keys = c("pval", "name")) %>% 
  set_header_labels(pval = "P-value", name = "GO Term") %>%
  set_formatter(pval = function(x){
    formatC(x, format = "e", digits = 3)
  }) %>%
  theme_table() %>% 
  align(i = ~!is.na(name), align = "left", part = "body") %>% 
  bg(i = ~ module %in% "midnightblue", bg = "#E0FFFF") %>%
  hline(i = rep(c(TRUE, TRUE, TRUE, TRUE), length = nrow_part(.)))
mt

midnightblue_GO <-gen_grob(mt)

tables <- plot_grid(green_GO, purple_GO, tan_GO, cyan_GO, ncol = 2, rel_heights = c(1, 1,1,1), scale = 0.9)
tables

