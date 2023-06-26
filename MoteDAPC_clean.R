### MOTE GE EXPERIMENT ####
### DESeq2 for DAPC Analysis - Following Quigley & Strader (2022)
### https://github.com/LaserKate/AGF18_RNAseq/blob/main/DESeq2_Juveniles_GE_AGF_NC.R

# October 2022

#setting working directory and loading required packages
setwd("~/Desktop/USC/Research Projects/Mote GE Analysis")

library(DESeq2)
library(readxl)
library(vegan)
library(adegenet)
library(ape)
library(dplyr)
library(stringr)
library(tidyverse)
library(readr)
library(Rmisc)
library(ggplot2)
library(ggridges)
library(cowplot)
library(RColorBrewer)
library(pals)


## import counts and build dataframe including treatment and genotype of each sample
## removing the problem samples and low expression isogroups identified in previous DESEq2 analysis which followed Maria's script

#read in counts data
countsHost=read.table("Cleaned New Ref/AllCountsHost.txt",header=TRUE,row.names=1)%>%
  dplyr::rename("G44_R2_T10"="G44_R2_T210")%>%
  select(everything(), -"G31_R3_T10")%>%
  select(sort(names(.)))

head(countsHost)
length(countsHost[,1]) # 34842isogroups
names(countsHost)

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

#check the sample names in counts and metadata match up
idx<-which(key$Sample %in% colnames(countsHost))
key$Sample[idx]==names(countsHost) ## should all return TRUE or something is out of order

#before creating the matrix, remove the samples with b, 21, etc.
#we only need unique names for sample names, for the matrix we just need genotype info
key$Genotype <- sub("G3.b", "G3", key$Genotype)
key$Genotype <- sub("G44.20", "G44", key$Genotype)
key$Genotype <- sub("G62.21", "G62", key$Genotype)
key$Genotype <- sub("G44.b", "G44", key$Genotype)
key$Genotype <- sub("G34.21", "G34", key$Genotype)
key$Genotype <- sub("G63", "G62", key$Genotype)
key$Genotype <- sub("^G1$", "G01", key$Genotype)
key$Genotype <- sub("^G3$", "G03", key$Genotype)
key$Genotype <- sub("^G7$", "G07", key$Genotype)

genotype<-key$Genotype
treatment<-key$Treatment

conditions=data.frame(cbind(genotype,treatment))%>%
  mutate(as.factor(genotype))%>%
  mutate(as.factor(treatment))

####all the data you ever wanted about quality control
library(arrayQualityMetrics)

cds <- DESeqDataSetFromMatrix(countData = countsHost,
                              colData = conditions,
                              design = ~ genotype + treatment)

cds = estimateSizeFactors(cds)
cds=estimateDispersions(cds)

par(mar=c(1, 1, 1, 1))
quartz()
plotDispEsts(cds)
vsdBlind=varianceStabilizingTransformation(cds)

v="~/Desktop/USC/Research Projects/Mote GE Analysis/AQM_clean"
e=ExpressionSet(assay(vsdBlind), AnnotatedDataFrame(as.data.frame(colData(vsdBlind))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treatment","genotype"),force=TRUE) #check .html output file in new folder

## sample 133 ticked 2 outlier boxes, removing and re-running.

###############Remove outliers as detected; repeat arrayQualityMetrics after regenerating newCountDataSet
counts.nobad1=countsHost[,-c(133)]
conditions.nobad1=conditions[-c(133),]

#DESeq will calculate the log2fold and Wald test p-value for last variable in design formula
cds2 <- DESeqDataSetFromMatrix(countData = counts.nobad1,
                               colData = conditions.nobad1,
                               design = ~ genotype + treatment)

cds2=estimateSizeFactors(cds2)
cds2=estimateDispersions(cds2)
vsdBlind2=varianceStabilizingTransformation(cds2)

v="~/Desktop/USC/Research Projects/Mote GE Analysis/AQM_clean2"
e2=ExpressionSet(assay(vsdBlind2), AnnotatedDataFrame(as.data.frame(colData(vsdBlind2))))
arrayQualityMetrics(e2,outdir=v,intgroup=c("treatment", "genotype"),force=TRUE)
# now this looks good!

# remove low expression isogroups
counts.nobad1$low = apply(counts.nobad1[,1:189],1,function(x){sum(x<=10)})  #making new column counting number of samples with counts <=10 within each isogroup (host)
counts.nobad1<-counts.nobad1[-which(counts.nobad1$low>170),1:189] #170 is 90% of 189 samples
nrow(counts.nobad1) #16125 high expression isogroups
head(counts.nobad1)

counts = counts.nobad1
colData = conditions.nobad1

#this looks good so moving ahead with the analysis!

### add in grouping by different factors to colData
colData$group1=factor(paste0(colData$genotype,colData$treatment))

#all treatments, no interactions ###

dds<-DESeqDataSetFromMatrix(counts,
                            colData = colData, 
                            design = ~genotype+treatment)
vst=vst(dds)
vstAll.df = assay(vst)

######## test for combined vs control treatment effect 
dds<-DESeq(dds, minReplicatesForReplace=Inf) 
resCombined<-results(dds, contrast=c('treatment', 'Combined', 'Control')) #here is where the two contrasting conditions get defined
mcols(resCombined,use.names=TRUE)
summary(resCombined)
#out of 16125 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 4540, 28%
#LFC < 0 (down)     : 5232, 32%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resCombined$padj < 0.05)
#FALSE  TRUE 
# 7344  8781  
res=data.frame(cbind("gene"=row.names(resCombined),"stat"= resCombined$stat))
head(res)
write.csv(res,file="resCombined_stat_2023.csv",quote=F, row.names=F)

save(vstAll.df, resCombined, colData,file="MoteGE_WaldTreatmentGeno_clean.Rdata")

###FULL MODEL - INTERACTIONS ###

dds<-DESeqDataSetFromMatrix(counts,
                            colData = colData, 
                            design = ~group1)

vstInteract=vst(dds)
vstInteract.df = assay(vstInteract)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 

#how much of the vsriance is explained by treatment and geno in the full model?
adFull = adonis(t(vstInteract.df)~ treatment+genotype, data = colData, method="manhattan")
adFull$aov.tab

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#treatment   3 1004695075 334898358  8.1864 0.10204  0.001 ***
#  genotype   10 1682281490 168228149  4.1122 0.17086  0.001 ***
#  Residuals 175 7159105679  40909175         0.72710           
#Total     188 9846082244                   1.00000     

# Moving on to DAPC

#### COMBINED TREATMENT #####
# need to create a counts table and colData that only includes the control and combined treatment samples

colData<-colData%>%mutate_at("treatment", as.factor)
colDataComp=colData[colData$treatment=="Control"|colData$treatment=="Combined",]
colDataComp$group2=factor(paste0(colDataComp$genotype,colDataComp$treatment))

#create a new key with only the samples of interest
keyComp<-key%>%filter(Treatment =="Control"|Treatment =="Combined")

#use this key to filter the counts table 
countsComp<-counts[, which((names(counts) %in% keyComp$Sample)==TRUE)]

dds<-DESeqDataSetFromMatrix(countsComp,
                            colData = colDataComp, 
                            design = ~group2)

vstComp=vst(dds)
vstComp.df = assay(vstComp)

#how much variance explained by treatment vs geno now?
adFull = adonis(t(vstComp.df)~ treatment+genotype, data = colDataComp, method="manhattan")
adFull$aov.tab

#Permutation: free

#Number of permutations: 999

#Terms added sequentially (first to last)

#Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#treatment   1  588160752 588160752 16.6064 0.12375  0.001 ***
#  genotype   10 1012442993 101244299  2.8586 0.21302  0.001 ***
#  Residuals  89 3152183480  35417792         0.66323           
#Total     100 4752787225                   1.00000           


dim(vstComp.df) #16125   101

# building DAPC 
dapc1 <- dapc(t(vstComp.df), colDataComp$group2)
## choose to retain 55 PCs, 5 discriminant functions
dapc1
## this explains 92.9% of the conserved variance

## for the vst, they suggest retaining 17 PCs, which explains 72% of the conserved variance
temp <- optim.a.score(dapc1, n.sim = 1) 
combined.whole.dapc <- dapc(t(vstComp.df), colDataComp$group2, n.da=4, n.pca=17)
combined.whole.dapc

#quartz()
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=FALSE,solid=.4)
#scatter(dapc,1,1,bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

combined.whole.varexpl <- round((combined.whole.dapc$eig/sum(combined.whole.dapc$eig))[1:2] * 100, 1) #66 19.4

combined.whole.dapc1 <- tibble(sample = rownames(combined.whole.dapc$ind.coord),
                grp = combined.whole.dapc$grp,
                LD1 = combined.whole.dapc$ind.coord[,1],
                LD2 = combined.whole.dapc$ind.coord[,2])
combined.whole.dapc2 <- combined.whole.dapc1 %>%
  group_by(grp)%>%
  dplyr::summarize(c1 = mean(LD1),
                   c2 = mean(LD2), .groups="keep")%>%
  full_join(combined.whole.dapc1)

### figures - combined whole genome

combined.whole.dapc.fig <-   ggplot(combined.whole.dapc2, aes(fill = factor(str_sub(grp, 1,3)), 
                                shape = factor(str_sub(grp, 4)))) +
  geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
  geom_point(aes(x = c1, y = c2), size = 3) +
  geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
  scale_shape_manual(name="Treatment", labels = c("Combined","Control"), values = c(24,21)) +
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_continuous(limits = c(-15,25)) +
  scale_y_continuous(limits = c(-10,20)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("LD1 [", combined.whole.varexpl[1],"%]"), y = paste0("LD2 [", combined.whole.varexpl[2],"%]")) +
  theme_bw()+
  theme(legend.position = "none")
combined.whole.dapc.fig

combined.whole.ridges=ggplot(combined.whole.dapc1, aes(x = LD2, fill = factor(str_sub(grp, 4)), y = factor(str_sub(grp, 1,3)), height = ..density..)) +
  geom_density_ridges(scale = 1, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(limits = c(-15,20)) +
  scale_fill_manual(name="Treatment", values =c("#E6E6FA","#D3D3D3")) +
  ylab("Genotype")+
  theme_ridges() 
combined.whole.ridges

# get the genes that are most important to the LD2 from DAPC
combined.whole.geneloadings <-as.data.frame(combined.whole.dapc$pca.loadings)%>%
  select("PCA-pa.2")

#want to know if any of these are annotated
iso2gene=read.table('~/Desktop/USC/Research Projects/Mote GE Analysis/New Ref Mapping/Acer_iso2gene.tab',sep = '\t', colClasses = "character")

#select top 20 on either end
combined.whole.pos_genes<-combined.whole.geneloadings%>%filter(.>0.02412101)
combined.whole.neg_genes<-combined.whole.geneloadings%>%filter(.< -0.03603515)

combined.whole.PC2_pos_genes=iso2gene[iso2gene$V1 %in% rownames(combined.whole.pos_genes),] 
combined.whole.PC2_pos_genes$V2
# 6 annotated
#[1] "PiggyBac transposable element-derived protein 4 OS=Homo sapiens OX=9606 GN=PGBD4 PE=2 SV=3 E(blastx)=3e-44"
#[2] "G-protein-signaling modulator 2 OS=Homo sapiens OX=9606 GN=GPSM2 PE=1 SV=3 E(blastx)=4e-30"    
#[3] "G-protein-signaling modulator 2 OS=Mus musculus OX=10090 GN=Gpsm2 PE=1 SV=2 E(blastx)=5e-35"  
#[4] "Betaine--homocysteine S-methyltransferase 1 OS=Danio rerio OX=7955 GN=bhmt PE=2 SV=1 E(blastx)=1e-154"     
#[5] "G-protein-signaling modulator 2 OS=Homo sapiens OX=9606 GN=GPSM2 PE=1 SV=3 E(blastx)=8e-34"                
#[6] "Tetratricopeptide repeat protein 28 OS=Mus musculus OX=10090 GN=Ttc28 PE=1 SV=3 E(blastx)=3e-21" 

combined.whole.PC2_neg_genes=iso2gene[iso2gene$V1 %in% rownames(combined.whole.neg_genes),] 
combined.whole.PC2_neg_genes$V2
# 3 annotated
#[1] "40S ribosomal protein S15a OS=Rattus norvegicus OX=10116 GN=Rps15a PE=1 SV=2 E(blastx)=2e-76"
#[2] "RING-box protein 1 OS=Salmo salar OX=8030 GN=rbx1 PE=2 SV=2 E(blastx)=2e-57"                 
#[3] "Transmembrane protein 138 OS=Bos taurus OX=9913 GN=TMEM138 PE=2 SV=1 E(blastx)=5e-48"  

### creating dataframe with the physiological data from combined vs. control treatments
##making experimental trait df 
traits<-read.csv("Datasheets/Phenotype_metrics_measured.csv",header=T)%>%
  dplyr::rename("Genotype"="Genotype..")%>%
  mutate(sample = paste0(Genotype,sep="_",Tank))%>%
  select(-Tank)%>%
  filter(Treatment=="Control"|Treatment=="Combined"|Treatment=="High Temperature")

distribution_ridges=ggplot(traits, aes(x = Symbionts.cm2, fill = Treatment, y = Genotype, height = ..density..)) +
  geom_density_ridges(scale = 1.5, stat = "density") +
  #scale_y_discrete(expand = c(0.01, 0)) +
  #scale_x_reverse()+
  scale_fill_manual(values =c("#E6E6FA","#D3D3D3", "#F88379")) +
  ylab("Genotype")+
  theme_ridges()
  #theme(legend.position="none")
distribution_ridges

#rename genotypes
traits$Genotype <- sub("G63", "G62", traits$Genotype)
traits$Genotype <- sub("G3.b", "G3", traits$Genotype)
traits$Genotype <- sub("G34.21", "G34", traits$Genotype)
traits$Genotype <- sub("G44.b", "G44", traits$Genotype)
traits$Genotype <- sub("G44.20", "G44", traits$Genotype)
traits$Genotype <- sub("G62.21", "G62", traits$Genotype)
traits$Genotype <- sub("^G1$", "G01", traits$Genotype)
traits$Genotype <- sub("^G3$", "G03", traits$Genotype)
traits$Genotype <- sub("^G7$", "G07", traits$Genotype)

# summarize trait data by genotype in the combined treatment
combined_symbiont_traits<-traits%>%select(Genotype,Treatment,Symbionts.cm2)%>%
  group_by(Genotype,Treatment)%>%
  dplyr::summarize(mean=mean(Symbionts.cm2, na.rm=TRUE), .groups="keep")%>%
  spread(Treatment, mean)%>%
  mutate(symb.norm=Combined/Control)%>%
  mutate(symb.diff=Combined-Control)%>%
  select(Genotype, symb.norm, symb.diff, Control)

# keeping trait data unique to fragment level, but normalizing by the mean of the control.
combined_individual_traits<-traits%>%select(Genotype,Treatment,Symbionts.cm2,sample)%>%
  filter(Treatment == "Combined")%>%
  spread(Treatment, Symbionts.cm2)%>%
  group_by(Genotype)%>%
  full_join(combined_symbiont_traits, by = "Genotype")%>%
  select(sample, Genotype, Combined, Control)%>%
  dplyr::rename("Control.mean"="Control")%>%
  mutate(Symb.normalized= Combined/Control.mean)%>%
  mutate(Sym.difference = Combined-Control.mean)

# next, looking at extracting plasticity and frontloading values from LD2, combine with symbiont df
# plasticity = distance between distributions for each genotype
# baseline = mean of the control distribution for each genotype
combined.whole.expression<- combined.whole.dapc2%>%select(grp,c2)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  ungroup()%>%
  select(treatment, genotype, c2)%>%
  unique()%>%
  spread(treatment, c2)%>%
  mutate(plasticity = (Combined-Control))%>%
  dplyr::rename("baseline"="Control")%>%
  dplyr::rename("Genotype"="genotype")%>%
  select(everything(), -"Combined")%>%
  full_join(combined_symbiont_traits, by = "Genotype")

#individual sample plasticity dataframe
combined.whole.frag <-combined.whole.dapc2%>%select(grp,sample,LD2)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  filter(treatment == "Combined")%>%
  ungroup()%>%
  select(sample, genotype, LD2)%>%
  dplyr::rename("combined.expression"="LD2")%>%
  dplyr::rename("Genotype"='genotype')%>%
  full_join(combined.whole.expression, by = "Genotype")%>%
  select(sample, Genotype, combined.expression, baseline)%>%
  mutate(frag.plasticity = combined.expression-baseline)%>%
  full_join(combined_individual_traits, by = c("Genotype"="Genotype", "sample"="sample"))%>%
  drop_na()%>%
  filter(Genotype!="G34")

### linear models

## symbiont difference - linear regression
combined.whole.baseline_diff<-lm(symb.diff~baseline, data = combined.whole.expression)
summary(combined.whole.baseline_diff)
# NS: p = 0.311

combined.whole.plasticity_diff<-lm(symb.diff~plasticity, data = combined.whole.expression)
summary(combined.whole.plasticity_diff)
#NS: p = 0.2972

#individual frag data
individual_whole_diff_pl<-lm(Sym.difference~frag.plasticity, data = combined.whole.frag)
summary(individual_whole_diff_pl)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)      -315067     151880  -2.074   0.0438 *
#  frag.plasticity   -78443      34136  -2.298   0.0263 *

#Residual standard error: 379800 on 45 degrees of freedom
#Multiple R-squared:  0.105,	Adjusted R-squared:  0.08514 
#F-statistic: 5.281 on 1 and 45 DF,  p-value: 0.02627

wholecomb.pl.plot <- ggplot(combined.whole.frag, colour = Genotype, aes(x=frag.plasticity, y = Sym.difference))+
  geom_smooth(method = "lm", se = TRUE, colour = "grey60")+
  geom_point(shape = 21, aes(x=frag.plasticity, y = Sym.difference, fill = Genotype))+
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  ylab("Difference in Symbiont Density/cm2")+
  xlab("Plasticity")+
  theme_bw()
wholecomb.pl.plot

plot(individual_whole_diff_pl)

individual_whole_diff_bl<-lm(Sym.difference~baseline, data = combined.whole.frag)
summary(individual_whole_diff_bl)
# NS: 0.0555

## normalized symbs - quasibinomial
combined.whole.baseline_norm<-glm(symb.norm~baseline, data = combined.whole.expression, family = quasibinomial(link = 'logit'))
summary(combined.whole.baseline_norm)
#NS: p = 0.3984

combined.whole.plasticity_norm<-glm(symb.norm~plasticity, data = combined.whole.expression, family = quasibinomial(link = 'logit'))
summary(combined.whole.plasticity_norm)
#NS: p = 0.203

combined.whole.plasticity_diff<-lm(symb.diff~plasticity, data = combined.whole.expression)
summary(combined.whole.plasticity_diff)
#NS: p = 0.2972

#individual frag data
individual_whole_norm_pl<-glm(Sym.difference~frag.plasticity, data = combined.whole.frag, family=gaussian(link='log'))
summary(individual_whole_diff_pl)

### MODULE YELLOW - COMBINED TREATMENT ####
## first, load in the gene list of green module output from WGCNA
modYellow<-read.csv("Host GO New Ref/unmerged_hostGO_MM_yellow_numeric.csv",header=T)%>%
  filter(kME > 0)

# now filter the counts data by isogroups that are present in the module
yellowCounts=counts[rownames(counts) %in% modYellow$ProbeID,]

# filter the yellowCounts to be only combined and control using keyComp
yellowComp<-yellowCounts[, which((names(yellowCounts) %in% keyComp$Sample)==TRUE)]

dds<-DESeqDataSetFromMatrix(yellowComp,
                            colData = colDataComp, 
                            design = ~group2)

vstYellow=varianceStabilizingTransformation(dds)
vstYellow.df = assay(vstYellow)

dim(vstYellow.df) #1045  101

#how much variation in this module is explained by treatment vs geno?
adYellow = adonis(t(vstYellow.df)~ treatment+genotype, data = colDataComp, method="manhattan")
adYellow$aov.tab

#Permutation: free
#Number of permutations: 999

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treatment   1   8876184 8876184 27.6116 0.19994  0.001 ***
 # genotype   10   6906565  690657  2.1485 0.15558  0.001 ***
#  Residuals  89  28610496  321466         0.64448           
#Total     100  44393245                 1.00000       

dapcY <- dapc(t(vstYellow.df), colDataComp$group2)
## choose to retain 60 PCs, 5 discriminant functions
dapcY
## this explains 94.1% of the conserved variance

## for the vst, they suggest retaining 12 PCs, which explains 61.8% of the conserved variance
temp <- optim.a.score(dapcY, n.sim = 1) 
combined.yellow.dapc <- dapc(t(vstYellow.df), colDataComp$group2, n.da=4, n.pca=12)
combined.yellow.dapc
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=FALSE,solid=.4)
#scatter(dapc,1,1,bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

combined.yellow.varexpl <- round((combined.yellow.dapc$eig/sum(combined.yellow.dapc$eig))[1:2] * 100, 1) #66 19.4

combined.yellow.dapc1 <- tibble(sample = rownames(combined.yellow.dapc$ind.coord),
                 grp = combined.yellow.dapc$grp,
                 LD1 = combined.yellow.dapc$ind.coord[,1],
                 LD2 = combined.yellow.dapc$ind.coord[,2])
combined.yellow.dapc2 <- combined.yellow.dapc1 %>%
  group_by(grp)%>%
  dplyr::summarize(c1 = mean(LD1),
                   c2 = mean(LD2), .groups="keep")%>%
  full_join(combined.yellow.dapc1)

## visualization
yellow.dapc<-   ggplot(combined.yellow.dapc2, aes(fill = factor(str_sub(grp, 1,3)), 
                                   shape = factor(str_sub(grp, 4)))) +
  geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
  geom_point(aes(x = c1, y = c2), size = 3) +
  geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
  scale_shape_manual(name="Treatment", labels = c("Combined","Control"), values = c(24,21)) +
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_continuous(limits = c(-15,25)) +
  scale_y_continuous(limits = c(-10,20)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("LD1 [", combined.yellow.varexpl[1],"%]"), y = paste0("LD2 [", combined.yellow.varexpl[2],"%]")) +
  theme_bw()
yellow.dapc

yellow_ridges=ggplot(combined.yellow.dapc2, aes(x = LD1, fill = factor(str_sub(grp, 4)), y = factor(str_sub(grp, 1,3)), height = ..density..)) +
  geom_density_ridges(scale = 1.5, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_reverse()+
  scale_fill_manual(name="Treatment", values =c("#E6E6FA","#D3D3D3")) +
  ylab("Genotype")+
  theme_ridges() +
  theme(legend.position="none")
yellow_ridges

#pulling out sig gene loadings for PC2 (treatment axis) in the yellow module
geneloadings_yellow <-as.data.frame(combined.yellow.dapc$pca.loadings)%>%
  select("PCA-pa.1")

#select top 20 on either end
yellow.pos_genes<-geneloadings_yellow%>%filter(.>	0.08539806)
yellow.neg_genes<-geneloadings_yellow%>%filter(.< -0.02130230)

yellowPC1_pos_genes=iso2gene[iso2gene$V1 %in% rownames(yellow.pos_genes),] 
yellowPC1_pos_genes$V2
# 2 annotated
#[1] "Secreted acidic protein 1A OS=Acropora millepora OX=45264 PE=1 SV=1 E(blastx)=5e-54"
#[2] "Calumenin OS=Mus musculus OX=10090 GN=Calu PE=1 SV=1 E(blastx)=1e-23"  

yellowPC1_neg_genes=iso2gene[iso2gene$V1 %in% rownames(yellow.neg_genes),] 
yellowPC1_neg_genes$V2
# 12 annotated
#[1] "Coiled-coil domain-containing protein 85C OS=Xenopus tropicalis OX=8364 GN=ccdc85c PE=2 SV=1 E(blastx)=2e-36"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
#[2] "Heterochromatin protein 1 OS=Drosophila virilis OX=7244 GN=HP1A PE=3 SV=1 E(blastx)=2e-29"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
#[3] "Endoplasmic reticulum transmembrane helix translocase OS=Mus musculus OX=10090 GN=Atp13a1 PE=1 SV=2 E(blastx)=4e-95"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
#[4] "Polynucleotide 5-hydroxyl-kinase NOL9 OS=Mus musculus OX=10090 GN=Nol9 PE=1 SV=1 E(blastx)=2e-72\nisogroupDN72532_c9_g2\tAcetoacetyl-CoA synthetase OS=Mus musculus OX=10090 GN=Aacs PE=1 SV=1 E(blastx)=1e-178\nisogroupDN72485_c2_g1\tGlutamine-rich protein 2 OS=Mus musculus OX=10090 GN=Qrich2 PE=2 SV=1 E(blastx)=4e-06\nisogroupDN78416_c0_g1\tHypoxia-inducible factor 1-alpha inhibitor OS=Danio rerio OX=7955 GN=hif1an PE=2 SV=2 E(blastx)=1e-10\nisogroupDN13462_c0_g1\tProtoporphyrinogen oxidase OS=Homo sapiens OX=9606 GN=PPOX PE=1 SV=1 E(blastx)=5e-122\nisogroupDN62126_c0_g2\tMembrane progestin receptor delta OS=Homo sapiens OX=9606 GN=PAQR6 PE=1 SV=2 E(blastx)=4e-28\nisogroupDN69528_c5_g1\t26S proteasome non-ATPase regulatory subunit 8 (Fragment) OS=Pongo abelii OX=9601 GN=PSMD8 PE=2 SV=2 E(blastx)=4e-93\nisogroupDN71526_c0_g1\tSmall RNA 2-O-methyltransferase OS=Mus musculus OX=10090 GN=Henmt1 PE=1 SV=1 E(blastx)=2e-63"
#[5] "Protein FRG1 OS=Takifugu rubripes OX=31033 GN=frg1 PE=2 SV=1 E(blastx)=2e-60"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
#[6] "Deleted in malignant brain tumors 1 protein OS=Homo sapiens OX=9606 GN=DMBT1 PE=1 SV=2 E(blastx)=6e-102"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
#[7] "PRKCA-binding protein OS=Rattus norvegicus OX=10116 GN=Pick1 PE=1 SV=1 E(blastx)=2e-175"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
#[8] "Ribonuclease P/MRP protein subunit POP5 OS=Mus musculus OX=10090 GN=Pop5 PE=1 SV=1 E(blastx)=1e-32"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
#[9] "Protein regulator of cytokinesis 1 OS=Homo sapiens OX=9606 GN=PRC1 PE=1 SV=2 E(blastx)=4e-107"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
#[10] "DNA polymerase delta subunit 3 OS=Gallus gallus OX=9031 GN=POLD3 PE=1 SV=2 E(blastx)=1e-36"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
#[11] "Extracellular tyrosine-protein kinase PKDCC OS=Homo sapiens OX=9606 GN=PKDCC PE=1 SV=2 E(blastx)=2e-74"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
#[12] "Histone H3-like centromeric protein A OS=Danio rerio OX=7955 GN=cenpa PE=1 SV=1 E(blastx)=1e-40"   

# next, looking at extracting plasticity and frontloading values from LD2, combine with symbiont df
# plasticity = distance between distributions for each genotype
# baseline = mean of the control distribution for each genotype
combined.yellow.expression<- combined.yellow.dapc2%>%select(grp,c1)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  ungroup()%>%
  select(treatment, genotype, c1)%>%
  unique()%>%
  spread(treatment, c1)%>%
  mutate(plasticity = (Combined-Control))%>%
  dplyr::rename("baseline"="Control")%>%
  dplyr::rename("Genotype"="genotype")%>%
  select(everything(), -"Combined")%>%
  full_join(combined_symbiont_traits, by = "Genotype")

#individual frag info
combined.yellow.frag <-combined.yellow.dapc2%>%select(grp,sample,LD1)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  filter(treatment == "Combined")%>%
  ungroup()%>%
  select(sample, genotype, LD1)%>%
  dplyr::rename("combined.expression"="LD1")%>%
  dplyr::rename("Genotype"='genotype')%>%
  full_join(combined.yellow.expression, by = "Genotype")%>%
  select(sample, Genotype, combined.expression, baseline)%>%
  mutate(frag.plasticity = combined.expression-baseline)%>%
  full_join(combined_individual_traits, by = c("Genotype"="Genotype", "sample"="sample"))%>%
  drop_na()

### linear models

## symbiont difference - linear regression
combined.yellow.baseline_diff<-lm(symb.diff~baseline, data = combined.yellow.expression)
summary(combined.yellow.baseline_diff)
# NS: p = 0.859

combined.yellow.plasticity_diff<-lm(symb.diff~plasticity, data = combined.yellow.expression)
summary(combined.yellow.plasticity_diff)
#NS: p = 0.991

#individual frag data
individual_yellow_diff_pl<-lm(Sym.difference~frag.plasticity, data = combined.yellow.frag)
summary(individual_yellow_diff_pl)
#NS: 0.152

individual_yellow_diff_bl<-lm(Sym.difference~baseline, data = combined.yellow.frag)
summary(individual_yellow_diff_bl)
#NS: 0.314

## normalized symbs - quasibinomial
combined.yellow.baseline_norm<-glm(symb.norm~baseline, data = combined.yellow.expression, family = quasibinomial(link = 'logit'))
summary(combined.yellow.baseline_norm)
#NS: p = 0.1198

combined.yellow.plasticity_norm<-glm(symb.norm~plasticity, data = combined.yellow.expression, family = quasibinomial(link = 'logit'))
summary(combined.yellow.plasticity_norm)
#NS: p = 0.481

### MODULE PURPLE - COMBINED TREATMENT ####
## first, load in the gene list of green module output from WGCNA
modPurple<-read.csv("Host GO New Ref/unmerged_hostGO_MM_purple_numeric.csv",header=T)%>%
  filter(kME > 0)

# now filter the counts data by isogroups that are present in the module
purpleCounts=counts[rownames(counts) %in% modPurple$ProbeID,]

# filter the purpleCounts to be only combined and control using keyComp
purpleComp<-purpleCounts[, which((names(purpleCounts) %in% keyComp$Sample)==TRUE)]

dds<-DESeqDataSetFromMatrix(purpleComp,
                            colData = colDataComp, 
                            design = ~group2)

vstPurple=varianceStabilizingTransformation(dds)
vstPurple.df = assay(vstPurple)

dim(vstPurple.df) #403  101

#how much variation in this module is explained by treatment vs geno?
adPurple = adonis(t(vstPurple.df)~ treatment+genotype, data = colDataComp, method="manhattan")
adPurple$aov.tab

#Permutation: free
#Number of permutations: 999

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treatment   1    741821  741821 29.7954 0.21270  0.001 ***
#  genotype   10    530004   53000  2.1288 0.15196  0.001 ***
#  Residuals  89   2215852   24897         0.63534           
#Total     100   3487677                 1.00000      

dapcP <- dapc(t(vstPurple.df), colDataComp$group2)
## choose to retain 60 PCs, 5 discriminant functions
dapcP
## this explains 96.4% of the conserved variance

## for the vst, they suggest retaining 22 PCs, which explains 61.8% of the conserved variance
temp <- optim.a.score(dapcP, n.sim = 1) 
combined.purple.dapc <- dapc(t(vstPurple.df), colDataComp$group2, n.da=4, n.pca=22)
combined.purple.dapc
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=FALSE,solid=.4)
#scatter(dapc,1,1,bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

combined.purple.varexpl <- round((combined.purple.dapc$eig/sum(combined.purple.dapc$eig))[1:2] * 100, 1) #66 19.4

combined.purple.dapc1 <- tibble(sample = rownames(combined.purple.dapc$ind.coord),
                                grp = combined.purple.dapc$grp,
                                LD1 = combined.purple.dapc$ind.coord[,1],
                                LD2 = combined.purple.dapc$ind.coord[,2])
combined.purple.dapc2 <- combined.purple.dapc1 %>%
  group_by(grp)%>%
  dplyr::summarize(c1 = mean(LD1),
                   c2 = mean(LD2), .groups="keep")%>%
  full_join(combined.purple.dapc1)

## visualization
purple.dapc<-   ggplot(combined.purple.dapc2, aes(fill = factor(str_sub(grp, 1,3)), 
                                                  shape = factor(str_sub(grp, 4)))) +
  geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
  geom_point(aes(x = c1, y = c2), size = 3) +
  geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
  scale_shape_manual(name="Treatment", labels = c("Combined","Control"), values = c(24,21)) +
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_continuous(limits = c(-15,25)) +
  scale_y_continuous(limits = c(-10,20)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("LD1 [", combined.purple.varexpl[1],"%]"), y = paste0("LD2 [", combined.purple.varexpl[2],"%]")) +
  theme_bw()
purple.dapc

purple_ridges=ggplot(combined.purple.dapc2, aes(x = LD1, fill = factor(str_sub(grp, 4)), y = factor(str_sub(grp, 1,3)), height = ..density..)) +
  geom_density_ridges(scale = 2, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_reverse(limits = c(12,-11))+
  scale_fill_manual(name="Treatment", values =c("#E6E6FA","#D3D3D3")) +
  ylab("Genotype")+
  theme_ridges() +
  theme(legend.text=element_text(size=10),
        legend.title=element_text(size=11))
purple_ridges

#pulling out sig gene loadings for PC2 (treatment axis) in the yellow module
geneloadings_purple <-as.data.frame(combined.purple.dapc$pca.loadings)%>%
  select("PCA-pa.1")

#select top 20 on either end
purple.pos_genes<-geneloadings_purple%>%filter(.>	2.403210e-02)
purple.neg_genes<-geneloadings_purple%>%filter(.< -0.10605473)

purplePC1_pos_genes=iso2gene[iso2gene$V1 %in% rownames(purple.pos_genes),] 
purplePC1_pos_genes$V2
# 7 annotated
#[1] "Dual specificity protein phosphatase 7 OS=Homo sapiens OX=9606 GN=DUSP7 PE=1 SV=4 E(blastx)=5e-85"        
#[2] "Protein phosphatase 1 regulatory subunit 27 OS=Mus musculus OX=10090 GN=Ppp1r27 PE=2 SV=1 E(blastx)=8e-20"
#[3] "Ribitol-5-phosphate transferase FKTN OS=Macaca fascicularis OX=9541 GN=FKTN PE=1 SV=1 E(blastx)=2e-104"   
#[4] "E3 ubiquitin ligase RNF157 OS=Homo sapiens OX=9606 GN=RNF157 PE=1 SV=3 E(blastx)=4e-103"                  
#[5] "Sorting nexin-19 OS=Mus musculus OX=10090 GN=Snx19 PE=1 SV=1 E(blastx)=2e-17"                             
#[6] "Glutamate carboxypeptidase 2 OS=Rattus norvegicus OX=10116 GN=Folh1 PE=1 SV=1 E(blastx)=2e-163"           
#[7] "Probable ATP-dependent RNA helicase DDX20 OS=Danio rerio OX=7955 GN=ddx20 PE=3 SV=1 E(blastx)=2e-169" 

purplePC1_neg_genes=iso2gene[iso2gene$V1 %in% rownames(purple.neg_genes),] 
purplePC1_neg_genes$V2
# 3 annotated
#[1] "Aldehyde dehydrogenase family 3 member B1 OS=Homo sapiens OX=9606 GN=ALDH3B1 PE=1 SV=1 E(blastx)=5e-54"
#[2] "Betaine--homocysteine S-methyltransferase 1 OS=Danio rerio OX=7955 GN=bhmt PE=2 SV=1 E(blastx)=1e-154" 
#[3] "Epidermal retinol dehydrogenase 2 OS=Mus musculus OX=10090 GN=Sdr16c5 PE=2 SV=1 E(blastx)=1e-100"  

#create table for ms
library(gt)

purple.pos <-purple.pos_genes%>%
  rownames_to_column()%>%
  dplyr::rename("V1"="rowname")%>%
  full_join(purplePC1_pos_genes)%>%
  na.omit()%>%
  dplyr::rename("Isogroup"="V1")%>%
  dplyr::rename("Annotation"="V2")%>%
  dplyr::rename("LD1 Loading"="PCA-pa.1")

purple.neg <-purple.neg_genes%>%
  rownames_to_column()%>%
  dplyr::rename("V1"="rowname")%>%
  full_join(purplePC1_neg_genes)%>%
  na.omit()%>%
  dplyr::rename("Isogroup"="V1")%>%
  dplyr::rename("Annotation"="V2")%>%
  dplyr::rename("LD1 Loading"="PCA-pa.1")%>%
  full_join(purple.pos)

purple_table <- gt(purple.neg)
purple_table<-
  purple_table|>
  tab_row_group(label= "Positive", rows = 4:10)|>
  tab_row_group(label= "Negative",rows = 1:3)|>
  tab_style(
    style=cell_fill(color = "grey90"), locations = cells_row_groups(groups =c("Positive", "Negative"))
  )
purple_table

gtsave(purple_table, filename="purple_genes.pdf")


# next, looking at extracting plasticity and frontloading values from LD2, combine with symbiont df
# plasticity = distance between distributions for each genotype
# baseline = mean of the control distribution for each genotype
combined.purple.expression<- combined.purple.dapc2%>%select(grp,c1)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  ungroup()%>%
  select(treatment, genotype, c1)%>%
  unique()%>%
  spread(treatment, c1)%>%
  mutate(plasticity = (Combined-Control))%>%
  dplyr::rename("baseline"="Control")%>%
  dplyr::rename("Genotype"="genotype")%>%
  select(everything(), -"Combined")%>%
  full_join(combined_symbiont_traits, by = "Genotype")

#individual frag info
combined.purple.frag <-combined.purple.dapc2%>%select(grp,sample,LD1)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  filter(treatment == "Combined")%>%
  ungroup()%>%
  select(sample, genotype, LD1)%>%
  dplyr::rename("combined.expression"="LD1")%>%
  dplyr::rename("Genotype"='genotype')%>%
  full_join(combined.purple.expression, by = "Genotype")%>%
  select(sample, Genotype, combined.expression, baseline)%>%
  mutate(frag.plasticity = combined.expression-baseline)%>%
  full_join(combined_individual_traits, by = c("Genotype"="Genotype", "sample"="sample"))%>%
  drop_na()

### linear models

## symbiont difference - linear regression
combined.purple.baseline_diff<-lm(symb.diff~baseline, data = combined.purple.expression)
summary(combined.purple.baseline_diff)
# NS: p = 0.1346

plot(combined.purple.baseline_diff)

combined.purple.plasticity_diff<-lm(symb.diff~plasticity, data = combined.purple.expression)
summary(combined.purple.plasticity_diff)
#NS: p = 0.0843

#individual frag data
individual_purple_diff_pl<-lm(Sym.difference~frag.plasticity, data = combined.purple.frag)
summary(individual_purple_diff_pl)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept)       439344     362428   1.212  0.23175   
#frag.plasticity   102040      33888   3.011  0.00426 **

#Residual standard error: 366300 on 45 degrees of freedom
#Multiple R-squared:  0.1677,	Adjusted R-squared:  0.1492 
#F-statistic: 9.067 on 1 and 45 DF,  p-value: 0.004259

plot(individual_purple_diff_pl)

purple.pl.plot <- ggplot(combined.purple.frag, colour = Genotype, aes(x=frag.plasticity, y = Sym.difference))+
  geom_smooth(method = "lm", se = TRUE, colour = "grey60")+
  geom_point(shape = 21, size = 3, aes(x=frag.plasticity, y = Sym.difference, fill = Genotype))+
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_reverse()+
  xlab("Plasticity")+
  theme_bw(base_size = 14)+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
purple.pl.plot

individual_purple_diff_bl<-lm(~Genotype*Sym.difference, data = combined.purple.frag)
summary(individual_purple_diff_bl)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  -411696     122949  -3.349  0.00165 **
#  baseline      -44831      21498  -2.085  0.04273 * 

#Residual standard error: 383400 on 45 degrees of freedom
#Multiple R-squared:  0.08813,	Adjusted R-squared:  0.06786 
#F-statistic: 4.349 on 1 and 45 DF,  p-value: 0.04273

plot(individual_purple_diff_bl)

purple.bl.plot <- ggplot(combined.purple.frag, colour = Genotype, aes(x=baseline, y = Sym.difference))+
  geom_smooth(method = "lm", se = TRUE, colour = "grey60")+
  geom_point(shape = 21, size= 3, aes(x=baseline, y = Sym.difference, fill = Genotype))+
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_reverse()+
  ylab("Difference in Symbiont Density/cm2")+
  xlab("Baseline Expression")+
  theme_bw(base_size=14)+
  theme(legend.position="none")
purple.bl.plot

## normalized symbs - quasibinomial
combined.purple.baseline_norm<-glm(symb.norm~baseline, data = combined.purple.expression, family = quasibinomial(link = 'logit'))
summary(combined.purple.baseline_norm)
#NS: p = 0.106

combined.purple.plasticity_norm<-glm(symb.norm~plasticity, data = combined.purple.expression, family = quasibinomial(link = 'logit'))
summary(combined.purple.plasticity_norm)
#NS: p = 0.0667

# combine baseline and plasticity figs
purple.symb.plot <- plot_grid(purple.bl.plot, purple.pl.plot, labels=c("B","C"), rel_widths = c(0.9,1))
purple.symb.plot

#module purple plot
purple_figure <- plot_grid(purple_ridges, purple.symb.plot, nrow = 2, labels = c("A"), rel_widths = c(1, 0.9))
purple_figure

#### HEAT TREATMENT #####
# need to create a counts table and colData that only includes the control and heat treatment samples
colData<-colData%>%mutate_at("treatment", as.factor)
colDataHeat=colData[colData$treatment=="Control"|colData$treatment=="Heat",]
colDataHeat$group2=factor(paste0(colDataHeat$genotype,colDataHeat$treatment))

#create a new key with only the samples of interest
keyHeat<-key%>%filter(Treatment =="Control"|Treatment =="Heat")

#use this key to filter the counts table 
countsHeat<-counts[, which((names(counts) %in% keyHeat$Sample)==TRUE)]


#whole genome response to heat
dds<-DESeqDataSetFromMatrix(countsHeat,
                            colData = colDataHeat, 
                            design = ~group2)

vstHeat=vst(dds)
vstHeat.df = assay(vstHeat)

#how much of thi variance is explained by geno vs treatment?
adHeat = adonis(t(vstHeat.df)~ treatment+genotype, data = colDataHeat, method="manhattan")
adHeat$aov.tab

#Permutation: free
#Number of permutations: 999

#Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#treatment  1  414943292 414943292  9.1886 0.07596  0.001 ***
#  genotype  10 1073423372 107342337  2.3770 0.19651  0.001 ***
#  Residuals 88 3973941158  45158422         0.72752           
#Total     99 5462307823                   1.00000  

dim(vstHeat.df) #16988  100

dapc1 <- dapc(t(vstHeat.df), colDataHeat$group2)
## choose to retain 50 PCs, 6 discriminant functions
dapc1
## this explains 92.4% of the conserved variance

## for the vst, they suggest retaining 17 PCs, which explains 72.2% of the conserved variance
temp<- optim.a.score(dapc1, n.sim = 1) 
heat.whole.dapc <- dapc(t(vstHeat.df), colDataHeat$group2, n.da=4, n.pca=17)
heat.whole.dapc

#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=FALSE,solid=.4)
#scatter(dapc,1,1,bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

heat.whole.varexpl <- round((heat.whole.dapc$eig/sum(heat.whole.dapc$eig))[1:2] * 100, 1) 

heat.whole.dapc1 <- tibble(sample = rownames(heat.whole.dapc$ind.coord),
                grp = heat.whole.dapc$grp,
                LD1 = heat.whole.dapc$ind.coord[,1],
                LD2 = heat.whole.dapc$ind.coord[,2])
heat.whole.dapc2 <- heat.whole.dapc1 %>%
  group_by(grp)%>%
  dplyr::summarize(c1 = mean(LD1),
                   c2 = mean(LD2), .groups="keep")%>%
  full_join(heat.whole.dapc1)

### DAPC FIGURES

heat.whole.dapc.fig <-   ggplot(heat.whole.dapc2, aes(fill = factor(str_sub(grp, 1,3)), 
                                 shape = factor(str_sub(grp, 4)))) +
  geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
  geom_point(aes(x = c1, y = c2), size = 3) +
  geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
  scale_shape_manual(name="Treatment", labels = c("Control","Heat"), values = c(21,23)) +
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_continuous(limits = c(-20,10)) +
  scale_y_continuous(limits = c(-10,10)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("LD1 [", heat.whole.varexpl[1],"%]"), y = paste0("LD2 [", heat.whole.varexpl[2],"%]")) +
  theme_bw()
heat.whole.dapc.fig

heat.ridges=ggplot(heat.whole.dapc2, aes(x = LD1, fill = factor(str_sub(grp, 4)), y = factor(str_sub(grp, 1,3)), height = ..density..)) +
  geom_density_ridges(scale = 2, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(limits = c(-12,15)) +
  scale_fill_manual(name="Treatment", values =c("#D3D3D3","#F88379")) +
  ylab("Genotype")+
  theme_ridges() +
  theme(legend.position = "top",
        legend.title=element_text(size=11),
        legend.text=element_text(size =10))
heat.ridges

# get gene loadings from LD1
heat.whole.geneloadings <-as.data.frame(heat.whole.dapc$pca.loadings)%>%
  select("PCA-pa.1")

#select top 20 on either end
heat.whole.pos_genes<-heat.whole.geneloadings%>%filter(.>0.02885172)
heat.whole.neg_genes<-heat.whole.geneloadings%>%filter(.< -0.01957875)

heat.whole.PC1_pos_genes=iso2gene[iso2gene$V1 %in% rownames(heat.whole.pos_genes),] 
heat.whole.PC1_pos_genes$V2
# 4 annotated
#[1] "40S ribosomal protein S15a OS=Rattus norvegicus OX=10116 GN=Rps15a PE=1 SV=2 E(blastx)=2e-76"
#[2] "Plasma alpha-L-fucosidase OS=Rattus norvegicus OX=10116 GN=Fuca2 PE=2 SV=1 E(blastx)=2e-175" 
#[3] "RING-box protein 1 OS=Salmo salar OX=8030 GN=rbx1 PE=2 SV=2 E(blastx)=2e-57"                 
#[4] "Transmembrane protein 138 OS=Bos taurus OX=9913 GN=TMEM138 PE=2 SV=1 E(blastx)=5e-48"  

heat.whole.PC1_neg_genes=iso2gene[iso2gene$V1 %in% rownames(heat.whole.neg_genes),] 
heat.whole.PC1_neg_genes$V2
# 2 annotated
#[1] "Actin-101 OS=Solanum tuberosum OX=4113 GN=AC101 PE=3 SV=1 E(blastx)=9e-87"                     
#[2] "E3 ubiquitin-protein ligase RNF220 OS=Homo sapiens OX=9606 GN=RNF220 PE=1 SV=1 E(blastx)=3e-58"

heat.pos <-heat.whole.pos_genes%>%
  rownames_to_column()%>%
  dplyr::rename("V1"="rowname")%>%
  full_join(heat.whole.PC1_pos_genes)%>%
  na.omit()%>%
  dplyr::rename("Isogroup"="V1")%>%
  dplyr::rename("Annotation"="V2")%>%
  dplyr::rename("LD1 Loading"="PCA-pa.1")

heat.neg <-heat.whole.neg_genes%>%
  rownames_to_column()%>%
  dplyr::rename("V1"="rowname")%>%
  full_join(heat.whole.PC1_neg_genes)%>%
  na.omit()%>%
  dplyr::rename("Isogroup"="V1")%>%
  dplyr::rename("Annotation"="V2")%>%
  dplyr::rename("LD1 Loading"="PCA-pa.1")%>%
  full_join(heat.pos)

heat_table <- gt(heat.neg)
heat_table<-
  heat_table|>
  tab_row_group(label= "Positive", rows = 3:6)|>
  tab_row_group(label= "Negative",rows = 1:2)|>
  tab_style(
    style=cell_fill(color = "grey90"), locations = cells_row_groups(groups =c("Positive", "Negative"))
  )
heat_table

gtsave(heat_table, filename="heat_genes.pdf")

### creating dataframe with the physiological data from combined vs. control treatments
##making experimental trait df 
heat.traits<-read.csv("Datasheets/Phenotype_metrics_measured.csv",header=T)%>%
  dplyr::rename("Genotype"="Genotype..")%>%
  mutate(sample = paste0(Genotype,sep="_",Tank))%>%
  select(-Tank)%>%
  filter(Treatment=="Control"|Treatment=="High Temperature")

#rename genotypes
heat.traits$Genotype <- sub("G63", "G62", heat.traits$Genotype)
heat.traits$Genotype <- sub("G3.b", "G3", heat.traits$Genotype)
heat.traits$Genotype <- sub("G34.21", "G34", heat.traits$Genotype)
heat.traits$Genotype <- sub("G44.b", "G44", heat.traits$Genotype)
heat.traits$Genotype <- sub("G44.20", "G44", heat.traits$Genotype)
heat.traits$Genotype <- sub("G62.21", "G62", heat.traits$Genotype)
heat.traits$Genotype <- sub("^G1$", "G01", heat.traits$Genotype)
heat.traits$Genotype <- sub("^G3$", "G03", heat.traits$Genotype)
heat.traits$Genotype <- sub("^G7$", "G07", heat.traits$Genotype)
heat.traits$Treatment <- sub("High Temperature", "Heat", heat.traits$Treatment)

# summarize trait data by genotype in the heat treatment
heat_symbiont_traits<-heat.traits%>%
  select(Genotype,Treatment,Symbionts.cm2)%>%
  group_by(Genotype,Treatment)%>%
  dplyr::summarize(mean=mean(Symbionts.cm2, na.rm=TRUE), .groups="keep")%>%
  spread(Treatment, mean)%>%
  mutate(symb.norm=Heat/Control)%>%
  mutate(symb.diff=Heat-Control)%>%
  select(Genotype, symb.norm, symb.diff)
# can't use binomial models because the normalized values are above 1

# keeping trait data unique to fragment level, but normalizing by the mean of the control.
heat_individual_traits<-heat.traits%>%select(Genotype,Treatment,Symbionts.cm2,sample)%>%
  filter(Treatment == "Heat")%>%
  spread(Treatment, Symbionts.cm2)%>%
  group_by(Genotype)%>%
  full_join(combined_symbiont_traits, by = "Genotype")%>%
  select(sample, Genotype, Heat, Control)%>%
  dplyr::rename("Control.mean"="Control")%>%
  mutate(Symb.normalized= Heat/Control.mean)%>%
  mutate(Sym.difference = Heat-Control.mean)

# next, looking at extracting plasticity and frontloading values from LD1, combine with symbiont df
# plasticity = distance between distributions for each genotype
# baseline = mean of the control distribution for each genotype
heat.whole.expression<- heat.whole.dapc2%>%
  select(grp,c1)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  ungroup()%>%
  select(treatment, genotype, c1)%>%
  unique()%>%
  spread(treatment, c1)%>%
  mutate(plasticity = (Heat-Control))%>%
  dplyr::rename("baseline"="Control")%>%
  dplyr::rename("Genotype"="genotype")%>%
  select(everything(), -"Heat")%>%
  full_join(heat_symbiont_traits, by = "Genotype")

#individual frag info
heat.whole.frag <-heat.whole.dapc2%>%select(grp,sample,LD1)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  filter(treatment == "Heat")%>%
  ungroup()%>%
  select(sample, genotype, LD1)%>%
  dplyr::rename("heat.expression"="LD1")%>%
  dplyr::rename("Genotype"='genotype')%>%
  full_join(heat.whole.expression, by = "Genotype")%>%
  select(sample, Genotype, heat.expression, baseline)%>%
  mutate(frag.plasticity = heat.expression-baseline)%>%
  full_join(heat_individual_traits, by = c("Genotype"="Genotype", "sample"="sample"))%>%
  drop_na()

### linear models

## symbiont difference - linear regression
heat.whole.baseline_diff<-lm(symb.diff~baseline, data = heat.whole.expression)
summary(heat.whole.baseline_diff)
# NS: p = 0.0.2130

heat.whole.plasticity_diff<-lm(symb.diff~plasticity, data = heat.whole.expression)
summary(heat.whole.plasticity_diff)
#NS: p = 0.399

#individual frag data
individual_heat_diff_pl<-lm(Sym.difference~frag.plasticity, data = heat.whole.frag)
summary(individual_heat_diff_pl)

plot(individual_heat_diff_pl)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       142402     234582   0.607   0.5470  
#frag.plasticity   -66482      30388  -2.188   0.0342 *

#Residual standard error: 369400 on 43 degrees of freedom
#Multiple R-squared:  0.1002,	Adjusted R-squared:  0.07924 
#F-statistic: 4.787 on 1 and 43 DF,  p-value: 0.03417

heat.pl.plot <- ggplot(heat.whole.frag, colour = Genotype, aes(x=frag.plasticity, y = Sym.difference))+
  geom_smooth(method = "lm", se = TRUE, colour = "grey60")+
  geom_point(shape = 21, size = 3, aes(x=frag.plasticity, y = Sym.difference, fill = Genotype))+
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  ylab("Difference in Symbiont Density/cm2")+
  xlab("Plasticity")+
  theme_bw(base_size=14)
heat.pl.plot

individual_heat_diff_bl<-lm(Sym.difference~baseline, data = heat.whole.frag)
summary(individual_heat_diff_bl)
#NS: 0.088420

### composite plot for whole genome response in heat
whole_heat_fig <- plot_grid(heat.ridges, heat.pl.plot, labels = "AUTO", rel_widths =c(0.7,1))
whole_heat_fig

### MIDNIGHTBLUE MODULE ######
## first, load in the gene list of midnightblue module output from WGCNA
modMblue<-read.csv("Host GO New Ref/unmerged_hostGO_MM_midnightblue_numeric.csv",header=T)%>%
  filter(kME > 0)

# now filter the counts data by isogroups that are present in the module
mblueCounts=counts[rownames(counts) %in% modMblue$ProbeID,]

# filter the purpleCounts to be only combined and control using keyComp
mblueHeat<-mblueCounts[, which((names(mblueCounts) %in% keyHeat$Sample)==TRUE)]

dds<-DESeqDataSetFromMatrix(mblueHeat,
                            colData = colDataHeat, 
                            design = ~group2)

vstMblue=varianceStabilizingTransformation(dds)
vstMblue.df = assay(vstMblue)

dim(vstMblue.df) #291  100

#how much variation in this module is explained by treatment vs geno?
adMblue = adonis(t(vstMblue.df)~ treatment+genotype, data = colDataHeat, method="manhattan")
adMblue$aov.tab

#Permutation: free
#Number of permutations: 999

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treatment  1    188166  188166 14.6482 0.11125  0.001 ***
#  genotype  10    372728   37273  2.9016 0.22038  0.001 ***
#  Residuals 88   1130417   12846         0.66837           
#Total     99   1691311                 1.00000          

dapcMb <- dapc(t(vstMblue.df), colDataHeat$group2)
## choose to retain 60 PCs, 5 discriminant functions
dapcMb
## this explains 96.5% of the conserved variance

## for the vst, they suggest retaining 14 PCs, which explains 68.1% of the conserved variance
temp <- optim.a.score(dapcMb, n.sim = 1) 
heat.mblue.dapc <- dapc(t(vstMblue.df), colDataHeat$group2, n.da=4, n.pca=14)
heat.mblue.dapc

#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=FALSE,solid=.4)
#scatter(dapc,1,1,bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

heat.mblue.varexpl <- round((heat.mblue.dapc$eig/sum(heat.mblue.dapc$eig))[1:2] * 100, 1) 

heat.mblue.dapc1 <- tibble(sample = rownames(heat.mblue.dapc$ind.coord),
                           grp = heat.mblue.dapc$grp,
                           LD1 = heat.mblue.dapc$ind.coord[,1],
                           LD2 = heat.mblue.dapc$ind.coord[,2])
heat.mblue.dapc2 <- heat.mblue.dapc1 %>%
  group_by(grp)%>%
  dplyr::summarize(c1 = mean(LD1),
                   c2 = mean(LD2), .groups="keep")%>%
  full_join(heat.mblue.dapc1)

## visualization
mblue.dapc<-ggplot(heat.mblue.dapc2, aes(fill = factor(str_sub(grp, 1,3)), 
                                         shape = factor(str_sub(grp, 4)))) +
  geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
  geom_point(aes(x = c1, y = c2), size = 3) +
  geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
  scale_shape_manual(name="Treatment", labels = c("Control","Heat"), values = c(21,23)) +
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  scale_x_continuous(limits = c(-20,10)) +
  scale_y_continuous(limits = c(-10,10)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("LD1 [", heat.mblue.varexpl[1],"%]"), y = paste0("LD2 [", heat.mblue.varexpl[2],"%]")) +
  theme_bw()
mblue.dapc

mblue.ridges=ggplot(heat.mblue.dapc2, aes(x = LD2, fill = factor(str_sub(grp, 4)), y = factor(str_sub(grp, 1,3)), height = ..density..)) +
  geom_density_ridges(scale = 1.5, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_reverse(limits = c(10,-12)) +
  scale_fill_manual(name="Treatment", values =c("#D3D3D3","#F88379")) +
  ylab("Genotype")+
  theme_ridges() +
  theme(legend.position = "none")
mblue.ridges

#pulling out sig gene loadings for PC2 (treatment axis) in the yellow module
geneloadings_mblue <-as.data.frame(heat.mblue.dapc$pca.loadings)%>%
  select("PCA-pa.2")

#select top 20 on either end
mblue.pos_genes<-geneloadings_mblue%>%filter(.>	0.0733873634)
mblue.neg_genes<-geneloadings_mblue%>%filter(.< -0.093564680)

mbluePC2_pos_genes=iso2gene[iso2gene$V1 %in% rownames(mblue.pos_genes),] 
mbluePC2_pos_genes$V2
# 4 annotated
#[1] "Tyrosine-protein kinase transforming protein Abl OS=Feline sarcoma virus (strain Hardy-Zuckerman 2) OX=11776 GN=ABL PE=2 SV=1 E(blastx)=6e-29"
#[2] "Extracellular superoxide dismutase [Cu-Zn] OS=Drosophila melanogaster OX=7227 GN=Sod3 PE=1 SV=1 E(blastx)=1e-06"                              
#[3] "2-5-oligoadenylate synthase 1 OS=Sus scrofa OX=9823 GN=OAS1 PE=1 SV=3 E(blastx)=2e-46"                                                        
#[4] "Transcription factor MafA OS=Xenopus tropicalis OX=8364 GN=mafa PE=2 SV=1 E(blastx)=2e-12"  

mbluePC2_neg_genes=iso2gene[iso2gene$V1 %in% rownames(mblue.neg_genes),] 
mbluePC2_neg_genes$V2
# 4 annotated
#[1] "E3 ubiquitin-protein ligase SH3RF1 OS=Xenopus tropicalis OX=8364 GN=sh3rf1 PE=2 SV=1 E(blastx)=3e-23"    
#[2] "Ankyrin repeat domain-containing protein 50 OS=Homo sapiens OX=9606 GN=ANKRD50 PE=1 SV=4 E(blastx)=6e-32"
#[3] "Contactin-associated protein-like 2 OS=Mus musculus OX=10090 GN=Cntnap2 PE=1 SV=2 E(blastx)=2e-75"       
#[4] "Syntaxin-17 OS=Homo sapiens OX=9606 GN=STX17 PE=1 SV=2 E(blastx)=4e-38" 

# next, looking at extracting plasticity and frontloading values from LD2, combine with symbiont df
# plasticity = distance between distributions for each genotype
# baseline = mean of the control distribution for each genotype
heat.mblue.expression<- heat.mblue.dapc2%>%
  select(grp,c2)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  ungroup()%>%
  select(treatment, genotype, c2)%>%
  unique()%>%
  spread(treatment, c2)%>%
  mutate(plasticity = (Heat-Control))%>%
  dplyr::rename("baseline"="Control")%>%
  dplyr::rename("Genotype"="genotype")%>%
  select(everything(), -"Heat")%>%
  full_join(heat_symbiont_traits, by = "Genotype")

#individual frag info
heat.mblue.frag <-heat.mblue.dapc2%>%select(grp,sample,LD2)%>%
  mutate(genotype = factor(str_sub(grp,1,3)))%>%
  mutate(treatment = factor(str_sub(grp,4)))%>%
  filter(treatment == "Heat")%>%
  ungroup()%>%
  select(sample, genotype, LD2)%>%
  dplyr::rename("heat.expression"="LD2")%>%
  dplyr::rename("Genotype"='genotype')%>%
  full_join(heat.mblue.expression, by = "Genotype")%>%
  select(sample, Genotype, heat.expression, baseline)%>%
  mutate(frag.plasticity = heat.expression-baseline)%>%
  full_join(heat_individual_traits, by = c("Genotype"="Genotype", "sample"="sample"))%>%
  drop_na()

### linear models

## symbiont difference - linear regression
heat.mblue.baseline_diff<-lm(symb.diff~baseline, data = heat.mblue.expression)
summary(heat.mblue.baseline_diff)
# NS: p = 0.9721

heat.mblue.plasticity_diff<-lm(symb.diff~plasticity, data = heat.mblue.expression)
summary(heat.mblue.plasticity_diff)
#NS: p = 0.0837

#individual frag data
individual_mblue_diff_pl<-lm(Sym.difference~frag.plasticity, data = heat.mblue.frag)
summary(individual_mblue_diff_pl)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)        84802     148822   0.570  0.57176   
#frag.plasticity    87517      27635   3.167  0.00283 **

#Residual standard error: 350600 on 43 degrees of freedom
#Multiple R-squared:  0.1891,	Adjusted R-squared:  0.1703 
#F-statistic: 10.03 on 1 and 43 DF,  p-value: 0.002832

plot(individual_mblue_diff_pl)

mblue.pl.plot <- ggplot(heat.mblue.frag, colour = Genotype, aes(x=frag.plasticity, y = Sym.difference))+
  geom_smooth(method = "lm", se = TRUE, colour = "grey60")+
  geom_point(shape = 21, aes(x=frag.plasticity, y = Sym.difference, fill = Genotype))+
  scale_fill_brewer(name ="Genotype", palette = "Set3")+
  ylab("Difference in Symbiont Density/cm2")+
  xlab("Plasticity")+
  theme_bw()
mblue.pl.plot

individual_mblue_diff_bl<-lm(Sym.difference~baseline, data = heat.mblue.frag)
summary(individual_mblue_diff_bl)
#NS: p = 0.421312
