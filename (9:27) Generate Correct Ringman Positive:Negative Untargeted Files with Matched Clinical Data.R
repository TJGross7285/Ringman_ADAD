############################### (9/27) Generate Correct Ringman Positive/Negative Untargeted Files with Matched Clinical Data

library(Biobase)
library(sva)
library(WGCNA)


####Read in negative mode data
negative<-read.csv("20150513_RINGMAN_BEH_AMIDE_NEG_XCMSv1.csv",header=FALSE)
rt_neg<-as.numeric(as.character(negative[3:dim(negative)[1],3]))
mz_neg<-as.numeric(as.character(negative[3:dim(negative)[1],2]))
neg_SampleID<-t(negative[1,4:dim(negative)[2]])
neg_abunds<-as.data.frame(t(negative[3:dim(negative)[1],4:dim(negative)[2]]))
colnames(neg_abunds)<-mz_neg
negativeIDs<-cbind(as.data.frame(neg_SampleID),as.data.frame(neg_abunds))
colnames(negativeIDs)[1]<-"SampleID"

####Read in positive mode data 
positive<-read.csv("20150513_RINGMAN_BEH_AMIDE_POS_XCMSv1.csv",header=FALSE)
rt_pos<-as.numeric(as.character(positive[3:dim(positive)[1],3]))
mz_pos<-as.numeric(as.character(positive[3:dim(positive)[1],2]))
pos_SampleID<-t(positive[1,4:dim(positive)[2]])
pos_abunds<-as.data.frame(t(positive[3:dim(positive)[1],4:dim(positive)[2]]))
colnames(pos_abunds)<-mz_pos
positiveIDs<-cbind(as.data.frame(pos_SampleID),as.data.frame(pos_abunds))
colnames(positiveIDs)[1]<-"SampleID"

####Read in clinical data, join to abundances, and write to RData file
clinical<-read.csv("Revised_Ringman_Clinical_7_14.csv",check.names=FALSE)
positive_final<-dplyr::inner_join(clinical,positiveIDs)
colnames(positive_final)[2]<-"DiseaseState"
negative_final<-dplyr::inner_join(clinical,negativeIDs)
colnames(negative_final)[2]<-"DiseaseState"
levels(positive_final$DiseaseState) <- list("MCI"=c("MCI"), "AD"=c("Dementia"),"Carrier"=c("Asymptomatic"),"Control"=c("Noncarrier control"))
levels(negative_final$DiseaseState) <- list("MCI"=c("MCI"), "AD"=c("Dementia"),"Carrier"=c("Asymptomatic"),"Control"=c("Noncarrier control"))


identical(positive_final[,c(1:11)],negative_final[,c(1:11)])

pheno<-positive_final$DiseaseState
levels(pheno)<-list("Control"=c("Noncarrier control"),"Carrier"=c("Asymptomatic"),"MCI"=c("MCI"),"AD"=c("AD"))
pheno_meta<-cbind(as.data.frame(pheno),as.data.frame(positive_final[,1:11]))
colnames(pheno_meta)[1]<-"ClinStatus"
colnames(pheno_meta)[9:10]<-c("CDR_Global","CDR_Sum")

pheno_meta<-AnnotatedDataFrame(pheno_meta)

abunds<-apply(cbind(as.data.frame(positive_final[,-c(1:11)]),as.data.frame(negative_final[,-c(1:11)])),2,as.numeric)
abunds<-log2(abunds)

write.csv(abunds,file="Ringman_Abunds.csv")

abunds<-read.csv("Ringman_Abunds.csv",check.names=FALSE,na.strings="0")[,-c(1)]
colnames(abunds)<-seq(1,dim(abunds)[2])
index<-caret::nearZeroVar(abunds)
abunds<-abunds[,-index]
fix<-apply(abunds,2,is.infinite)
abunds[fix==TRUE]<-NA


mode<-as.factor(c(rep("ESI+",dim(positive_final[,-c(1:11)])[2]),rep("ESI-",dim(negative_final[,-c(1:11)])[2])))
feature_meta<-cbind(rbind(cbind(mz_pos,rt_pos),cbind(mz_neg,rt_neg)),as.data.frame(as.factor(mode)))
feature_meta<-feature_meta[-index,]
colnames(feature_meta)<-c("MZ","RT","Mode")
rownames(feature_meta)<-seq(1,dim(abunds)[2])

feature_meta_use<-cbind(as.data.frame(seq(1,dim(abunds)[2])),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

feature_meta<-AnnotatedDataFrame(feature_meta)



####Impute missing data with KNN imputation
library(impute)
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
rownames(edata)<-seq(1,dim(abunds)[2])

combined<-ExpressionSet(assayData=edata,phenoData=pheno_meta,featureData=feature_meta)

mod<-model.matrix(~as.factor(ClinStatus)+Gene,data=combined)
mod0<-model.matrix(~Gene,data=combined)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)
design_table<-cbind(model.matrix(~0+as.factor(ClinStatus),data=combined),as.data.frame(svobj$sv))
colnames(design_table)<-c("Control","Carrier","MCI","AD","SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10","SV11","SV12")

library(limma)
library(dplyr)
####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design_table)
cm <- makeContrasts(
	`Carrier-Control` = Carrier-Control,
	`AD-MCI` = AD-MCI,
	`MCI-Carrier`= MCI-Carrier, 
	levels=design_table)


cm1<- makeContrasts(
	`AD-MCI` = AD-MCI, 
	levels=design_table)
cm2<-makeContrasts(
	`Carrier-Control` = Carrier-Control, 
	levels=design_table)
cm3<-makeContrasts(
	`MCI-Carrier`= MCI-Carrier, 
	levels=design_table)



fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"

fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)
T<-topTableF(fit1_F,adjust="BH",number=100000)
U<-cbind(rownames(T),T)
colnames(U)[1]<-"Feature"

fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)
T<-topTableF(fit3_F,adjust="BH",number=100000)
V<-cbind(rownames(T),T)
colnames(V)[1]<-"Feature"





























fi

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"

feature_meta_use$Feature<-as.factor(feature_meta_use$Feature)
join<-dplyr::inner_join(as.data.frame(feature_meta_use),T,by="Feature")%>%arrange(P.Value)

lip_neg_abunds<-read.csv("mapstone_r56_lipidomics_NEG_sampleids.csv", header=FALSE)[-c(1),-c(1:2)]
lip_pos_abunds<-read.csv("mapstone_r56_lipidomics_POS_sampleids.csv", header=FALSE)[-c(1),-c(1:2)]
metab_neg_abunds<-read.csv("mapstone_r56_metabolomics_NEG_sampleids.csv", header=FALSE)[-c(1),-c(1:2)]
metab_pos_abunds<-read.csv("mapstone_r56_metabolomics_POS_sampleids.csv", header=FALSE)[-c(1),-c(1:2)]















fin<-cbind(as.data.frame(rownames(fData(combined))),as.data.frame(fData(combined)))
colnames(fin)[1]<-"Feature"
Carrier_Control_UP<-T%>%filter(P.Value<.05 & Carrier.Control >.1)%>%select(Feature,Carrier.Control)%>%inner_join(fin,by="Feature")
Carrier_Control_DOWN<-T%>%filter(P.Value<.05 & Carrier.Control < -.1)%>%select(Feature,Carrier.Control)%>%inner_join(fin,by="Feature")
MCIAD_Carrier_UP<-T%>%filter(P.Value<.05 & MCIAD.Carrier > .1)%>%select(Feature,MCIAD.Carrier)%>%inner_join(fin,by="Feature")
MCIAD_Carrier_DOWN<-T%>%filter(P.Value<.05 & MCIAD.Carrier < -.1)%>%select(Feature,MCIAD.Carrier)%>%inner_join(fin,by="Feature")

overlap_UP<-intersect(Carrier_Control_UP$Feature,MCIAD_Carrier_UP$Feature)
overlap_DOWN<-intersect(Carrier_Control_DOWN$Feature,MCIAD_Carrier_DOWN$Feature)
overlap_CROSS<-c(intersect(Carrier_Control_DOWN$Feature,MCIAD_Carrier_UP$Feature),intersect(Carrier_Control_UP$Feature,MCIAD_Carrier_DOWN$Feature))

UP_index<-T$Feature %in% overlap_UP
UP<-fin[UP_index==TRUE,]

#####NO FEATURES!!!!!!!
DOWN_index<-T$Feature %in% overlap_DOWN
DOWN<-fin[DOWN_index==TRUE,]

CROSS_index<-T$Feature %in% overlap_CROSS
CROSS<-fin[CROSS_index==TRUE,]

pdf(file="PACE_DE_Collapsed_Venn_Ringman_6_25.pdf")

GOVenn(list1,list2,label=c("Carrier-Control",
								"MCIAD-Carrier"),
								lfc.col=c("seagreen","gold","firebrick2"))

dev.off()











colnames(table)<-c("p-Value","q-Value","MZ","RT","Mode","Log2(MCIAD/Carrier)","Log2(Carrier/Control)","Log2(MCIAD/Control)")

write.csv(table,file="Ringman_DE.csv")

































####Get non carrier control negative mode
noncarrier_NEG<-dplyr::filter(negative_final,DiseaseState=="Noncarrier Control")

####Get MCIAD negative mode
mciad_NEG<-dplyr::filter(negative_final,DiseaseState=="MCIAD")

####Get Asymptomatic negative mode
asym_NEG<-dplyr::filter(negative_final,DiseaseState=="Asymptomatic")

####Get noncarrier control positive mode
noncarrier_POS<-dplyr::filter(positive_final,DiseaseState=="Noncarrier Control")

####Get MCIAD positive mode
mciad_POS<-dplyr::filter(positive_final,DiseaseState=="MCIAD")

####Get Asymptomatic positive mode
asym_POS<-dplyr::filter(positive_final,DiseaseState=="Asymptomatic")

####Write positive mode metadata 
pos_metadata<-t(rbind(rt_pos,mz_pos))

####Write negative mode metadata 
neg_metadata<-t(rbind(rt_neg,mz_neg))

write.csv(noncarrier_NEG, file="noncarrier_NEG.csv")
write.csv(mciad_NEG, file="mciad_NEG.csv")
write.csv(asym_NEG, file="asym_NEG.csv")
write.csv(noncarrier_POS, file="noncarrier_POS.csv")
write.csv(mciad_POS,file="mciad_POS.csv")
write.csv(asym_POS,file="asym_POS.csv")
write.csv(pos_metadata,file="POS_metadata.csv")
write.csv(neg_metadata, file="NEG_metadata.csv")

###################################
###################################
setwd("~/Desktop/Ringman")
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
mciad_NEG<-read.csv("mciad_NEG.csv",check.names=FALSE)
asym_NEG<-read.csv("asym_NEG.csv",check.names=FALSE)
noncarrier_POS<-read.csv("noncarrier_POS.csv",check.names=FALSE)
mciad_POS<-read.csv("mciad_POS.csv",check.names=FALSE)
asym_POS<-read.csv("asym_POS.csv",check.names=FALSE)
neg_metadata<-read.csv("NEG_metadata.csv")

####MCIAD vs Asymptomatic contrast negative mode
setwd("~/Desktop/Mummichog_Input")
neg_asym_mciad<-rbind(asym_NEG, mciad_NEG)
pheno<-as.factor(neg_asym_mciad$DiseaseState)
abunds<-neg_asym_mciad[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~pheno, exact=FALSE)$p.val)
median_case<-apply(abunds[pheno=="MCIAD",],2,median)
median_control<-apply(abunds[pheno=="Asymptomatic",],2,median)
log2FC<-median_case-median_control
mz<-colnames(abunds)
complete<-cbind(as.data.frame(neg_metadata$mz_neg),as.data.frame(neg_metadata$rt_neg),as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","rt","p-value","statistic")
write.table(complete, sep = "\t", file="NEG_asym_mciad.txt")

###################################
###################################
setwd("~/Desktop/Ringman")
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
mciad_NEG<-read.csv("mciad_NEG.csv",check.names=FALSE)
asym_NEG<-read.csv("asym_NEG.csv",check.names=FALSE)
noncarrier_POS<-read.csv("noncarrier_POS.csv",check.names=FALSE)
mciad_POS<-read.csv("mciad_POS.csv",check.names=FALSE)
asym_POS<-read.csv("asym_POS.csv",check.names=FALSE)
neg_metadata<-read.csv("NEG_metadata.csv")

####Asymptomatic vs non carrier contrast negative mode 
setwd("~/Desktop/Mummichog_Input") 
neg_asym_noncar<-rbind(asym_NEG, noncarrier_NEG)
pheno<-as.factor(neg_asym_noncar$DiseaseState)
abunds<-neg_asym_noncar[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~pheno, exact=FALSE)$p.val)
median_case<-apply(abunds[pheno=="Asymptomatic",],2,median)
median_control<-apply(abunds[pheno=="Noncarrier Control",],2,median)
log2FC<-median_case-median_control
mz<-colnames(abunds)
complete<-cbind(as.data.frame(neg_metadata$mz_neg),as.data.frame(neg_metadata$rt_neg), as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","rt","p-value","statistic")
write.table(complete, sep = "\t", file="NEG_asym_noncar.txt")

###################################
###################################
setwd("~/Desktop/Ringman")
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
mciad_NEG<-read.csv("mciad_NEG.csv",check.names=FALSE)
asym_NEG<-read.csv("asym_NEG.csv",check.names=FALSE)
noncarrier_POS<-read.csv("noncarrier_POS.csv",check.names=FALSE)
mciad_POS<-read.csv("mciad_POS.csv",check.names=FALSE)
asym_POS<-read.csv("asym_POS.csv",check.names=FALSE)
pos_metadata<-read.csv("POS_metadata.csv")


####MCIAD vs Asymptomatic contrast positive mode
setwd("~/Desktop/Mummichog_Input")
pos_asym_mciad<-rbind(asym_POS, mciad_POS)
pheno<-as.factor(pos_asym_mciad$DiseaseState)
abunds<-pos_asym_mciad[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~pheno, exact=FALSE)$p.val)
median_case<-apply(abunds[pheno=="MCIAD",],2,median)
median_control<-apply(abunds[pheno=="Asymptomatic",],2,median)
log2FC<-median_case-median_control
mz<-gsub("_.+$","", colnames(abunds))
complete<-cbind(as.data.frame(pos_metadata$mz_pos),as.data.frame(pos_metadata$rt_pos), as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","rt","p-value","statistic")
write.table(complete, sep = "\t", file="POS_asym_mciad.txt")

###################################
###################################
setwd("~/Desktop/Ringman")
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
mciad_NEG<-read.csv("mciad_NEG.csv",check.names=FALSE)
asym_NEG<-read.csv("asym_NEG.csv",check.names=FALSE)
noncarrier_POS<-read.csv("noncarrier_POS.csv",check.names=FALSE)
mciad_POS<-read.csv("mciad_POS.csv",check.names=FALSE)
asym_POS<-read.csv("asym_POS.csv",check.names=FALSE)
pos_metadata<-read.csv("POS_metadata.csv")


####Asymptomatic vs non carrier contrast positive mode 
setwd("~/Desktop/Mummichog_Input") 
pos_asym_noncar<-rbind(asym_POS, noncarrier_POS)
pheno<-as.factor(pos_asym_noncar$DiseaseState)
abunds<-pos_asym_noncar[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~pheno, exact=FALSE)$p.val)
median_case<-apply(abunds[pheno=="Asymptomatic",],2,median)
median_control<-apply(abunds[pheno=="Noncarrier Control",],2,median)
log2FC<-median_case-median_control
mz<-colnames(abunds)
complete<-cbind(as.data.frame(pos_metadata$mz_pos),as.data.frame(pos_metadata$rt_pos), as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","rt","p-value","statistic")
write.table(complete, sep = "\t", file="POS_asym_noncar.txt")

###################################
###################################
setwd("~/Desktop/Ringman")
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
mciad_NEG<-read.csv("mciad_NEG.csv",check.names=FALSE)
asym_NEG<-read.csv("asym_NEG.csv",check.names=FALSE)
noncarrier_POS<-read.csv("noncarrier_POS.csv",check.names=FALSE)
mciad_POS<-read.csv("mciad_POS.csv",check.names=FALSE)
asym_POS<-read.csv("asym_POS.csv",check.names=FALSE)
neg_metadata<-read.csv("NEG_metadata.csv")


####MCIAD vs non carrier contrast negative mode
setwd("~/Desktop/Mummichog_Input")
neg_noncar_mciad<-rbind(noncarrier_NEG, mciad_NEG)
pheno<-as.factor(neg_noncar_mciad$DiseaseState)
abunds<-neg_noncar_mciad[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~pheno, exact=FALSE)$p.val)
median_case<-apply(abunds[pheno=="MCIAD",],2,median)
median_control<-apply(abunds[pheno=="Noncarrier Control",],2,median)
log2FC<-median_case-median_control
mz<-gsub("_.+$","", colnames(abunds))
complete<-cbind(as.data.frame(neg_metadata$mz_neg),as.data.frame(neg_metadata$rt_neg), as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","rt","p-value","statistic")
write.table(complete, sep = "\t", file="NEG_noncar_mciad.txt")

###################################
###################################
setwd("~/Desktop/Ringman")
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
noncarrier_NEG<-read.csv("noncarrier_NEG.csv",check.names=FALSE)
mciad_NEG<-read.csv("mciad_NEG.csv",check.names=FALSE)
asym_NEG<-read.csv("asym_NEG.csv",check.names=FALSE)
noncarrier_POS<-read.csv("noncarrier_POS.csv",check.names=FALSE)
mciad_POS<-read.csv("mciad_POS.csv",check.names=FALSE)
asym_POS<-read.csv("asym_POS.csv",check.names=FALSE)
pos_metadata<-read.csv("POS_metadata.csv")


####MCIAD vs non carrier contrast positive mode
setwd("~/Desktop/Mummichog_Input")
pos_noncar_mciad<-rbind(noncarrier_POS, mciad_POS)
pheno<-as.factor(pos_noncar_mciad$DiseaseState)
abunds<-pos_noncar_mciad[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~pheno, exact=FALSE)$p.val)
median_case<-apply(abunds[pheno=="MCIAD",],2,median)
median_control<-apply(abunds[pheno=="Noncarrier Control",],2,median)
log2FC<-median_case-median_control
mz<-gsub("_.+$","", colnames(abunds))
complete<-cbind(as.data.frame(pos_metadata$mz_pos),as.data.frame(pos_metadata$rt_pos), as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","rt","p-value","statistic")
write.table(complete, sep = "\t", file="POS_noncar_mciad.txt")
























































save(rt_pos,rt_neg, file="Ringman_RT_9_27.RData")
asym_pos<-positive_final[positive_final$ClinicalGrouping=="Asymptomatic",]
asym_neg<-negative_final[negative_final$ClinicalGrouping=="Asymptomatic",]
write.csv(asym_pos, file="POS_Ringman_Final.csv")
write.csv(asym_neg, file="NEG_Ringman_Final.csv")





















####Load in processed Ringman data and subset to include only asymptomatic patients
load("Ringman_RT_9_27.RData",verbose=TRUE)
pos<-read.csv("POS_Ringman_Final.csv",check.names=FALSE)
neg<-read.csv("NEG_Ringman_Final.csv",check.names=FALSE)

####Visualize distribution of time relative to parental Dx *****MEDIAN = -13.5 years; Identical == TRUE 
hist(asym_pos$`Age Rel to PDx`, breaks=20)
median(asym_pos$`Age Rel to PDx`)
identical(asym_pos$`Age Rel to PDx`,asym_neg$`Age Rel to PDx`)

####Set up median split phenotypic variable for statistical contrast 
key<-pos$`Age Rel to PDx`> median(pos$`Age Rel to PDx`)
median_split<-as.factor(ifelse(key, yes="Soon", no="Later"))

####Soon-Later contrast for positive mode Ringman data 
abunds<-pos[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~median_split, exact=FALSE)$p.val)
median_case<-apply(abunds[median_split=="Soon",],2,median)
median_control<-apply(abunds[median_split=="Later",],2,median)
log2FC<-median_case-median_control
complete<-cbind(as.data.frame(colnames(abunds)),as.data.frame(rt_pos),as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","retention time","p-value","statistic score")
write.table(complete, sep = "\t", file="POS_Ringman_Soon_Later.txt")

####Soon-Later contrast for negative mode Ringman data 
abunds<-neg[,-c(1:11)]
stat<-apply(abunds,2,function(x) wilcox.test(x~median_split, exact=FALSE)$p.val)
median_case<-apply(abunds[median_split=="Soon",],2,median)
median_control<-apply(abunds[median_split=="Later",],2,median)
log2FC<-median_case-median_control
complete<-cbind(as.data.frame(colnames(abunds)),as.data.frame(rt_pos),as.data.frame(stat),as.data.frame(log2FC))
rownames(complete)<-NULL
colnames(complete)<-c("m/z","retention time","p-value","statistic score")
write.table(complete, sep = "\t", file="NEG_Ringman_Soon_Later.txt")