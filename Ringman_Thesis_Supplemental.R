library(Biobase)
library(sva)
library(impute)
library(limma)
library(dplyr)
library(GridOnClusters)
library(bnlearn)

###########START HERE TO REPRODUCE FINDINGS!!!
####Read in negative mode data
negative<-read.csv("20150513_RINGMAN_BEH_AMIDE_NEG_XCMSv1.csv",header=FALSE,na.strings="0")
rt_neg<-as.numeric(as.character(negative[3:dim(negative)[1],3]))
mz_neg<-as.numeric(as.character(negative[3:dim(negative)[1],2]))
neg_SampleID<-t(negative[1,4:dim(negative)[2]])
neg_abunds<-as.data.frame(t(negative[3:dim(negative)[1],4:dim(negative)[2]]))
colnames(neg_abunds)<-mz_neg
negativeIDs<-cbind(as.data.frame(neg_SampleID),as.data.frame(neg_abunds))
colnames(negativeIDs)[1]<-"SampleID"

####Read in positive mode data 
positive<-read.csv("20150513_RINGMAN_BEH_AMIDE_POS_XCMSv1.csv",header=FALSE,na.strings="0")
rt_pos<-as.numeric(as.character(positive[3:dim(positive)[1],3]))
mz_pos<-as.numeric(as.character(positive[3:dim(positive)[1],2]))
pos_SampleID<-t(positive[1,4:dim(positive)[2]])
pos_abunds<-as.data.frame(t(positive[3:dim(positive)[1],4:dim(positive)[2]]))
colnames(pos_abunds)<-mz_pos
positiveIDs<-cbind(as.data.frame(pos_SampleID),as.data.frame(pos_abunds))
colnames(positiveIDs)[1]<-"SampleID"

clinical<-read.csv("Ring_Test_3_30.csv",check.names=FALSE)###########%>%filter(ClinicalGrouping!="Noncarrier control")
positive_final<-dplyr::inner_join(clinical,positiveIDs)
colnames(positive_final)[2]<-"DiseaseState"
negative_final<-dplyr::inner_join(clinical,negativeIDs)
colnames(negative_final)[2]<-"DiseaseState"
positive_final$DiseaseState<-as.factor(positive_final$DiseaseState)
negative_final$DiseaseState<-as.factor(negative_final$DiseaseState)
levels(positive_final$DiseaseState)<-list("Control"=c("Noncarrier control"),"Carrier"=c("Asymptomatic"),"MCI/AD"=c("MCI","Dementia"))
levels(negative_final$DiseaseState)<-list("Control"=c("Noncarrier control"),"Carrier"=c("Asymptomatic"),"MCI/AD"=c("MCI","Dementia"))

####Should evaluate TRUE 
identical(positive_final[,c(1:27)],negative_final[,c(1:27)])

pheno<-positive_final$DiseaseState
pheno_meta<-cbind(as.data.frame(pheno),as.data.frame(positive_final[,1:27]))
colnames(pheno_meta)[1]<-"ClinStatus"
colnames(pheno_meta)[10:11]<-c("CDR_Global","CDR_Sum")
pheno_plot<-pheno_meta
pheno_meta<-AnnotatedDataFrame(pheno_meta)


abunds<-apply(cbind(as.data.frame(positive_final[,-c(1:27)]),as.data.frame(negative_final[,-c(1:27)])),2,as.numeric)
colnames(abunds)<-seq(1,dim(abunds)[2])
index<-caret::nearZeroVar(abunds)
T_abunds<-t(abunds[,-index])
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
edata<-log2(edata)
rownames(edata)<-rownames(T_abunds)

mode<-as.factor(c(rep("positive",dim(positive_final[,-c(1:27)])[2]),rep("negative",dim(negative_final[,-c(1:27)])[2])))
feature_meta<-cbind(rbind(cbind(mz_pos,rt_pos),cbind(mz_neg,rt_neg)),as.data.frame(as.factor(mode)))
feature_meta<-feature_meta[-index,]
rownames(feature_meta)<-rownames(T_abunds)
colnames(feature_meta)<-c("MZ","RT","Mode")

####Features for Output 
feature_meta_use<-cbind(rownames(T_abunds),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

####Complete Data Object Construction 
feature_meta<-AnnotatedDataFrame(feature_meta)
combined<-ExpressionSet(assayData=edata,phenoData=pheno_meta,featureData=feature_meta)

########LEEK = 0 sig SVs
mod<-model.matrix(~as.factor(ClinStatus),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)
design_table<-cbind(model.matrix(~0+as.factor(ClinStatus),data=combined),as.data.frame(svobj$sv))
colnames(design_table)<-c("Control","Carrier","MCIAD","SV1","SV2","SV3",
							"SV4","SV5","SV6","SV7",
							"SV8","SV9","SV10","SV11",
							"SV12","SV13")

####Fit linear model for pairwise contrasts 
arrayw<-arrayWeights(edata, design=design_table)
fit1<-lmFit(edata,design_table,weights=arrayw)
cm1<- makeContrasts(
	`Carrier-MCIAD` = Carrier-MCIAD, 
	levels=design_table)

fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Carrier.MCIAD"),
			sep="\t",file="Ringman_Carrier_MCIAD_4_1.POS.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Carrier.MCIAD"),
			sep="\t",file="Ringman_Carrier_MCIAD_4_1.NEG.txt",row.names=FALSE)


#######PIUmet

write.table(joinT%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="PIUmet_Ringman_Carrier_MCIAD_9_2.txt",row.names=FALSE)





age_Rel_cor<- numeric(nrow(edata)) 
P.Value<- numeric(nrow(edata))
for(i in 1:nrow(edata)){
         age_Rel_cor[i]<-cor.test(edata[i,],pheno_plot$`Age Rel to PDx`,method="spearman",exact=FALSE)$estimate
         P.Value[i]<-cor.test(edata[i,],pheno_plot$`Age Rel to PDx`,method="spearman",exact=FALSE)$p.value
}
A<-cbind(feature_meta_use,age_Rel_cor,P.Value)
write.table(A%>%filter(Mode=="positive")%>%select("MZ","RT","P.Value","age_Rel_cor"),
			sep="\t",file="Ringman_CORR_9_18.POS_MUMMI.txt",row.names=FALSE)
write.table(A%>%filter(Mode=="negative")%>%select("MZ","RT","P.Value","age_Rel_cor"),
			sep="\t",file="Ringman_CORR_9_18.NEG_MUMMI.txt",row.names=FALSE)


write.table(A%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="Ringman_CORR_9_18_PIUmet.txt",row.names=FALSE)



age_Rel_cor<- numeric(nrow(edata)) 
P.Value<- numeric(nrow(edata)) 
for(i in 1:nrow(design_table)){
	     set.seed(122)
         age_Rel_cor[i]<- bcov.test(edata[i,],pheno_plot$`Age Rel to PDx`, method="limit")$statistic
         P.Value[i]<-bcov.test(edata[i,],pheno_plot$`Age Rel to PDx`,  method="limit")$p.value
}
B<-cbind(feature_meta_use,age_Rel_cor,P.Value)
B[B$P.Value<0.00000000]<-0.000000001
write.table(B%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="Ringman_CORR_9_18_PIUmet_bcov.txt",row.names=FALSE)






cd /Volumes/TJGross_Remote/318557/Active_Projects_4_9_2021/Ringman_4_2021

mummichog -f Ringman_CORR_9_18.POS_MUMMI.txt -o Ringman_CORR_9_18.POS_MUMMI -m positive
mummichog -f Ringman_CORR_9_18.NEG_MUMMI.txt -o Ringman_CORR_9_18.NEG_MUMMI -m negative  












fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)
T<-topTableF(fit2_F,adjust="BH",number=100000)
U<-cbind(rownames(T),T)
colnames(U)[1]<-"Feature"
joinU<-dplyr::inner_join(feature_meta_use,U)
write.table(joinU%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Carrier.Control"),
			sep="\t",file="Ringman_Carrier_Control_4_1.POS.txt",row.names=FALSE)
write.table(joinU%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Carrier.Control"),
			sep="\t",file="Ringman_Carrier_Control_4_1.NEG.txt",row.names=FALSE)

fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)
T<-topTableF(fit3_F,adjust="BH",number=100000)
V<-cbind(rownames(T),T)
colnames(V)[1]<-"Feature"
joinV<-dplyr::inner_join(feature_meta_use,V)
write.table(joinV%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCI.Carrier"),
			sep="\t",file="Ringman_MCI_Carrier_4_1.POS.txt",row.names=FALSE)
write.table(joinV%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCI.Carrier"),
			sep="\t",file="Ringman_MCI_Carrier_4_1.NEG.txt",row.names=FALSE)


###############
cd /Volumes/TJGross_Remote/318557/Active_Projects_4_9_2021/Ringman_4_2021

mummichog -f Ringman_Carrier_MCIAD_4_1.NEG.txt -o Ringman_Carrier_MCIAD_4_1.NEG -m negative
mummichog -f Ringman_Carrier_MCIAD_4_1.POS.txt -o Ringman_Carrier_MCIAD_4_1.POS -m positive  

mummichog -f Ringman_Carrier_Control_4_1.NEG.txt -o Ringman_Carrier_Control_4_1.NEG -m negative
mummichog -f Ringman_Carrier_Control_4_1.POS.txt -o Ringman_Carrier_Control_4_1.POS -m positive 

mummichog -f Ringman_AD_MCI_4_1.NEG.txt -o Ringman_AD_MCI_4_1.NEG -m negative
mummichog -f Ringman_AD_MCI_4_1.POS.txt -o Ringman_AD_MCI_4_1.POS -m positive 


##############
sex_SVs_nominalP<-apply(design_table[,5:17],2,function(x) wilcox.test(x~as.factor(pheno_plot$Sex), exact=FALSE)$p.val)      
apoe_SVs_nominalP<-apply(design_table[,5:17],2,function(x) kruskal.test(x, g=as.factor(pheno_plot$APOE.x), exact=FALSE)$p.val)  
gene_SVs_nominalP<-apply(design_table[,5:17],2,function(x) wilcox.test(x~as.factor(pheno_plot$Gene.x), exact=FALSE)$p.val)  
bound<-cbind(colnames(design_table[,5:17]),sex_SVs_nominalP,apoe_SVs_nominalP,gene_SVs_nominalP)
colnames(bound)[1]<-c("Surrogate_Variables")



age_Rel_cor<- numeric(ncol(design_table[,4:16])) 
age_Rel_sig<- numeric(ncol(design_table[,4:16]))
for(i in 4:ncol(design_table)){
         age_Rel_cor[i]<- cor.test(design_table[,i],pheno_plot$`Age Rel to PDx`,method="spearman",exact=FALSE)$estimate
         age_Rel_sig[i]<-cor.test(design_table[,i],pheno_plot$`Age Rel to PDx`,method="spearman",exact=FALSE)$p.value
}
cbind(age_Rel_cor,age_Rel_sig)

age_cor<- numeric(ncol(design_table[,5:17])) 
age_sig<- numeric(ncol(design_table[,5:17]))
for(i in 5:ncol(design_table)){
         age_cor[i]<- cor.test(design_table[,i],pheno_plot$Age,method="spearman",exact=FALSE)$estimate
         age_sig[i]<-cor.test(design_table[,i],pheno_plot$Age,method="spearman",exact=FALSE)$p.value
}
cbind(age_cor,age_sig)

MMSE_cor<- numeric(ncol(design_table[,5:17])) 
MMSE_sig<- numeric(ncol(design_table[,5:17]))
for(i in 5:ncol(design_table)){
         MMSE_cor[i]<- cor.test(design_table[,i],pheno_plot$MMSE,method="spearman",exact=FALSE)$estimate
         MMSE_sig[i]<-cor.test(design_table[,i],pheno_plot$MMSE,method="spearman",exact=FALSE)$p.value
}
cbind(MMSE_cor,MMSE_sig)

bound<-cbind(bound,age_Rel_sig[-c(1:4)],age_sig[-c(1:4)],MMSE_sig[-c(1:4)])
colnames(bound)[5:7]<-c("age_Rel_nominalP","age_nominalP","mmse_nominalP")
write.csv(bound,file="Conventional_NonParametric_RingmanSVs.csv")


















sex_SVs_nominalP<-apply(design_table[,5:17],2,function(x) bd.test(x~as.factor(pheno_plot$Sex), exact=FALSE)$p.val)      
apoe_SVs_nominalP<-apply(design_table[,5:17],2,function(x) bd.test(x~as.factor(pheno_plot$APOE.x), exact=FALSE)$p.val)  
gene_SVs_nominalP<-apply(design_table[,5:17],2,function(x) bd.test(x~as.factor(pheno_plot$Gene.x), exact=FALSE)$p.val)  
bound_ball<-cbind(colnames(design_table[,5:17]),sex_SVs_nominalP,apoe_SVs_nominalP,gene_SVs_nominalP)
colnames(bound_ball)[1]<-c("Surrogate_Variables")



age_Rel_cor<- numeric(ncol(design_table[,5:17])) 
age_Rel_sig<- numeric(ncol(design_table[,5:17]))
for(i in 5:ncol(design_table)){
	     set.seed(122)
         age_Rel_cor[i]<- bcov.test(design_table[,i],pheno_plot$`Age Rel to PDx`)$statistic
         age_Rel_sig[i]<-bcov.test(design_table[,i],pheno_plot$`Age Rel to PDx`)$p.value
}
cbind(age_Rel_cor,age_Rel_sig)

age_cor<- numeric(ncol(design_table[,5:17])) 
age_sig<- numeric(ncol(design_table[,5:17]))
for(i in 5:ncol(design_table)){
	     set.seed(122)
         age_cor[i]<- bcov.test(design_table[,i],pheno_plot$Age)$statistic
         age_sig[i]<-bcov.test(design_table[,i],pheno_plot$Age)$p.value
}
cbind(age_cor,age_sig)

MMSE_cor<- numeric(ncol(design_table[,5:17])) 
MMSE_sig<- numeric(ncol(design_table[,5:17]))
for(i in 5:ncol(design_table)){
	     set.seed(122)
         MMSE_cor[i]<- bcov.test(design_table[,i],pheno_plot$MMSE)$statistic
         MMSE_sig[i]<-bcov.test(design_table[,i],pheno_plot$MMSE)$p.value
}
cbind(MMSE_cor,MMSE_sig)

bound_ball<-cbind(bound_ball,age_Rel_sig[-c(1:4)],age_sig[-c(1:4)],MMSE_sig[-c(1:4)])
colnames(bound_ball)[5:7]<-c("age_Rel_nominalP","age_nominalP","mmse_nominalP")
write.csv(bound_ball,file="Ball_RingmanSVs.csv")



























################
all_together<-cbind(design_table,pheno_plot)
pheno_factors<-apply(all_together[,c(19,21,24,28)],2,as.factor)
pheno_cont<-apply(all_together[,c(4:16,22:23,25:26)],2,as.numeric)
discrete<-cbind(as.data.frame(discretize.jointly(pheno_cont,k=c(2:50))$D),as.data.frame(pheno_factors))
colnames(discrete)[c(19,21)]<-c("APOE","Mutation")
discrete[,1:dim(discrete)[2]]<-lapply(discrete[,1:dim(discrete)[2]],as.factor)
discrete<-discrete[,-caret::nearZeroVar(discrete)]

########SV8-SEX
G<-h2pc(discrete)
H<-arc.strength(G,discrete)
strength.plot(G,H)