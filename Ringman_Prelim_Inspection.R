################## Are features differentially distributed according to clinical group in Ringman? PCA with outlier sample[11]

expression<-t(exprs(combined))[-c(11),]
pheno<-combined$ClinStatus[-c(11)]

String_P<- numeric(ncol(expression)) 
for(i in 1:ncol(expression)){
               String_P[i] <- bd.test(expression[,i]~pheno)$p.value
}

String_S<-numeric(ncol(expression)) 
for(i in 1:ncol(expression)){
               String_S[i] <- bd.test(expression[,i]~pheno)$statistic
}


187 features nominal p < .05; plot them 

A<-as.numeric(colnames(expression[,String<.05]))

pdf(file="Ring.pdf")
for(i in 1:length(A)){
 	vioplot(expression[,A[i]]~pheno,main=A[i])
}
dev.off()


################## Are clusters of features formed that are differentially distributed in Ringman? 

datExpr<-t(exprs(combined))[-c(11),]
####ALL OK   gsg<-goodSamplesGenes(datExpr, verbose = 3)

### Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

###Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf(file="choice_Ringman.pdf")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###Define soft thresholding power to be used
softPower <- 9
adjacency <- adjacency(datExpr, power = softPower,type="signed")

###Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency);
dissTOM <- 1-TOM

### Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

###Define module size
minModuleSize <- 50

### Identify modules using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

### Convert numeric labels into colors
set.seed(122)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

###Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
colors<-MEList$validColors

###
###

String<- numeric(ncol(MEs)) 
for(i in 1:ncol(MEs)){
               String[i] <- bd.test(MEs[,i]~pheno)$p.value
}

##################two eigenfeatures trending at nominal p=.07; RED and SALMON  

mode<-as.factor(c(rep("ESI+",dim(positive_final[,-c(1:11)])[2]),rep("ESI-",dim(negative_final[,-c(1:11)])[2])))
feature_meta<-cbind(rbind(cbind(mz_pos,rt_pos),cbind(mz_neg,rt_neg)),as.data.frame(as.factor(mode)))
feature_meta<-feature_meta[-index,]
colnames(feature_meta)<-c("MZ","RT","Mode")
rownames(feature_meta)<-seq(1,dim(abunds)[2])

feature_meta_use<-cbind(as.data.frame(seq(1,dim(abunds)[2])),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

write.csv(feature_meta[colors=="red",]%>%filter(Mode=="ESI+"),file="red_pos.csv")
write.csv(feature_meta[colors=="red",]%>%filter(Mode=="ESI-"),file="red_neg.csv")

write.csv(feature_meta[colors=="salmon",]%>%filter(Mode=="ESI+"),file="salmon_pos.csv")
write.csv(feature_meta[colors=="salmon",]%>%filter(Mode=="ESI-"),file="salmon_neg.csv")



G<-bcorsis(x=expression,y=as.numeric(pheno_meta$`Age Rel to PDx`)[-c(11)])
H<-G$ix

pdf(file="Ring_corsis.pdf")
for(i in 1:length(H)){
 	plot(as.numeric(pheno_meta$`Age Rel to PDx`)[-c(11)],expression[,H[i]],,main=H[i])
}
dev.off()


F<-bcorsis(x=MEs,y=as.numeric(pheno_meta$`Age Rel to PDx`)[-c(11)])
H<-F$ix

pdf(file="Ring_corsis_WGCNA.pdf")
for(i in 1:length(H)){
 	plot(as.numeric(pheno_meta$`Age Rel to PDx`)[-c(11)],MEs[,i],,main=colnames(MEs[i]))
}
dev.off()



expression<-t(exprs(combined))[-c(11),]
pheno<-combined$ClinStatus[-c(11)]

String<- numeric(ncol(expression)) 
for(i in 1:ncol(expression)){
               String[i] <- bd.test(expression[,i]~pheno)$p.value
}





