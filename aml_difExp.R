# Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia, N Engl J Med. 2013 May 30; 368(22): 2059–2074. 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("Biobase", "limma"))

library(Biobase)
library(limma)

#####################################################
# Clean data 
#####################################################
# get feature data
x<-read.csv(file="https://raw.githubusercontent.com/YunMai-SPS/RNA_difExp/master/laml_tcga_pub/data_RNA_Seq_v2_expression_median.txt",header=T,sep="\t")
x <-x[!duplicated(x$Entrez_Gene_Id),] 

f<-x[,1:2]
rownames(f)<-f[,2]

p<-read.delim(file="https://raw.githubusercontent.com/YunMai-SPS/RNA_difExp/master/laml_tcga_pub/data_clinical_patient.txt",sep="\t",header=T,stringsAsFactors=F)
p<- p[5:nrow(p),]
p[,c(5,6,7,15,20)] <- apply(p[,c(5,6,7,15,20)], 2, function(x) as.numeric(x))
colnames(p)[colnames(p)=='X.Patient.Identifier'] <- "patientID"
p <-p[!duplicated(p$patientID),] 

rownames(x)<-x[,2]
x<-x[,-c(1,2)]
colnames(x) <- gsub("\\.03","",colnames(x))
colnames(x) <- gsub("\\.","-",colnames(x))
x<-as.matrix(x)

p<-p[p$patientID %in% colnames(x),]
x<-x[,c(p$patientID)]

# make sure the samples in assaydata x and in patient info are in same order
p$patientID == colnames(x)
rownames(p)<-p[,1]
p<-p[,-1]

# make sure the order of rownames of feature set and assaydata x are in the same order
table(rownames(f)==rownames(x))

write.csv(f,"f_AML.csv")
write.csv(x,"x_AML.csv")
write.csv(p,"p_AML.csv")

# visualize gene expression between different sample with survival status 
boxplot(x[1000,] ~ p[,'Overall.Survival.Status'],main=f[1000,"Hugo_Symbol"])

############################################
# creat ExpressionSet object
############################################
eset <- ExpressionSet(assayData = x,
phenoData = AnnotatedDataFrame(p),
featureData = AnnotatedDataFrame(f))

######################################################
# preprocess data
######################################################
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, legend = FALSE)

############################################
# linear model
############################################
design <- model.matrix(~Overall.Survival.Status, data = pData(eset))
table(pData(eset)[,'Overall.Survival.Status'])
fit<-lmFit(eset,design)
fit<-eBayes(fit)
head(fit$coffecients,3)
results<-decideTests(fit[,'Overall.Survival.StatusLIVING'])
summary(results)

############################################
# flexible linear model
############################################
design<-model.matrix(~0+Overall.Survival.Status, data = pData(eset))
colSums(design)
cm<-makeContrasts(status=Overall.Survival.StatusLIVING-Overall.Survival.StatusDECEASED,
levels=design)
fit<-lmFit(eset,design)
head(fit$coefficients, 3)
fit2<-contrasts.fit(fit,contrasts=cm)
head(fit2$coffecients,3)
fit2<-eBayes(fit2)
results<-decideTests(fit2)
summary(results)

################################################################
# visualize the data

topTable(fit2, number = 3)
stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
hist(stats[, "P.Value"])
volcanoplot(fit2, highlight = 5, names = fit2$genes[, "Hugo_Symbol"])

# enrichment testing, KEGG
entrez <- fit2$genes[, "Entrez_Gene_Id"]
enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg, number = 3)

# enrichment testing, GO
enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP", number = 3)