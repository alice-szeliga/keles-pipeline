setwd("/p/keles/CAGI2015/volumeA/Data/After_Nov5/")
trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
dim(trainingD)
#[1] 3044   25
length(unique(trainingD$eQTL_gene))
#[1] 1052



scoreList <- vector("list", 13)
names(scoreList) <- c("fb_minPval", "fb_maxScore", "fb_numOcc", "fb_sumScore",
                      "en_minPval", "en_maxScore", "en_numOcc", "en_sumScore",
                      "ja_minPval", "ja_maxScore", "ja_numOcc", "ja_sumScore", "db_score")




# Code for exploratory analysis and extracting scores
#########################################
#                                       #
# Exploratory analysis with FactorBook  #
# PWM scores                            #
#                                       #
#########################################
setwd("/p/keles/CAGI2015/volumeB/Data/")
trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)


#setwd("/Users/keles/Research_Projects/CAGI2015/Fimo/")
setwd("/p/keles/CAGI2015/volumeB/Fimo")
load("Training_ref_thres0.001_FactorBook_minPval.RData")
load("Training_ref_thres0.001_FactorBook_maxScore.RData")
load("Training_ref_thres0.001_FactorBook_sumScore.RData")
load("Training_ref_thres0.001_FactorBook_numOcc.RData")


fimo.minPval <- out.1
fimo.maxScore <- out.2
fimo.numOcc <- out.3
fimo.sumScore <- out.4


fimoS <- fimo.minPval
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(1, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["fb_minPval"]] <- fimoS.ordered


fimoS <- fimo.maxScore
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["fb_maxScore"]] <- fimoS.ordered


fimoS <- fimo.numOcc
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["fb_numOcc"]] <- fimoS.ordered

fimoS <- fimo.sumScore
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["fb_sumScore"]] <- fimoS.ordered


# Simple correlations
fimoS.ordered <- scoreList[["fb_minPval"]]
fimoS.ordered.binary <- apply(fimoS.ordered, 2, function(x){res <- x; res[x>=0.05]<-0; res[x <0.05] <-1; return(res)}) #PWM match crude
corRes.binary <- apply(fimoS.ordered.binary, 2, function(x){ifelse(length(unique(x)) > 1, wilcox.test(trainingD$C.A.log2FC~x)$p.value,1)})
#order(corRes.binary, decreasing = FALSE)
cat("minPVal\n")
summary(-log10(corRes.binary))

cat("maxScore\n")
fimoS.ordered <- scoreList[["fb_maxScore"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
summary(abs(corRes))

cat("sumScore\n")
fimoS.ordered <- scoreList[["fb_sumScore"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
summary(abs(corRes))

cat("numOcc\n")
fimoS.ordered <- scoreList[["fb_numOcc"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
summary(abs(corRes))


#####################################
#                                   #
# Exploratory analysis with Deepbind#
# scores                            #
#                                   #
#####################################
setwd("/p/keles/CAGI2015/volumeA/Data/After_Nov5")
trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)

setwd("/p/keles/CAGI2015/volumeB/DeepBind/deepbind")
db <- read.table("ref_SS.v1_deepbind.txt", header = TRUE)

scoreList[["db_score"]] <- db

## Plots
setwd("/Users/keles/Research_Projects/CAGI2015/Data/After_Nov5")
trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)
setwd("/Users/keles/Research_Projects/CAGI2015/DeepBind/")
db <- read.table("Training_ref_3044_deepbind_927.txt", header = TRUE)
corRes <- apply(db, 2, function(x){cor(x, trainingD$C.A.log2FC)})
index <- order(abs(corRes), decreasing = TRUE)[1:50]

require(ggplot2)
require(reshape2)

db.sub <- db[, index]
colnames(db.sub) <- paste("db_", 1:50, sep = "")

dd <- cbind(trainingD, db.sub)
ddd <- dd
ddd$grp <- paste(ddd$Regulatory_Hit, ddd$emVar_Hit, sep = "")
ggplot(ddd, aes(x = db_2, y =C.A.log2FC)) + geom_smooth(method='lm',formula=y~x) + geom_point(aes(colour=grp)) + facet_wrap( ~ grp, nrow = 1)

#########################################
#                                       #
# Exploratory analysis with ENCODE      #
# PWM scores                            #
#                                       #
#########################################
setwd("/p/keles/CAGI2015/volumeA/Data/After_Nov5")
trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)


#setwd("/Users/keles/Research_Projects/CAGI2015/Fimo/")
setwd("/p/keles/CAGI2015/volumeB/Fimo")
load("Training_ref_thres0.001_ENCODE_minPval.RData")
load("Training_ref_thres0.001_ENCODE_maxScore.RData")
load("Training_ref_thres0.001_ENCODE_sumScore.RData")
load("Training_ref_thres0.001_ENCODE_numOcc.RData")


fimo.minPval <- out.1
fimo.maxScore <- out.2
fimo.numOcc <- out.3
fimo.sumScore <- out.4


fimoS <- fimo.minPval
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(1, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["en_minPval"]] <- fimoS.ordered


fimoS <- fimo.maxScore
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["en_maxScore"]] <- fimoS.ordered


fimoS <- fimo.numOcc
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["en_numOcc"]] <- fimoS.ordered

fimoS <- fimo.sumScore
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["en_sumScore"]] <- fimoS.ordered


# Simple correlations
fimoS.ordered <- scoreList[["en_minPval"]]
fimoS.ordered.binary <- apply(fimoS.ordered, 2, function(x){res <- x; res[x>=0.05]<-0; res[x <0.05] <-1; return(res)}) #PWM match crude
corRes.binary <- apply(fimoS.ordered.binary, 2, function(x){ifelse(length(unique(x)) > 1, wilcox.test(trainingD$C.A.log2FC~x)$p.value,1)})
#order(corRes.binary, decreasing = FALSE)
cat("minPVal\n")
summary(-log10(corRes.binary))


fimoS.ordered <- scoreList[["en_maxScore"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
cat("maxScore\n")
summary(abs(corRes))

fimoS.ordered <- scoreList[["en_sumScore"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
cat("sumScore\n")
summary(abs(corRes))


fimoS.ordered <- scoreList[["en_numOcc"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
cat("numOcc\n")
summary(abs(corRes))


#########################################
#                                       #
# Exploratory analysis with Jaspar      #
# PWM scores                            #
#                                       #
#########################################
setwd("/p/keles/CAGI2015/volumeB/Data/")
trainingD <- read.table("4-eQTL-causal_SNPs_sample_v2.txt", header = T)


#setwd("/Users/keles/Research_Projects/CAGI2015/Fimo/")
setwd("/p/keles/CAGI2015/volumeB/Fimo")
load("Training_ref_thres0.001_Jaspar_minPval.RData")
load("Training_ref_thres0.001_Jaspar_maxScore.RData")
load("Training_ref_thres0.001_Jaspar_sumScore.RData")
load("Training_ref_thres0.001_Jaspar_numOcc.RData")


fimo.minPval <- out.1
fimo.maxScore <- out.2
fimo.numOcc <- out.3
fimo.sumScore <- out.4


fimoS <- fimo.minPval
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(1, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["ja_minPval"]] <- fimoS.ordered


fimoS <- fimo.maxScore
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["ja_maxScore"]] <- fimoS.ordered


fimoS <- fimo.numOcc
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["ja_numOcc"]] <- fimoS.ordered

fimoS <- fimo.sumScore
oIndex <- as.numeric(apply(as.matrix(rownames(fimoS)), 1, function(x){strsplit(x, "_SNP_")[[1]][2]}))
fimoS.ordered <- matrix(0, nrow(trainingD), ncol(fimoS))
colnames(fimoS.ordered) <- colnames(fimoS)
fimoS.ordered[oIndex, ] <- fimoS
scoreList[["ja_sumScore"]] <- fimoS.ordered


# Simple correlations
fimoS.ordered <- scoreList[["ja_minPval"]]
fimoS.ordered.binary <- apply(fimoS.ordered, 2, function(x){res <- x; res[x>=0.05]<-0; res[x <0.05] <-1; return(res)}) #PWM match crude
corRes.binary <- apply(fimoS.ordered.binary, 2, function(x){ifelse(length(unique(x)) > 1, wilcox.test(trainingD$C.A.log2FC~x)$p.value,1)})
#order(corRes.binary, decreasing = FALSE)
cat("minPVal\n")
summary(-log10(corRes.binary))

cat("maxScore\n")
fimoS.ordered <- scoreList[["ja_maxScore"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
summary(abs(corRes))

cat("sumScore\n")
fimoS.ordered <- scoreList[["ja_sumScore"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
summary(abs(corRes))

cat("numOcc\n")
fimoS.ordered <- scoreList[["ja_numOcc"]]
corRes <- apply(fimoS.ordered, 2, function(x){cor(x, trainingD$C.A.log2FC)})
#order(abs(corRes), decreasing = TRUE)
#sort(abs(corRes), decreasing = TRUE)[1:10]
summary(abs(corRes))




#########################################
#                                       #
# Output scores                         #
#                                       #
#                                       #
#########################################
unlist(lapply(scoreList, ncol))
fb_minPval fb_maxScore   fb_numOcc fb_sumScore  en_minPval en_maxScore
79          79          79          79        2065        2065
en_numOcc en_sumScore  ja_minPval ja_maxScore   ja_numOcc ja_sumScore
2065        2065        1082        1082        1082        1082
db_score
927

setwd("/p/keles/CAGI2015/volumeB/Predictors/")
save(scoreList, file = "scoreList_Training_ref_thres0.001_pwm.RData")








