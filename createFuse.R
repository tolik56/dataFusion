library(ROCit)
library(ROCR)
library(plotROC)
library(caret)
library(glmnet)
library(dplyr)
library(readxl)
library(rlist)
# LASSO adds intersect to the list of coefficients
# before using the fist element should be excluded and col number shifted by one
#
setwd("~/Proj/model")

dvMut <- read_excel("Validation43.xlsx", sheet="Mut")
dvMet <- read_excel("Validation43.xlsx", sheet="Met")
dvCNV <- read_excel("Validation43.xlsx", sheet="CNV")
dvEnd <- read_excel("Validation43.xlsx", sheet="End")
dvPro <- read_excel("Validation43.xlsx", sheet="Pro")
dv <- Reduce(function(x, y) merge(x, y, by = "Sample"), list(dvMut, dvPro, dvMet, dvCNV, dvEnd))

daMut <- read_excel("DataRaw175.xlsx", sheet="Mut"); namesMut <- names(daMut)[4]
daMet <- read_excel("DataRaw175.xlsx", sheet="Met"); namesMet <- names(daMet)[-c(1:3)]
daCNV <- read_excel("DataRaw175.xlsx", sheet="CNV"); namesCNV <- names(daCNV)[-c(1:3)]
daEnd <- read_excel("DataRaw175.xlsx", sheet="End"); namesEnd <- names(daEnd)[-c(1:3)]
daPro <- read_excel("DataRaw175.xlsx", sheet="Pro"); namesPro <- names(daPro)[-1]
da <- Reduce(function(x, y) merge(x, y, by = "Sample"), list(daMut[c(1:4)], daPro, daMet[-c(2:3)], daCNV[-c(2:3)], daEnd[-c(2:3)]))

da$Condition <- as.factor(da$Condition)
fuse <- da[1]
probma <- dv[1]

bestSeeds <-c(5,186,246,299)
seed <- 5
CUTOFF <- 0.47

for(i in seq(3)) {
  set.seed(seed)
  if(i==1) {labS = "Fixed_"; index <- grep("train",da$Set)}
  if(i==2) {labS = "Random_"; index <- index <- createDataPartition(da$Condition, p=0.5, list=FALSE)}
  if(i==3) {labS = "Total_"; index <- seq(nrow(da))}
  modelSet <- list()
  
  #### Methylation #########################################################################################
  #
  namesMet <- c("B3GALT4", "RUNX2", "EMX1", "MKL1", "PSD", "TSPYL5", "CFTR", "TCF24", "BCL2L11", "DCDC2", "ZFHX3", "OPLAH", "PTPN7", 
                "LINC00538PROX1", "TIAM1", "NFIX", "FAM78A", "FAM55C", "DNAH10")
  dd <- da[namesMet]
  
  set.seed(seed)
  model <- train(x = dd[index,], y = da$Condition[index], 
                 method="glmnet",
                 verbose= FALSE,
                 trControl = trainControl(method = "none",classProbs = TRUE),
                 tuneGrid = expand.grid(alpha= 0.8,lambda=0.006615977),
                 metric= "ROC")
  
  cutoff <- 0.47
  pred <- predict(model, newdata = dd, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  confusionMatrix(predClass, da$Condition, positive = "HCC")
  (AUC <- round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))
  
  fuse$Met <- pred$HCC
  modelSet <- list.append(modelSet, Met=list(fit=model, varli=names(dd)))
  
  #### End Motif #########################################################################################
  #
  namesEnd <- c("ACGA", "CACT", "GCGT", "CATT", "GAGG", "AACT", "ACGT", "CGAT", "CGAC", "ACGG", 
                "AATC", "GTGG", "CGAG", "CGAA", "TCCG", "GTGC", "GCGC", "TACG", "TACC", "GCGA")
  dd <- da[namesEnd]
  
  set.seed(seed)
  model <- train(x = dd[index,], y = da$Condition[index], 
                 method="glmnet",
                 verbose= FALSE,
                 trControl = trainControl(method = "none",classProbs = TRUE),
                 tuneGrid = expand.grid(alpha= 0.1,lambda=0.1965686),
                 metric= "ROC")
  
  cutoff <- 0.5
  pred <- predict(model, newdata = dd, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  confusionMatrix(predClass, da$Condition, positive = "HCC")
  (AUC <- round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))
  
  fuse$End <- pred$HCC
  modelSet <- list.append(modelSet, End=list(fit=model, varli=names(dd)))
  
  #### CNV #########################################################################################
  #
  dd <- da[namesCNV]
  
  set.seed(seed)
  model <- train(x = dd[index,], y = da$Condition[index], 
                 method="glmnet",
                 verbose= FALSE,
                 trControl = trainControl(method = "none",classProbs = TRUE),
                 tuneGrid = expand.grid(alpha= 0.1,lambda=0.005767354),
                 metric= "ROC")

  cutoff <- 0.5
  pred <- predict(model, newdata = dd, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  confusionMatrix(predClass, da$Condition, positive = "HCC")
  (AUC <- round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))
  
  fuse$CNV <- pred$HCC
  modelSet <- list.append(modelSet, CNV=list(fit=model, varli=names(dd)))
  
  #### Validator ######################################################
  #
  for( j in seq(2)) {
    if(j==1){labL <- "MMCE_Fusion"; varli <- namesMut}
    if(j==2){labL <- "MMCEP_Fusion"; varli <- c(namesMut, namesPro)}
    
    tr <- cbind(da[varli], fuse[-1])
    
    model <- train(x= tr[index,], y= da$Condition[index], 
                   method="glmnet",
                   verbose= FALSE,
                   trControl = trainControl(method = "none",classProbs = TRUE),
                   tuneGrid = expand.grid(alpha= 0.7,lambda=0.01139748),
                   metric= "ROC")
    
    dd <- cbind(dv[varli], predict(modelSet$Met$fit, newdata = dv[,modelSet$Met$varli], 'prob')[2],
                predict(modelSet$End$fit, newdata = dv[,modelSet$End$varli], 'prob')[2],
                predict(modelSet$CNV$fit, newdata = dv[,modelSet$CNV$varli], 'prob')[2])
    names(dd) <- names(tr)
    cutoff <- CUTOFF
    pred <- predict(model, newdata = dd, 'prob')[2]
    predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
    probma[,paste0(labS,labL)] <- pred
    
    # Prediction for original set of 175 samples
    # AUC <- unname(round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))
    # fpfn <- confusionMatrix(predClass, da$Condition, positive = "HCC")$table[c(2,3)]
    # print(paste(cutoff, labS, labL, fpfn[1], fpfn[2], AUC, sep=","))
    
}}
clipr::write_clip(probma)
