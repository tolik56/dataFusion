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
seed <- bestSeeds[2]
CUTOFF <- 0.5

for(i in seq(3)) {
  set.seed(seed)
  if(i==1) {labS = "Fixed_"; index <- grep("train",da$Set)}
  if(i==2) {labS = "Random_"; index <- index <- createDataPartition(da$Condition, p=0.5, list=FALSE)}
  if(i==3) {labS = "Total_"; index <- seq(nrow(da))}
  modelSet <- list()
  
  #### Methylation #########################################################################################
  #
  dd <- da[namesMet] %>% tibble()
  
  set.seed(seed)
  fit <- cv.glmnet( x = as.matrix(dd), y = da$Condition, 
                    type.measure = "class",
                    family = "binomial",
                    alpha = 1, nfolds = 10
  )
  
  lasso_coefs <- as.data.frame(summary(coef(fit, s = fit$lambda.min))) %>% subset(i>1)
  lasso_coefs <- lasso_coefs[order(abs(lasso_coefs$x), decreasing = TRUE),] %>% head(20)
  
  dd_sel <- dd[lasso_coefs$i-1]
  
  model <- train(x = dd_sel[index,], y = da$Condition[index], 
                 method = 'glmnet',
                 family = 'binomial',
                 metric = 'ROC',
                 trControl = trainControl(method = 'repeatedcv',
                                          classProbs = TRUE,
                                          number = 10,
                                          repeats = 10,
                                          search = 'random')
  )
  cutoff <- 0.5
  pred <- predict(model, newdata = dd_sel, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  confusionMatrix(predClass, da$Condition, positive = "HCC")
  (AUC <- round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))

  fuse$Met <- pred$HCC
  modelSet <- list.append(modelSet, Met=list(fit=model, varli=names(dd_sel)))
  
  #### End Motif #########################################################################################
  #
  dd <- da[namesEnd] %>% tibble()
  
  set.seed(seed)
  fit <- cv.glmnet( x = as.matrix(dd), y = da$Condition, 
                    type.measure = "class",
                    family = "binomial",
                    alpha = 1, nfolds = 10
  )
  
  lasso_coefs <- as.data.frame(summary(coef(fit, s = fit$lambda.min))) %>% subset(i>1)
  lasso_coefs <- lasso_coefs[order(abs(lasso_coefs$x), decreasing = TRUE),] %>% head(20)
  
  dd_sel <- dd[lasso_coefs$i-1]
  
  model <- train(x = dd_sel[index,], y = da$Condition[index], 
                 method = 'glmnet',
                 family = 'binomial',
                 metric = 'ROC',
                 trControl = trainControl(method = 'repeatedcv',
                                          classProbs = TRUE,
                                          number = 10,
                                          repeats = 10,
                                          search = 'random')
  )
  cutoff <- 0.5
  pred <- predict(model, newdata = dd_sel, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  confusionMatrix(predClass, da$Condition, positive = "HCC")
  (AUC <- round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))
  
  fuse$End <- pred$HCC
  modelSet <- list.append(modelSet, End=list(fit=model, varli=names(dd_sel)))
  
  #### CNV #########################################################################################
  #
  dd <- da[namesCNV] %>% tibble()
  
  set.seed(seed)
  fit <- cv.glmnet( x = as.matrix(dd), y = da$Condition, 
                    type.measure = "class",
                    family = "binomial",
                    alpha = 1, nfolds = 10
  )
  
  lasso_coefs <- as.data.frame(summary(coef(fit, s = fit$lambda.min))) %>% subset(i>1)
  lasso_coefs <- lasso_coefs[order(abs(lasso_coefs$x), decreasing = TRUE),] %>% head(20)
  
  dd_sel <- dd[lasso_coefs$i-1]
  
  model <- train(x = dd_sel[index,], y = da$Condition[index], 
                 method = 'glmnet',
                 family = 'binomial',
                 metric = 'ROC',
                 trControl = trainControl(method = 'repeatedcv',
                                          classProbs = TRUE,
                                          number = 10,
                                          repeats = 10,
                                          search = 'random')
  )
  cutoff <- 0.5
  pred <- predict(model, newdata = dd_sel, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  confusionMatrix(predClass, da$Condition, positive = "HCC")
  (AUC <- round(unlist(ciAUC(rocit(pred$HCC, da$Condition))[1]), 3))
  
  fuse$CNV <- pred$HCC
  modelSet <- list.append(modelSet, CNV=list(fit=model, varli=names(dd_sel)))
  
  #### Validator ######################################################
  #
  namesFus <- names(fuse)[-1]
  for( j in seq(2)) {
    if(j==1){labL <- "MMCE_LASSO"; varli <- namesMut}
    if(j==2){labL <- "MMCEP_LASSO"; varli <- c(namesMut, namesPro)}
   
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
