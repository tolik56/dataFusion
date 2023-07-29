library(ROCit)
library(ROCR)
library(plotROC)
library(kernlab)
library(caret)
library(dplyr)
library(readxl)
library(ggfortify)
library(randomForest)

# clipr::write_clip(cbind(pred, predClass))

setwd("~/Proj/model")
df <- read_excel("DataBook175.xlsx", sheet="all")
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],  as.factor)


##### Main model ############################################################
#
index <- grep("train",df$Set)
varli <- c(3:7,28:30)
tr <- df[index, varli]
set.seed(111)
model <- train(x= tr[-1] , y = tr$Condition, 
               method="glmnet",
               verbose= FALSE,
               trControl = trainControl(method = "none",classProbs = TRUE),
               tuneGrid = expand.grid(alpha= 0.7,lambda=0.01139748),
               metric= "ROC" )

cutoff <- 0.47
pred <- predict(model, newdata = df[varli[-1]], 'prob')[2]
predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(tr)[-1]); confusionMatrix(predClass, df$Condition, positive = "HCC")
round(unlist(ciAUC(rocit(pred$HCC, df$Condition))[c(1,5,6)]), 3)

#### GALAD ##################################################################
#
predClass <- factor(ifelse(df$GALAD > -0.62, "HCC", "control"))
confusionMatrix(predClass, df$Condition, positive = "HCC")
round(unlist(ciAUC(rocit(df$GALAD, df$Condition))[c(1,5,6)]), 3)


da <- subset(df, Stage == "0" )
da <- subset(df, Stage == "A" )
da <- subset(df, Stage == "B" )
da <- subset(df, Stage == "C" )
da <- subset(df, Group == "ASH")
da <- subset(df, Group == "NASH")
da <- subset(df, Group == "VIR")

pred <- predict(model, newdata = da[varli], 'prob')[2]
predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
confusionMatrix(predClass, da$Condition, positive = "HCC")

clipr::write_clip(cbind(pred, predClass))
#
#### Loop to get individual models AUCs ###########################################
#
index <- grep("train",df$Set)
tr <- df[index,]
te <- df[-index,]
for(i in 4:7)  {
  a <- as.numeric(round(unlist(ciAUC(rocit(unlist(tr[i]), tr$Condition))[c(1,5,6)]), 3))
  print(paste0("tr ",a[1]," (95% CI:",a[2],"-",a[3],")"))
  a <- as.numeric(round(unlist(ciAUC(rocit(unlist(te[i]), te$Condition))[c(1,5,6)]), 3))
  print(paste0("te ",a[1]," (95% CI:",a[2],"-",a[3],")"))
}
li <- list(c(3,28:30), c(3:7,28:30))
for(i in 1:2) {
  varli <- li[[i]]
  tr <- df[index, varli]
  te <- df[-index, varli]
  model <- suppressWarnings(train(x= tr[-1] , y = tr$Condition, 
                                  method="glmnet",
                                  verbose= FALSE,
                                  trControl = trainControl(method = "none",classProbs = TRUE),
                                  tuneGrid = expand.grid(alpha= 0.7,lambda=0.01139748),
                                  metric= "ROC" ))
  
  pred <- predict(model, newdata = tr, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  a <- as.numeric(round(unlist(ciAUC(rocit(pred$HCC, tr$Condition))[c(1,5,6)]), 3))
  print(paste0("tr ",a[1]," (95% CI:",a[2],"-",a[3],")"))
  
  pred <- predict(model, newdata = te, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "conteol"))
  a <- as.numeric(round(unlist(ciAUC(rocit(pred$HCC, te$Condition))[c(1,5,6)]), 3))
  print(paste0("te ",a[1]," (95% CI:",a[2],"-",a[3],")"))
}

#
#################################################################
varli <- c(2:5)
varli <- c(2:5,7:9)
tr <- df[index, varli]
set.seed(111)
model <- train(x= tr[-1] , y = tr$Condition, 
               method="glmnet",
               verbose= FALSE,
               trControl = trainControl(method = "none",classProbs = TRUE),
               tuneGrid = expand.grid(alpha= 0.8,lambda=0.006615977),
               metric= "ROC" )

cutoff <- 0.4 #83/77 #56
pred <- predict(model, newdata = df[varli], 'prob')[2]
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(tr)[-1]); confusionMatrix(pred, df$Condition, positive = "HCC")

#
#################################################################
varli <- c(2:4)
varli <- c(2:4,7:9)
tr <- df[index, varli]
set.seed(111)
model <- train(x= tr[-1] , y = tr$Condition, 
               method="glmnet",
               verbose= FALSE,
               trControl = trainControl(method = "none",classProbs = TRUE),
               tuneGrid = expand.grid(alpha= 0.8,lambda=0.01139748),
               metric= "ROC" )
cutoff <- 0.6 #
pred <- predict(model, newdata = df[varli], 'prob')[2]
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(tr)[-1]); confusionMatrix(pred, df$Condition, positive = "HCC")


Sum7 <- readRDS("Sum7.rds")
df <- read_excel("DataBook175bis.xlsx", sheet="finalModel")
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],  as.factor)

d <- df[,2:4] # RACE
dd <- df[,c(2:4,16:18)]
dd <- df[,c(2:4,16:18,6,7)]

dd <- df[,c(2:4,8)] # RACE + all
dd <- df[,c(2:4,8,16:18)]
dd <- df[,c(2:4,8,16:18,6,7)]

dd <- df[,c(2:4,8)] # RACE + motif
dd <- df[,c(2:4,8,16:18)]
dd <- df[,c(2:4,8,16:18,6,7)]

dd <- df[,c(2:4,8,19)] # RACE + all + motif
dd <- df[,c(2:4,8,19,16:18)]
dd <- df[,c(2:4,8,19,16:18,6,7)]

dd <- cbind (df[,2:4], Sum7)
dd <- cbind (df[,c(2:4,5:7)], Sum7)
dd <- cbind (df[,c(2:4,16:18)], Sum7)
dd <- cbind (df[,c(2:4,16:18,5:7)], Sum7)

model <- glm(Condition ~ ., family = "binomial", data=dd)

cutoff <- 0.5
pred <- predict(model, newdata = dd, type = "response")
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(dd[-1])); confusionMatrix(pred, dd$Condition, positive="HCC")



set.seed(55)
index <- createDataPartition(df$Condition, p=0.5, list=FALSE)
tr <- dd[index,]

model <- train(x= tr[-1] , y = tr$Condition, 
               method = 'glmnet',
               family = 'binomial',
               trControl = trainControl(method = 'repeatedcv',
                                        number = 5,
                                        repeats =  5,
                                        search = 'random' )
)

cutoff <- 0.5
pred <- predict(model, newdata = dd, 'prob')[2]
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(dd)[-1]); confusionMatrix(pred, dd$Condition, positive = "HCC")






set.seed(22)
index <- createDataPartition(df$Condition, p=0.5, list=FALSE)

varli <- c(2,3,4,7,16:18)
tr <- df[index, varli]
te <- df[-index, varli]

model <- train(x= tr[-1] , y = tr$Condition, 
               method = 'glmnet',
               family = 'binomial',
               trControl = trainControl(method = 'repeatedcv',
                                        number = 5,
                                        repeats =  5,
                                        search = 'random' )
)

da <- subset(df, Stage == "0" | Stage == "A" | is.na(Stage))
da <- subset(df, Stage == "B" | Stage == "C" | is.na(Stage))
da <- subset(df, Group == "ASH")
da <- subset(df, Group == "NASH")
da <- subset(df, Group == "VIR")

cutoff <- 0.5
pred <- predict(model, da[varli][-1], 'prob')[2]
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(da[-1])); confusionMatrix(pred, da$Condition, positive="HCC")





# Use PC variables
library(ggbiplot)

pred <- prcomp(df[,-c(1:11)], center = TRUE, scale. = TRUE, rank.=7)
ggbiplot(pred, groups=df$Condition, var.axes=F)
dd <- cbind(df[,c(2:7)], predict(pred, df[,-c(1:11)]))

model <- glm(Condition ~ ., family = "binomial", data=dd)

cutoff <- 0.5
pred <- predict(model, newdata = dd, type = "response")
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(dd[-1])); confusionMatrix(pred, dd$Condition, positive="HCC")



library(rattle)
rattle()
######### DataBook 175 ##########################################################

setwd("~/Proj/model")

df <- read_excel("DataBook175.xlsx", sheet="final")
df$predClass[df$predClass=="nonHCC"] <- "control"
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],  as.factor)
df <- df %>% select_if(~ !any(is.na(.)))

da <- subset(df, subset != "out")
# da <- subset(df, subset=="train")

set.seed(12)
index <- createDataPartition(da$Condition, p=0.5, list=FALSE)
# tr <- da[index,c(2,10,11,18,19,20,25,27,28)]
# te <- da[-index,c(2,10,11,18,19,20,25,27,28)]

tr <- da[index,c(2,69,70,11,18,19,20,25,27,28)]
te <- da[-index,c(2,69,70,11,18,19,20,25,27,28)]

model <- glm(Condition ~., family = "binomial", data=tr)

cutoff <- 0.35
pred <- predict(model, newdata = te, type = "response")
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(te[-1])); confusionMatrix(pred, te$Condition, positive="HCC")

pred <- predict(model, newdata = tr, type = "response")
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(te[-1])); confusionMatrix(pred, tr$Condition, positive="HCC")

dd <- melt_roc(tibble(event = convertclass(te$Condition), GALAD = da$GALAD[-index], 
                      'Test' = predict(model, newdata = te, type = "response")),
               "event", c("GALAD", "Test"))

ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()
calc_auc(ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(event = convertclass(tr$Condition), GALAD = da$GALAD[index], 
                      'Train' = predict(model, newdata = tr, type = "response")),
               "event", c("GALAD", "Train"))

ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()
calc_auc(ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

### Train
da <- subset(df, subset == "train")
race <- pmin(100,da$tScore +
               ifelse(da$L3>15, 0.1, ifelse(da$L3>30, 1, 0)) +
               ifelse(da$Age>66, 0.1, 0) +
               ifelse(da$Gender=="M", 0.1, 0) +
               # ifelse(da$Group=="ASH", 0.1, 0) +
               # ifelse(grepl("ASH",da$Group), 0.1, 0) +
               0
)

dd <- melt_roc(tibble(event = convertclass(da$Condition), GALAD = da$GALAD, 'Train' = race), "event", c("GALAD", "Train"))
ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()
calc_auc(ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

cutoff <- 0.4
pred <- factor(ifelse(race >= cutoff, "HCC", "control"))
confusionMatrix(pred, da$Condition, positive="HCC")

### Test
da <- subset(df, subset == "test")
race <- pmin(100,da$tScore +
               ifelse(da$L3>15, 0.1, ifelse(da$L3>30, 1, 0)) +
               ifelse(da$Age>66, 0.1, 0) +
               ifelse(da$Gender=="M", 0.1, 0) +
               # ifelse(da$Group=="ASH", 0.1, 0) +
               # ifelse(grepl("ASH",da$Group), 0.1, 0) +
               0
)

dd <- melt_roc(tibble(event = convertclass(da$Condition), GALAD = da$GALAD, 'Test' = race), "event", c("GALAD", "Test"))
ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()
calc_auc(ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

cutoff <- 0.4
pred <- factor(ifelse(race >= cutoff, "HCC", "control"))
confusionMatrix(pred, da$Condition, positive="HCC")

### GLM
dd <- da[c(2,10,11,18,19,20,25,27,28)]
dd <- da[c(2,69,70,11,18,19,20,25,27,28)]

model <- glm(Condition ~., family = "binomial", data=dd)

cutoff <- 0.35
pred <- predict(model, newdata = dd, type = "response")
pred <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(dd[-1])); confusionMatrix(pred, dd$Condition, positive="HCC")

dd <- melt_roc(tibble(event = convertclass(da$Condition), GALAD = da$GALAD, 
                      'Model' = predict(model, newdata = dd, type = "response")),
               "event", c("GALAD", "Model"))

ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()
calc_auc(ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())


th <- 1
race <- as.factor(ifelse(0.52*da$mut_level+0.48*da$meth_level+0.52*da$afp_lev+0.48*da$dcp_lev-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

confusionMatrix(da$predClass, da$Condition, positive="HCC")

galad <- da$GALAD + 7.2*da$Mutation.level + 13.3*da$methylation.prediction.score
# cutoff <- ksplot(rocit(galad, convertclass(da$Condition)))$`KS Cutoff`
cutoff <- 6
galad <- as.factor(ifelse(galad > cutoff, "HCC", "control"))
confusionMatrix(galad, da$Condition, positive="HCC")

galad2 <- da$GALAD + 7.2*da$Mutation.level + 13.3*da$methylation.prediction.score
rocplot <- ggplot(da, aes(m = L3, d = Condition))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") + annotate("text", x = .75, y = .25, 
    label = paste("AUC =", round(calc_auc(rocplot)$AUC, 3)))

dd <- melt_roc(tibble(event = convertclass(da$Condition), GALAD = da$GALAD, 
                      'GALADMM' = c(da$GALAD + 7.2*da$Mutation.level + 13.3*da$methylation.prediction.score)),
               "event", c("GALAD", "GALADMM"))

ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()
calc_auc(ggplot(dd, aes(d = D.event, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

### Cross-validation
da <- subset(df, subset != "out")
dd <- da[c(2,10,11,18,19,20,25,27,28)]
dd <- da[c(2,69,70,11,18,19,20,25,27,28)]

ctrl <- trainControl(method = "cv", number = 5)
parModel <- train(Condition ~., data = dd, method = "glmnet", family = 'binomial', trControl = ctrl)
dV <- dummyVars(Condition ~., data = dd)

do { 

  library(tidyverse)  # for data manipulation
  library(dlstats)    # for package download stats
  library(pkgsearch)  # for searching packages
  
  rocPkg <-  pkg_search(query="ROC",size=200)
  
  rocPkgShort <- rocPkg %>% 
    filter(maintainer_name != "ORPHANED", score > 190) %>%
    select(score, package, downloads_last_month) %>%
    arrange(desc(downloads_last_month))
  head(rocPkgShort)
  
  library(dlstats)
  shortList <- c("pROC","precrec","ROCit", "PRROC","ROCR","plotROC")
  downloads <- cran_stats(shortList)
  ggplot(downloads, aes(end, downloads, group=package, color=package)) +
    geom_line() + geom_point(aes(shape=package)) +
    scale_y_continuous(trans = 'log2')
  
table(da$Condition)

th <- 1
race <- as.factor(ifelse(0.52*da$mut_level+0.48*da$meth_level-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

race <- as.factor(ifelse(0.52*da$mut_level+0.48*da$meth_level+0.52*da$afp_lev+0.48*da$dcp_lev-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

cutoff <- -0.52
galad <- as.factor(ifelse(da$GALAD>cutoff, "HCC", "control"))
confusionMatrix(galad, da$Condition, positive="HCC")

clipr::write_clip(race)
clipr::write_clip(galad)

th <- ksplot(rocit(c(da$mut_level-da$sex_lev), convertclass(da$Condition)))$`KS Cutoff`
race <- as.factor(ifelse(da$mut_level-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

th <- ksplot(rocit(c(da$meth_level-da$sex_lev), convertclass(da$Condition)))$`KS Cutoff`
race <- as.factor(ifelse(da$meth_level-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

th <- ksplot(rocit(c(da$meth_avg-0.2*da$RF_train), convertclass(da$Condition)))$`KS Cutoff`
race <- as.factor(ifelse(da$meth_avg-0.2*da$RF_train>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

th <- ksplot(rocit(c(0.52*da$afp_lev+0.48*da$dcp_lev-da$sex_lev), convertclass(df$Condition)))$`KS Cutoff`
race <- as.factor(ifelse(0.52*da$afp_lev+0.48*da$dcp_lev-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

th <- ksplot(rocit(c(2.34*da$lnAFP+1.33*da$lnDCP), convertclass(df$Condition)))$`KS Cutoff`
race <- as.factor(ifelse(2.34*da$lnAFP+1.33*da$lnDCP>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

df$Stage[is.na(df$Stage)] <- "X"
df$Stage[df$Stage=="B"] <- "C"

da <- subset(df, Stage != "C")
da <- subset(df, Stage == "C" | Stage == "X")
da <- subset(df, Group == "ASH")
da <- subset(df, Group == "NASH")
da <- subset(df, Group == "VIR")

da <- subset(df, subset == "train")
da <- subset(df, subset == "test")
da <- subset(df, subset != "out")

da <- subset(df, subset != "out" & Condition=="HCC" & tumor_size<2)
da <- subset(df, subset != "out" & Condition=="HCC" & tumor_size>=2 & tumor_size<3)
da <- subset(df, subset != "out" & Condition=="HCC" & tumor_size>=3 & tumor_size<5)
da <- subset(df, subset != "out" & Condition=="HCC" & tumor_size>=5)

da <- subset(df, subset != "out" & Condition=="HCC" & Stage == "0")
da <- subset(df, subset != "out" & Condition=="HCC" & Stage == "A")
da <- subset(df, subset != "out" & Condition=="HCC" & Stage == "B")
da <- subset(df, subset != "out" & Condition=="HCC" & Stage == "C")

da <- subset(df, subset != "out" & Gender == "M")
da <- subset(df, subset != "out" & Gender == "F")

da <- subset(df, subset != "out" & Age < 50)
da <- subset(df, subset != "out" & Age >= 50 & Age < 60)
da <- subset(df, subset != "out" & Age >= 60 & Age < 70)
da <- subset(df, subset != "out" & Age >= 70)

da <- subset(df, subset != "out" &  Race == "H")
da <- subset(df, subset != "out" &  Race == "NH")
da <- subset(df, subset != "out" &  Race == "A")
da <- subset(df, subset != "out" &  Race == "B")

da <- subset(df, subset != "out" &  Group == "ASH")
da <- subset(df, subset != "out" &  Group == "NASH")
da <- subset(df, subset != "out" &  Group == "OTHER")
da <- subset(df, subset != "out" &  Group == "VIR")

galad <- as.factor(ifelse(da$GALAD>cutoff, "HCC", "control"))
confusionMatrix(galad, da$Condition, positive="HCC")

# cutoff <- ksplot(rocit(da$GALAD, convertclass(da$Condition)))$`KS Cutoff`

th <- 1
race <- as.factor(ifelse(0.52*da$mut_level+0.48*da$meth_level+0.52*da$afp_lev+0.48*da$dcp_lev-da$sex_lev>=th,"HCC","control"))
confusionMatrix(race, da$Condition, positive="HCC")

confusionMatrix(da$predClass, da$Condition, positive="HCC")

cutoff <- -0.5242
galad <- as.factor(ifelse(da$GALAD>cutoff, "HCC", "control"))
confusionMatrix(galad, da$Condition, positive="HCC")

}

library(ROCit)
library(ROCR)
library(plotROC)
library(kernlab)
library(caret)
library(dplyr)
library(readxl)
library(synthpop)
library(ggfortify)
library(randomForest)

setwd("~/Proj/model")
df <- read_excel("DataBook175.xlsx", sheet="all")
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],  as.factor)

##### Main model ############################################################
#
index <- grep("train",df$Set)
varli <- c(3:7)
cutoff <- 0.47
tr <- df[index, varli]
te <- df[-index, varli]
c(cutoff, names(tr)[-1]);

model <- train(x= tr[-1] , y = tr$Condition, 
               method="glmnet",
               verbose= FALSE,
               trControl = trainControl(method = "none",classProbs = TRUE),
               tuneGrid = expand.grid(alpha= 0.7,lambda=0.01139748),
               metric= "ROC" )

cutoff <- 0.4
pred <- predict(model, newdata = df[varli], 'prob')[2]
predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
c(cutoff, names(tr)[-1]); confusionMatrix(predClass, df$Condition, positive = "HCC")
round(unlist(ciAUC(rocit(pred$HCC, df$Condition))[c(1,5,6)]), 3)

#### Test 22 Samples #####################################################
#
for(i in 1:1000) {
  set.seed(i)
  # te <- sample_n(df[-index,varli], 44)
  synData <- rbind(sample_n(subset(te, Condition=="HCC"), 22),
                   sample_n(subset(te, Condition=="control"), 22))
  pred <- predict(model, newdata = synData, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  print(unname(round(c(unlist(ciAUC(rocit(pred$HCC, synData$Condition))[1]),
                       confusionMatrix(predClass, synData$Condition, positive = "HCC")$byClass[c(11,1,2)]), 3)))
}
#### Synthetic 200 Samples #####################################################
#
synData <- te
for(i in 1:10) {
  mark <- sample(200:999,1)
  synData <- rbind(synData, cbind(Condition=df$Condition,
                                  rbind(syn(subset(df, Condition=="HCC", varli[2:8]), seed=mark)$syn,
                                        syn(subset(df, Condition=="control", varli[2:8]), seed=mark)$syn)))
}
for(i in 1:1100) {
  mark <- sample(200:99999,1)
  set.seed(mark)
  te <- sample_n(synData, 200)
  pred <- predict(model, newdata = te, 'prob')[2]
  predClass <- factor(ifelse(pred > cutoff, "HCC", "control"))
  print(unname(round(c(unlist(ciAUC(rocit(pred$HCC, te$Condition))[1]),
                       confusionMatrix(predClass, te$Condition, positive = "HCC")$byClass[c(11,1,2)]), 3)))
}

#### EOF ####

