library(caret)
library(readxl)
library(dplyr)
library(plotROC)
library(ROCit)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(pheatmap)
library(dendsort)

setwd("~/Proj/model")

df <- read_excel("DataBook175.xlsx", sheet="all")
df$Stage[is.na(df$Stage)] <- "control"
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],  as.factor)

#### AUC Plot ####
for(i in c("ASH", "NASH", "VIR")){
  da <- df %>% select('Condition', 'MMCEP', 'Group') %>% subset(Group == i)
  print(round(unlist(ciAUC(rocit(da$MMCEP, da$Condition))[c(1,5,6)]), 3))
}
sapply(df[c(4:7,14)], function(x){round(unlist(ciAUC(rocit(x, df$Condition))[c(1,5,6)]), 3)})

dd <- melt_roc(tibble(cond = convertclass(df$Condition), MMAD=df$MMAD, MMP=df$MMP, MMCP=df$MMCP, MMCEP=df$MMCEP,
               MM=df$MM, MMC=df$MMC, MMCE=df$MMCE, GALAD=df$GALAD),
               "cond", c("MMAD" ,"MM", "MMC", "MMCE", "MMP", "MMCP", "MMCEP","GALAD"))
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(cond = convertclass(df$Condition), MM=df$MM, MMC=df$MMC, MMCE=df$MMCE, GALAD=df$GALAD),
               "cond", c("MM", "MMC", "MMCE", "GALAD"))
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = FALSE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(cond = convertclass(df$Condition), MMP=df$MMP, MMCP=df$MMCP, MMCEP=df$MMCEP, GALAD=df$GALAD),
               "cond", c("MMP", "MMCP", "MMCEP", "GALAD"))
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = FALSE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(cond = convertclass(df$Condition), MM=df$MM, MMC=df$MMC, MMCE=df$MMCE),
               "cond", c("MM", "MMC", "MMCE",))
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = FALSE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(cond = convertclass(df$Condition), MMCEP=df$MMCEP, GALAD=df$GALAD),
               "cond", c("MMCEP", "GALAD"))
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = FALSE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(cond = convertclass(df$Condition), "Integrated Model"=df$MMCEP, GALAD=df$GALAD, Mutation=df$Mutation,
                      Methylation=df$Methylation, CNV=df$CNV, "5' End motif"=df$`EndMotif`),
               "cond", c("Integrated Model", "GALAD", "Methylation", "Mutation", "CNV", "5' End motif"))
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = TRUE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

lali <- c("Integrated Model", "Methylation", "Mutation", "CNV", "5' End motif", "GALAD")
index <- grep("train",df$Set)
tr <- df[index,]
te <- df[-index,]

dd <- melt_roc(tibble(cond = convertclass(tr$Condition), "Integrated Model"=tr$MMCEP, 
                      Mutation=tr$Mutation, Methylation=tr$Methylation, 
                      CNV=tr$CNV, "5' End motif"=tr$`EndMotif`, GALAD=tr$GALAD), "cond", lali)
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = FALSE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

dd <- melt_roc(tibble(cond = convertclass(te$Condition), "Integrated Model"=te$MMCEP, 
                      Mutation=te$Mutation, Methylation=te$Methylation, 
                      CNV=te$CNV, "5' End motif"=te$`EndMotif`, GALAD=te$GALAD), "cond", lali)
ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(show.legend = FALSE, n.cuts = 0) + 
  style_roc() + theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15))
calc_auc(ggplot(dd, aes(d = D.cond, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc())

#### Boxplot ####
ggplot(df, aes(x = Stage, y = MMAD, col = Stage)) +   # Change color of borders
  geom_boxplot(outlier.shape = NA, width=0.5) +
  scale_fill_viridis(discrete = TRUE, alpha=0.1) +
  geom_jitter(width = 0.1, alpha=0.9) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15),
    legend.position="none",
    plot.title = element_text(size=11)
  ) + ylab("Predicted Probabilities")

ggplot(df, aes(x = Group, y = MMAD, col = Condition)) +   # Change color of borders
  geom_boxplot(outlier.shape = NA, aes(fill=Condition), width=0.5) +
  scale_fill_viridis(discrete = TRUE, alpha=0.1) +
  geom_jitter(width = 0.05, alpha=0.9) +
  scale_color_manual(values = c("HCC"="red", "control"="blue")) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(hjust=0.5, size=15),
    axis.title.y = element_text(hjust=0.5, size=15),
#    legend.position="none",
    plot.title = element_text(size=11)
  ) + ylab("Predicted Probabilities")

#### Heatmap ####
mName <- "MMCEP"
da <- df %>% arrange(across(any_of(mName), desc), .by_group = TRUE)
dd <- da %>% select(c("Mutation", "Methylation", "CNV", "EndMotif", mName)) %>% data.matrix() %>% t()
rownames(dd)[which(rownames(dd) == mName)] <- "Integrated model"
rownames(dd)[which(rownames(dd) == "EndMotif")] <- "End Motif"
colnames(dd) = paste("C", 1:ncol(dd), sep = "")
ann <- da %>% select(L3, lnDCP, lnAFP, Stage, Condition) %>% data.frame()
rownames(ann) = paste("C", 1:ncol(dd), sep = "")
palet = list(Condition = c(HCC="firebrick", control="cadetblue1"),
             Stage = c('0'="beige", A="bisque",B="coral", C="brown", control="white"))
pheatmap(dd,border_color=NA, cellheight=12, gaps_row=4, cluster_rows=F, cluster_cols=F, 
         show_colnames=F ,annotation_col=ann, annotation_colors=palet, main = paste("Model", mName))


for(i in names(df)[grep("M.*[CDEMP]$",names(df))]) {
mName <- i
mName <- "MMCEP"
da <- df %>% group_by(desc(Condition)) %>% arrange(across(any_of(mName), desc), .by_group = TRUE) %>% ungroup()
dd <- da %>% select(c("Mutation", "Methylation", "CNV", "EndMotif", mName)) %>% data.matrix() %>% t()
rownames(dd)[which(rownames(dd) == mName)] <- "Integrated model"
rownames(dd)[which(rownames(dd) == "EndMotif")] <- "End Motif"
colnames(dd) = paste("C", 1:ncol(dd), sep = "")
ann <- da %>% select(L3, lnDCP, lnAFP, Stage, Condition) %>% data.frame()
rownames(ann) = paste("C", 1:ncol(dd), sep = "")
palet = list(Condition = c(HCC="firebrick", control="cadetblue1"),
             Stage = c('0'="beige", A="bisque",B="coral", C="brown", control="white"))
pheatmap(dd,border_color=NA, cellheight=12, gaps_row=4, cluster_rows=F, cluster_cols=F, 
         show_colnames=F ,annotation_col=ann, annotation_colors=palet, main = paste("Model", mName))
}

mName <- "MMCEP"
palet <- list(Condition = c(HCC="firebrick", control="cadetblue1"),
              Stage0 = c('0'="brown1", control="white"),
              StageA = c(A="brown2", control="white"),
              StageB = c(B="brown3", control="white"),
              StageC = c(C="brown4", control="white"),
              ASH = c(ASH="blue1", OTHER="white"),
              NASH = c(NASH="blue3", OTHER="white"),
              VIR = c(VIR="blue4", OTHER="white")
)
da <- df %>% group_by(desc(Condition)) %>% arrange(across(any_of(mName), desc), .by_group = TRUE) %>% ungroup()
da$Stage0 <- ifelse(da$Stage=="0","0","control")
da$StageA <- ifelse(da$Stage=="A","A","control")
da$StageB <- ifelse(da$Stage=="B","B","control")
da$StageC <- ifelse(da$Stage=="C","C","control")
da$ASH <- ifelse(da$Group=="ASH","ASH","OTHER")
da$NASH <- ifelse(da$Group=="NASH","NASH","OTHER")
da$VIR <- ifelse(da$Group=="VIR","VIR","OTHER")
dd <- da %>% select(c("Mutation", "Methylation", "CNV", "EndMotif", mName)) %>% data.matrix() %>% t()

rownames(dd)[which(rownames(dd) == mName)] <- "Integrated model"
rownames(dd)[which(rownames(dd) == "EndMotif")] <- "End Motif"
colnames(dd) = paste("C", 1:ncol(dd), sep = "")
ann <- da %>% select(
  # Stage0, 
  # StageA,
  # StageB,
  # StageC,
  # ASH,
  # NASH,
  VIR,
  Condition) %>% data.frame()
rownames(ann) = paste("C", 1:ncol(dd), sep = "")

pheatmap(dd,border_color=NA, cellheight=12, cluster_rows=F, cluster_cols=F, 
         show_colnames=F ,annotation_col=ann, annotation_colors=palet, main = paste("Model", mName))

ann <- da %>% select(Stage0, StageA, StageB, StageB, StageC, ASH, NASH, VIR, Condition) %>% data.frame()
rownames(ann) = paste("C", 1:ncol(dd), sep = "")
pheatmap(dd,border_color=NA, cellheight=20, gaps_row=4, cluster_rows=F, cluster_cols=F,
         show_colnames=F ,annotation_col=ann, annotation_colors=palet, main = paste("Model", mName))

mName <- "MMCEP"
palet = list(Condition = c(HCC="firebrick", control="cadetblue1"),
             Stage = c('0'="beige", A="bisque",B="coral", C="brown", control="white"),
             Group = c("ASH"="aquamarine", "NASH"="aquamarine3", "VIR"="aquamarine4"))
for(i in c("0", "A", "B", "C")) {
  stage <- i
  da <- df %>% group_by(desc(Condition)) %>% arrange(across(any_of(mName), desc), .by_group = TRUE) %>% ungroup()
  dd <- da %>% subset(Stage==stage) %>% select(c("Mutation", "Methylation", "CNV", "EndMotif", mName)) %>%  data.matrix() %>% t()
  rownames(dd)[which(rownames(dd) == mName)] <- "Integrated model"
  rownames(dd)[which(rownames(dd) == "EndMotif")] <- "End Motif"
  colnames(dd) = paste("C", 1:ncol(dd), sep = "")
  ann <- da %>% subset(Stage==stage) %>% select(Stage, Condition) %>% data.frame()
  rownames(ann) = paste("C", 1:ncol(dd), sep = "")
  pheatmap(dd,border_color=NA, cellheight=12, gaps_row=4, cluster_rows=F, cluster_cols=F, 
           show_colnames=F ,annotation_col=ann, annotation_colors=palet, main = paste("Model", mName))
}
for(i in c("ASH", "NASH", "VIR")) {
  group <- i
  da <- df %>% group_by(desc(Condition)) %>% arrange(across(any_of(mName), desc), .by_group = TRUE) %>% ungroup()
  dd <- da %>% subset(Group==group) %>% select(c("Mutation", "Methylation", "CNV", "EndMotif", mName)) %>%  data.matrix() %>% t()
  rownames(dd)[which(rownames(dd) == mName)] <- "Integrated model"
  rownames(dd)[which(rownames(dd) == "EndMotif")] <- "End Motif"
  colnames(dd) = paste("C", 1:ncol(dd), sep = "")
  ann <- da %>% subset(Group==group) %>% select(Group, Condition) %>% data.frame()
  rownames(ann) = paste("C", 1:ncol(dd), sep = "")
  pheatmap(dd,border_color=NA, cellheight=12, gaps_row=4, cluster_rows=F, cluster_cols=F, 
           show_colnames=F ,annotation_col=ann, annotation_colors=palet, main = paste("Model", mName))
}
graphics.off()
# ann <- df %>% select(Stage, Condition) %>% data.matrix()
# pheatmap(dd,border_color=NA, cellheight=12, gaps_row=4, cluster_rows=F,
#          cluster_cols=F, show_colnames=F ,annotation_col=ann, annotation_colors=palet)
# colL <- colorRampPalette(c("red", "orange", "blue"))(length(unique(df$L3)))

#### Sensitivity ####

cutoff <- 0.5
cList <- c(4:7) # mono models
  sapply(df[cList], function(x){round(unlist(
    confusionMatrix(factor(ifelse(x>cutoff, "HCC", "control")), df$Condition, positive="HCC")$byClass[c(1:2,11)]), 3)})

cList <- c(9,11,13,15,17,19,21) # multi models
gList <-sort(c("ASH", "NASH", "VIR"))
dd <- NULL
for(i in gList) {
  da <- subset(df, Group == i)
  dd <- rbind(dd, sapply(da[cList], function(x){round(unlist(
    confusionMatrix(x, da$Condition, positive="HCC")$byClass[c(1:2,11)]), 3)}))
}
rownames(dd) <- paste(sort(c(paste("Etiology", gList), paste("Etiology", gList), paste("Etiology", gList))), rownames(dd))
colnames(dd) <- sub("pred","",colnames(df[cList]))
dd

cList <- c(9,11,13,15,17,19,21)
gList <- sort(c("0", "A", "B", "C"))
dd <- NULL
for(i in gList) {
  da <- subset(df, Stage == i)
  dd <- rbind(dd, sapply(da[cList], function(x){round(unlist(
    confusionMatrix(x, da$Condition, positive="HCC")$byClass[1]), 3)}))
}
rownames(dd) <- paste(sort(c(paste("Stage", gList, "Sensitivity"))))
colnames(dd) <- sub("pred","",colnames(df[cList]))
dd

#### Final Model ####
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
pred <- predict(model, newdata = df[varli], 'prob')
predClass <- factor(ifelse(pred$HCC > cutoff, "HCC", "control"))
c(cutoff, names(tr)[-1]); confusionMatrix(predClass, df$Condition, positive = "HCC")

clipr::write_clip(cbind(pred, predClass))
round(unlist(ciAUC(rocit(pred$HCC, df$Condition))[c(1,5,6)]), 3)
