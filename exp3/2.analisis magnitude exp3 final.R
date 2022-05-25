
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(rstatix)
library(lmerTest)
library(MuMIn)
library(chutils)

source("../getPhysicalSimilarity.r")
outStatFileName <- "exp1StatsOut.txt"
psFileName <- "PhysicalSimilarity.txt"

appendState <- F
errorThreshold <- 0.15
grps <- c("sn", "Numbers", "Structure")
outlierTrim <- 0.025
cvOutlierTrim <- 0.95
cexSize <- 1.5
printUnusedParallelAnalyses <- FALSE

#the analysis formulas
numDistLmeFormula.N <- as.formula("mRT~ simSelf +  (1 |sn/Numbers)")
snarcLmeFormula.N <- as.formula("dRT~ Numbers + simSelf + (1|sn)")
numDistLmeFormula.M <- as.formula("mRT~ simSelf + (1 |sn/Numbers)")
snarcLmeFormula.M <- as.formula("dRT~ Numbers + simSelf + (1|sn)")

numDistFormula.N <- as.formula("mRT ~ simSelf")
snarcFormula.N <- as.formula("dRT ~ Numbers + simSelf ")
numDistFormula.M <- as.formula("mRT ~ simSelf")
snarcFormula.M <- as.formula("dRT ~ Numbers + simSelf")

#read in data
readHeaderData <- function(filename) {
  dat <- read.csv(filename,header=T, nrows = 1)
   return(dat)
}

#Read data file into dataframe and exclude first 3 lines
readDat <- function(filename) {
  tmp.list <- readLines(filename)
  exclude3 <- tmp.list[c(-1,-2,-3)]
  dat <- read.csv(textConnection(exclude3),header=T)
  dat$response <- ifelse(dat$response == "N/A", NA, as.character(dat$response))
  dat$response <- ifelse(dat$response == "timeout", NA, as.character(dat$response))
  return(dat)
  }

#List files
raw.data = list.files("./data",pattern=".csv",full.names=T, all.files=T)

#Merge data files, add sn column
dat.all <- NULL
head.all <- NULL
for(i in raw.data) {
  df <- readDat(i)
  df.2 <- readHeaderData(i)
  sn <- str_remove(i,"./data/")
  sn <- str_remove(sn,".csv")
  df$sn <- sn
  df.2$sn <- sn
  df$trial <- seq(1,nrow(df), 1)
  tmp.cond <- df$Condition[1]
  df$trialInBlock <- ifelse(df$Condition == tmp.cond, df$trial, df$trial-(nrow(df)/2))
  df$block <- ifelse(df$trial <= (nrow(df)/2), 1, 2)

  if(is.null(dat.all)) {
    dat.all <- df
    }
  else {
    dat.all <- rbind(dat.all,df)
  }
  if("gender" %in% colnames(df.2)) {
    names(df.2)[names(df.2) == 'gender'] <- 'sex'
  }
  if(is.null(head.all)) {
    head.all <- df.2[,c("sn","age", "sex", "handedness","link")]
  } else {
    head.all <- rbind(head.all,df.2[,c("sn","age", "sex", "handedness","link")])
  }
}

head.all$link <- NULL

#Remove unnecessary columns
dat.all <- subset(dat.all, select = - c(rowNo, type, stim1,stim2, stim3,presTime,ISI,responseWindow, feedback, feedbackTime,feedbackOptions,stimFormat,random,keyboard,cursor,presTime_ms,presTime_f,presTime_fDuration))
df <- merge(dat.all, head.all)

#download file Merge
write.table(dat.all, file="MergedData.txt", quote=F, sep="\t", row.names=F)

#remove sn with error > errorThreshold
df.err.raw <- data.frame(df %>% group_by(sn) %>% summarise(pCorrect = mean(correct, na.rm=T),mRT = mean(RT, na.rm=T), sdRT = sd(RT, na.rm=T), cv = sdRT/mRT))
numRawSns <- nrow(df.err.raw)
df.err.rm <- df.err.raw[df.err.raw$pCorrect < (1-errorThreshold), ]
df.err <- df.err.raw[df.err.raw$pCorrect >= (1-errorThreshold), ]
numFinalSns <- nrow(df.err)
numRemovedSns <- numRawSns - numFinalSns

#remove 5% of subjects with greatest coefficient of variation
df.cv<- ch.filterGrpByQuantile(df.err, "cv", grpCol = NULL, lowQuantileThreshold=0, highQuantileThreshold=cvOutlierTrim)
numRemovedcvSns <- df.cv$numRemoved

numFinalSns <- nrow(df.cv$datKept)

df <- df[df$sn %in% unique(df.err$sn), ]
df <- df[df$sn %in% unique(df.cv$datKept$sn), ]

#output to file
sink(outStatFileName, appendState)
  print("############# New Run #################")
  cat("\n Number of Raw Subjects:", numRawSns, "\n\n")
  print(df.err.raw)
  cat("\n Number of Removed Subjects because of error >", errorThreshold,":",  numRemovedSns ,"\n\n")
  print(df.err.rm)
  cat("\n Number of Removed Subjects because of CV >", cvOutlierTrim, "quintile:",  numRemovedcvSns , "\n")
  print(df.cv$datRemoved)
  cat("\n Number of Final Subjects:", numFinalSns)
sink(NULL)
appendState <- TRUE


#remove RT outliers
df.outList<- ch.filterGrpByQuantile(df, "RT", grpCol = c("Structure","Numbers"), lowQuantileThreshold=outlierTrim, highQuantileThreshold=(1-outlierTrim))
df.c1 <- df.outList$datKept
removedTrials <- df.outList$pRemoved

df.final <- df.c1

#get summary stats of final subjects
df.sum <- data.frame(df.final %>% group_by(sn) %>% summarise(mRT = mean(RT, na.rm=T), sdRT = sd(RT, na.rm=T), pCorrect = mean(correct, na.rm=T), cv = sdRT/mRT))

sink(outStatFileName, appendState)
  cat("\n Final Subjects:\n\n")
  print(df.sum)
  cat("\n\n Proportion removed trials because RTs exceeded criterion for each subject:", removedTrials, "\n\n")
sink(NULL)

#read in similarity information and merge with data
df.simAll <- read.table(psFileName, header=T, sep="\t")
df.final <- merge(df.final, df.simAll, by.x = c("Numbers"), by.y = c("probe"))
df.final$simTo5 <- ifelse(df.final$Structure == "NORMAL", df.final[["nn.compTo5"]], df.final[["nm.compTo5"]])
df.final$simRel5 <- ifelse(df.final$Structure == "NORMAL", df.final[["nn.relTo5"]], df.final[["nm.relTo5"]])

#generate df only with correct (df.c) and only with errors (df.err)
df.c <- df.final[df.final$correct==1,]
df.error <- df.final[df.final$correct==0,]

#similarity values subject, number, condition and structure
df.hand <- data.frame(df.c %>% group_by(sn, Numbers, Hand, Structure, Magnitude,order) %>% summarise(mRT = mean(RT, na.rm=T), wel = mean(welford, na.rm=T), psFive = mean(simTo5, na.rm=T), relFive = mean(simRel5, na.rm=T), relSelf = mean(nm.relSelf, na.rm=T), simSelf = mean(simSelf, na.rm=T)))

#run main mixed model analyses by structure
structs <- unique(df.hand$Structure)
df.LR <- NULL
df.hand.out <- NULL
for(i in structs) {
  df.hand.tmp <- df.hand[df.hand$Structure == i,]
  if(i == "NORMAL") {
    numDistLmeFormula = numDistLmeFormula.N
    snarcLmeFormula = snarcLmeFormula.N
  } else {
    numDistLmeFormula = numDistLmeFormula.M
    snarcLmeFormula = snarcLmeFormula.M
  }
  #run lmr for structure, welford, and similarity relative to 5
  df.wel.lmer <- lmer(numDistLmeFormula, control = lmerControl(check.conv.grad = .makeCC("warning", tol = 2e-2, relTol = NULL),optimizer ="Nelder_Mead"), data=df.hand.tmp)

  df.hand.tmp$res.RT <- df.hand.tmp$mRT

  #from long to wide (Hand) with residuals as DV
  df.LR.tmp <- spread(df.hand.tmp[, -which(names(df.hand.tmp) == "mRT")], Hand, res.RT)

  #ADD dRT (R - L)
  df.LR.tmp$dRT <- df.LR.tmp$Right - df.LR.tmp$Left

  df.snarc.lmer <- lmer(snarcLmeFormula, control = lmerControl(check.conv.grad = .makeCC("warning", tol = 2e-2, relTol = NULL),optimizer ="Nelder_Mead"), data=df.LR.tmp)

  sink(outStatFileName, appendState)
    cat("\n\n\n#############################",i,"#############################","\n")
    cat("\n\n ############### Mixed Model Regression assessing Similarity to Self Mirror vs Normal ###############\n\n")
    cat("\n ############ Full Similarity to Self Model ############n\n")
    cat("\n ##### summary output #####\n\n")
    print(summary(df.wel.lmer))
    cat("\n ##### anova output #####\n\n")
    print(anova(df.wel.lmer))
    cat("\n ## r-square ##\n\n")
    print(r.squaredGLMM(df.wel.lmer))

    if(printUnusedParallelAnalyses) {
      cat("\n\n\n ############### Mixed Model Regression assessing Snark ###############n\n")
      cat("\n ############ Full Snarc Model ############n\n")
      cat("\n ##### summary output #####\n\n")
      print(summary(df.snarc.lmer))
      cat("\n ##### anova output #####\n\n")
      print(anova(df.snarc.lmer))
      cat("\n ## r-square ##\n\n")
      print(r.squaredGLMM(df.snarc.lmer))
    }
  sink(NULL)
  df.LR <- ch.rbind(df.LR, df.LR.tmp)
  df.hand.out <- ch.rbind(df.hand.out,df.hand.tmp)
}
df.hand <- df.hand.out

#Plot numerical distance and snarc effect for each Subject
structs <- unique(df.LR$Structure)
subs <- unique(df.LR$sn)
df.slopes.n <- NULL
df.slopes.m <- NULL
for(i in subs) {
  for(j in structs) {

    struct <- j
    order <- unique(df.LR[df.LR$sn == i,"order"])

    df.lr.tmp <- df.LR[df.LR$sn == i & df.LR$Structure == j,]
    df.hand.tmp <- df.hand[df.hand$sn == i & df.hand$Structure == j,]

    if(struct == "NORMAL") {
      numDistFormula = numDistFormula.N
      snarcFormula = snarcFormula.N
    } else {
      numDistFormula = numDistFormula.M
      snarcFormula = snarcFormula.M
    }
    numDist.lm <- lm(numDistFormula, data =df.hand.tmp)
    snarc.lm <- lm(snarcFormula, data = df.lr.tmp)
    snarc.cr <- with(df.lr.tmp, cor(dRT,Numbers))
    numDist.cr <- with(df.hand.tmp, cor(mRT,wel))

    df.tmp.1 <- data.frame (sn = i, Structure = struct, order = order, wel_r2 = summary(numDist.lm)$r.squared, wel_r = numDist.cr, snc_r2 = summary(snarc.lm)$r.squared, snc_r = snarc.cr)
    df.tmp.wel <- data.frame(t(coef(numDist.lm)[2:length((coef(numDist.lm)))]))
    original_cols <- colnames(df.tmp.wel)
    colnames(df.tmp.wel) <- paste("simToSelf" ,original_cols,sep="_")
    df.tmp.snarc <- data.frame(t(coef(snarc.lm)[2:length((coef(snarc.lm)))]))
    original_cols <- colnames(df.tmp.snarc)
    colnames(df.tmp.snarc) <- paste("snarc" ,original_cols,sep="_")
    df.tmp <- cbind(df.tmp.1, df.tmp.wel, df.tmp.snarc)

    if(struct == "NORMAL") {
      df.slopes.n <- ch.rbind(df.slopes.n, df.tmp)
    } else {
      df.slopes.m <- ch.rbind(df.slopes.m, df.tmp)
    }
  }
}

df.sum <- data.frame(df.c %>% group_by(sn, Structure) %>% summarise(mRT = mean(RT, na.rm=T), sdRT = sd(RT, na.rm=T), pCorrect = mean(correct, na.rm=T), cv = sdRT/mRT))
df.slopes.n <- merge(df.slopes.n, df.sum[df.sum$Structure == "NORMAL",])
df.slopes.m <- merge(df.slopes.m, df.sum[df.sum$Structure == "MIRROR",])

#############  Overall Snarc (dRT) and Similrity to Self (mRT)  ################

  df.LR.MN <- df.LR

  df.num <- data.frame(df.c %>% group_by_at(grps) %>% summarise(mRT = mean(RT, na.rm=T), wel = mean(welford, na.rm=T), psFive = mean(simTo5, na.rm=T), relFive = mean(simRel5, na.rm=T), relSelf = mean(nm.relSelf, na.rm=T), simSelf = mean(simSelf, na.rm=T)))

pdf("summary plots.pdf")
  op <- par(bty="n")

yLims.ND <- c(550, 800)
yLims.SN <- c(-50,60)
for(i in structs) {
  if(i == "NORMAL") {
    numDistFormula = numDistFormula.N
    snarcFormula = snarcFormula.N
  } else {
    numDistFormula = numDistFormula.M
    snarcFormula = snarcFormula.M
  }

  df.tmp <- data.frame(df.LR[df.LR$Structure == i,] %>% group_by_at(grps) %>% summarise(dRT = mean(dRT, na.rm=T), psFive = mean(psFive, na.rm=T), relFive = mean(relFive, na.rm=T), relSelf = mean(relSelf, na.rm=T), simSelf = mean(simSelf, na.rm=T)))
  snarc.lm <- lm(snarcFormula, data = df.tmp)

  minY <- min(df.tmp$dRT, predict(snarc.lm)) - .1*abs(min(df.tmp$dRT, predict(snarc.lm)))
  maxY <- 1.1*max(df.tmp$dRT, predict(snarc.lm))
  with(df.tmp, plot(dRT~Numbers, main=paste(i, "SNARC"), col="grey", ylim=c(minY,maxY), xlim=c(0,10),pch = 16,cex = cexSize))
  with(df.tmp, points(predict(snarc.lm)~Numbers, pch = 16, col="black", cex = cexSize))

df.tmp.1 <- data.frame(df.tmp %>% group_by(Numbers) %>% summarise(dRT = mean(dRT, na.rm=T), psFive = mean(psFive, na.rm=T), relFive = mean(relFive, na.rm=T), relSelf = mean(relSelf, na.rm=T), simSelf = mean(simSelf, na.rm=T)))
minY <- min(df.tmp.1$dRT, predict(snarc.lm, df.tmp.1)) - .1*abs(min(df.tmp.1$dRT, predict(snarc.lm, df.tmp.1)))
maxY <- 1.1*max(df.tmp.1$dRT, predict(snarc.lm, df.tmp.1))
with(df.tmp.1, plot(dRT~Numbers, main=paste(i, "SNARC"), col="grey", ylim=yLims.SN, xlim=c(0,10),pch = 16,cex = cexSize))
with(df.tmp.1, points(predict(snarc.lm, df.tmp.1)~Numbers, pch = 16, col="black", cex = cexSize))



  numDist.lm <- lm(numDistFormula, data=df.num[df.num$Structure == i,])
  minY <- min(df.num[df.num$Structure == i,"mRT"], predict(numDist.lm)) - .1*abs(min(df.num[df.num$Structure == i,"mRT"], predict(numDist.lm)))
  maxY <- 1.1*max(df.num[df.num$Structure == i,"mRT"], predict(numDist.lm))
  with(df.num[df.num$Structure == i,], plot(mRT~Numbers, main=paste(i, "Similarity to Self"), col="grey", ylim=c(minY,maxY), xlim=c(0,10),pch = 16,cex = cexSize))
  with(df.num[df.num$Structure == i,], points(predict(numDist.lm)~Numbers, pch = 16, col="black", cex = cexSize))

df.num.1 <- data.frame(df.num[df.num$Structure == i,] %>% group_by(Numbers) %>% summarise(mRT = mean(mRT, na.rm=T), wel = mean(wel, na.rm=T), psFive = mean(psFive, na.rm=T), relFive = mean(relFive, na.rm=T), relSelf = mean(relSelf, na.rm=T), simSelf = mean(simSelf, na.rm=T)))
minY <- min(df.num.1$mRT, predict(numDist.lm, df.num.1)) - .1*abs(min(df.num.1$mRT, predict(numDist.lm, df.num.1)))
maxY <- 1.1*max(df.num.1$mRT, predict(numDist.lm, df.num.1))
with(df.num.1, plot(mRT~Numbers, main=paste(i, "Similarity to Self"), col="grey", ylim=yLims.ND, xlim=c(0,10),pch = 16,cex = cexSize))
with(df.num.1, points(predict(numDist.lm, df.num.1)~Numbers, pch = 16, col="black", cex = cexSize))


  sink(outStatFileName, appendState)
    cat("\n\n\n ############### ", i, " ###############n\n")
    cat("\n --------- SNARC \n")
    print(summary(snarc.lm))
    if(printUnusedParallelAnalyses) {
      cat("\n\n --------- Similarity To Self (mRT) \n")
      print(summary(numDist.lm))
    }
  sink(NULL)
}

dev.off()


 df.hand.1 <- df.hand[,c(1:6,length(df.hand))]
 df.num.a <- data.frame(df.hand.1 %>% group_by(sn, Hand, Structure, Magnitude,order) %>% summarise(mRT = mean(res.RT, na.rm=T)))

 #from long to wide (Hand)
 df.LR.a <- spread(df.num.a, Hand, mRT)
 df.LR.a$dRT <- df.LR.a$Right - df.LR.a$Left

yLim.Wel <- c(-50,150)
yLim.Snc <- c(-225,175)
pdf("snarc boxplot.pdf", width = 10, height = 7)
  op <- par(mfcol = c(1,2), bty="n")
  for(i in structs) {
    df.tmp.1 <- NULL
    if(i == "NORMAL") {
      df.tmp <- df.slopes.n
    } else {
      df.tmp <- df.slopes.m
    }
    df.tmp.1 <- data.frame(df.tmp[ , grepl( "snarc_" , names( df.tmp ) ) ])
    df.tmp.1 <-df.tmp.1 %>% rename_all(~stringr::str_replace(.,"^snarc_",""))
    boxplot(df.tmp.1, main = paste(i,"SNARC"),  ylim = yLim.Snc)
    abline(h=0, col = "grey")

    df.tmp.1 <- NULL
    df.tmp.1 <- data.frame(df.tmp[ , grepl( "simToSelf_" , names( df.tmp ) ) ])
    df.tmp.1 <-df.tmp.1 %>% rename_all(~stringr::str_replace(.,"^simToSelf_",""))
    boxplot(df.tmp.1, main = paste(i,"simToSelf"), ylim = yLim.Wel)
    abline(h=0, col = "grey")
  }

  par(op)
  dev.off()
