library(chutils)

getPhysicalSimilarity <- function (x, y, xParity = "normal", yParity = "normal", useProportionMethod = FALSE, optimizeOneToOneMirror = FALSE) {
  lines <- c("top", "topLeft", "topRight", "middle", "bottomLeft", "bottomRight", "bottom")
  v.1 <- c(F,F,T,F,F,T,F)
  v.2 <- c(T,F,T,T,T,F,T)
  v.3 <- c(T,F,T,T,F,T,T)
  v.4 <- c(F,T,T,T,F,T,F)
  v.5 <- c(T,T,F,T,F,T,T)
  v.6 <- c(T,T,F,T,T,T,T)
  v.7 <- c(T,F,T,F,F,T,F)
  v.8 <- c(T,T,T,T,T,T,T)
  v.9 <- c(T,T,T,T,F,T,F)

  lineMatrixNormal <- data.frame(v.1, v.2, v.3, v.4, v.5, v.6, v.7, v.8, v.9)
  row.names(lineMatrixNormal) <- lines

  #create mirror matrix
  lineMatricMirror <- lineMatrixNormal
  lineMatricMirror["topLeft",] <- lineMatrixNormal["topRight",]
  lineMatricMirror["topRight",] <- lineMatrixNormal["topLeft",]
  lineMatricMirror["bottomLeft",] <- lineMatrixNormal["bottomRight",]
  lineMatricMirror["bottomRight",] <- lineMatrixNormal["bottomLeft",]

  #make sure the x and y are equal lengths
  if(length(x) != length(y)) {
    stop ("the lengths of x and y must be equal")
  }

  #get the correct comparison matricies
  if(xParity == "normal") {
    xMatrix <- lineMatrixNormal
  } else {
    xMatrix <- lineMatricMirror
  }

  if(yParity == "normal") {
    yMatrix <- lineMatrixNormal
  } else {
    yMatrix <- lineMatricMirror
  }
  ps <- NULL
  for(i in 1 : length(x)) {
    x1 <- as.data.frame (table(xMatrix[[x[i]]],yMatrix[[y[i]]]))
    x1[[1]] <- factor(x1[[1]], levels=c(TRUE, FALSE))
    x1[[2]] <- factor(x1[[2]], levels=c(TRUE, FALSE))
    nOver <- x1[x1[[1]] == T & x1[[2]] == T,"Freq"]
    nDiff <- sum(x1[x1[[1]] != x1[[2]] ,"Freq"])

    if(optimizeOneToOneMirror & x[i] == 1 & y[i] == 1 & ((yParity == "normal" & xParity == "mirror") | (yParity == "mirror" & xParity == "normal")) ) {
      nOver <- 2 # two overlapping lines
      nDiff <- 2 #  a line at the top does not overlap
    }
    if(useProportionMethod) {
      ps[i] <- (nOver/(nDiff + nOver))^2
    } else {
      ps[i] <- nOver/nDiff
    }

  }
  df.out <- data.frame(x = x, y = y, xParity, yParity, ps = ps)
  return(df.out)
}

useProportionMethod = FALSE
optimizeOneToOneMirror = TRUE
probes <- c(1,3,4,6,7,9)
#comparison <- c(1,2,3,4,5,6,7,8,9)
comparison <- c(1,3,4,5,6,7,9)

xProbes <- rep(probes, each = length(comparison))
yComparisons <- rep(comparison, times = length(probes))

df.sim.mm <- getPhysicalSimilarity(xProbes, yComparisons, xParity="mirror", yParity="mirror",useProportionMethod = useProportionMethod, optimizeOneToOneMirror = optimizeOneToOneMirror)
df.sim.mm <- df.sim.mm[df.sim.mm$x != df.sim.mm$y,]
df.sim.nm <- getPhysicalSimilarity(xProbes, yComparisons, xParity="mirror", yParity="normal",useProportionMethod = useProportionMethod, optimizeOneToOneMirror = optimizeOneToOneMirror)
df.sim.nn <- getPhysicalSimilarity(xProbes, yComparisons, xParity="normal", yParity="normal",useProportionMethod = useProportionMethod, optimizeOneToOneMirror = optimizeOneToOneMirror)
df.sim.nn <- df.sim.nn[df.sim.nn$x != df.sim.nn$y,]

df.sim.1 <- rbind(df.sim.mm,df.sim.nm,df.sim.nn)
names(df.sim.1)[names(df.sim.1) == 'x'] <- 'probe'
names(df.sim.1)[names(df.sim.1) == 'y'] <- 'comparison'

df.sim.1$comp <- paste(df.sim.1$xParity, df.sim.1$yParity, sep="-")
df.sim.1 <- df.sim.1[,c(1:2,6,5)]

#read in similarity table
#define how numbers relate to five
df.sim.1$sideOfFive <- ifelse(df.sim.1$probe > 5 & df.sim.1$comparison > 5, "sameSideOfFive", ifelse(df.sim.1$probe < 5 & df.sim.1$comparison < 5, "sameSideOfFive",  ifelse(df.sim.1$probe == 5 | df.sim.1$comparison == 5, "compareToFive",  "contraFive")) )


#find similarity function for elements greater and less than 5
df.sim <- data.frame(df.sim.1 %>% group_by(probe, comp, sideOfFive) %>% summarize (mPS = mean(ps, na.rm=T)))
#create a wide dataset
df.sim <- df.sim  %>% spread(sideOfFive, mPS)
#cacluate similarity of items on the opposite side of 5 relative to those on the same side of 5.  This
#should tell us something about the likelihood of confusion with other numbers in the magnitude comparison task,
df.sim$simRelativeToFive <- df.sim$contraFive/df.sim$sameSideOfFive

for(i in 1:length(df.sim$probe)) {
  df.sim$welford[i] <- ifelse(df.sim$probe[i] > 5, log(df.sim$probe[i]/(df.sim$probe[i]-5)), log(5/(5-df.sim$probe[i])))
}

### rearrange to get comparison to 5
df.sim.wel <- unique(df.sim[,c("probe","welford")])
df.sim.compTo5 <- unique(df.sim[,c("probe","comp","compareToFive")])
df.sim.compTo5 <- df.sim.compTo5 %>% spread(comp, compareToFive)
names(df.sim.compTo5)[names(df.sim.compTo5) == 'mirror-mirror'] <- 'mm.compTo5'
names(df.sim.compTo5)[names(df.sim.compTo5) == 'normal-normal'] <- 'nn.compTo5'
names(df.sim.compTo5)[names(df.sim.compTo5) == 'mirror-normal'] <- 'nm.compTo5'

df.simAll <- merge(df.sim.wel, df.sim.compTo5)

### rearrange to get relative to 5
df.sim.relTo5 <- unique(df.sim[,c("probe","comp", "simRelativeToFive")])
df.sim.relTo5 <- df.sim.relTo5 %>% spread(comp, simRelativeToFive)
names(df.sim.relTo5)[names(df.sim.relTo5) == 'mirror-mirror'] <- 'mm.relTo5'
names(df.sim.relTo5)[names(df.sim.relTo5) == 'normal-normal'] <- 'nn.relTo5'
names(df.sim.relTo5)[names(df.sim.relTo5) == 'mirror-normal'] <- 'nm.relTo5'

df.simAll <- merge(df.simAll, df.sim.relTo5)

#get similarity to self vs other
df.sim.2 <- df.sim.1[df.sim.1$comp == "mirror-normal", ]
df.sim.2$sideOfFive <- ifelse(df.sim.2$probe == df.sim.2$comparison,"self", "other")
df.sim.2 <- data.frame(df.sim.2 %>% group_by(probe, comp, sideOfFive) %>% summarize (mPS = mean(ps, na.rm=T)))
df.sim.2 <- df.sim.2  %>% spread(sideOfFive, mPS)
df.sim.2$nm.relSelf <- df.sim.2$self/df.sim.2$other

df.simAll <- merge(df.simAll, df.sim.2[,c("probe", "self", "nm.relSelf")])
names(df.simAll)[names(df.simAll) == 'self'] <- 'simSelf'

write.table(df.simAll, file="PhysicalSimilarity.txt", quote=F, sep="\t", row.names=F)

cexSize <- 1.5
pdf("physical similarity functions.pdf", width = 10, height = 7)
  op <- par(mfcol = c(2,2), bty="n")
  yAxis <- c(0,ceiling(with(df.simAll, max (welford))))
  with(df.simAll, plot(welford ~ probe, ylim = yAxis, xlim = c(0,10), pch=16, cex = cexSize, main = "Welford"))
  yAxis <- c(0,ceiling(with(df.simAll, max (nn.relTo5))))
  with(df.simAll, plot(nn.relTo5 ~ probe, ylim = yAxis, xlim = c(0,10), pch=16, cex = cexSize, main = "Normal - Normal Similarity Relative to Five"))
  yAxis <- c(0,ceiling(with(df.simAll, max (simSelf))))
  with(df.simAll, plot(simSelf ~ probe, ylim = yAxis, xlim = c(0,10), pch=16, cex = cexSize, main = "Mirror - Normal Similarity Relative to Self"))
  yAxis <- c(0,ceiling(with(df.simAll, max (nm.relTo5))))
  with(df.simAll, plot(nm.relTo5 ~ probe, ylim = yAxis, xlim = c(0,10), pch=16, cex = cexSize, main = "Mirror - Normal Similarity Relative to Five"))
dev.off()
par(op)
