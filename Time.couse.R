## time couse code
setwd("C:/Users/xli/Desktop/TimeCouse/TimeCouse_MA_0606/RAWgene")
raw_3d7 <- read.csv("Raw_3D7.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)
filter1_3d7 = rowSums(raw_3d7 > 10) >=4
countMatrixFiltered_3d7 = raw_3d7[filter1_3d7, ]

## find genes expressed in all samples
name5 = row.names(countMatrixFiltered_3d7)
common2 <- intersect(name5, name6)
Commonfilter1_3d7 <- countMatrixFiltered_3d7[common2,]

##########################################################################################
## WGCNA
wgcnaInputMA <- cbind(Commonfilter1_3d7, Commonfilter1_mal12)
row.names(wgcnaInputMA) <- common2
write.table(wgcnaInputMA,"wgcnaInputMA.csv", row.names=common2, sep=",")

library(WGCNA)

data<-read.csv('wgcnaInputMA.csv',header=T)
datExpr=as.data.frame(t(data[,-1]))
names(datExpr)=data[,1]
rownames(datExpr)=names(data)[-1]


dim(datExpr)
[1]  114 4010

# Sample clustering to detect outliers
sampleTree = hclust(dist(datExpr), method = "average");
pdf('sampleTree_MA.pdf', width=14, height=7)
par(cex = 0.7)
par(mar = c(0,6,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Scale-free topology (sft) fit index as a function of the soft-thresholding power
powers = c(1:20)
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",verbose = 5, corFnc='bicor')
sft
pdf('sft.pdf',width=12,height=7)
par(mfrow = c(1,2));
cex1 = 0.6
plot(sft$fitIndices[,1], sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Gene dendrogram and module colors
softPower = 15
adjacency = adjacency(datExpr, power = softPower,type='signed')
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize =100
detectCutHeight=0.995
dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=detectCutHeight, deepSplit =F, minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
pdf('module_MA_0.995.pdf',width=12,height=7)
plotDendroAndColors(geneTree, dynamicColors, rowText=dynamicColors, dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors", cex.rowText=0.5)
dev.off()

# generate toptable加上P值
nSamples = nrow(datExpr);
probes = names(datExpr)
connectivity=intramodularConnectivity(adjacency, dynamicColors)
geneInfo0 = data.frame(probes,moduleColor = dynamicColors, connectivity[,2])
write.csv(geneInfo0,file='geneInfo_MA.csv',row.names=F)

MEs= moduleEigengenes(datExpr, dynamicColors)
modNames = substring(names(MEs$eigengenes), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs$eigengenes, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneInfo1 = data.frame(probes,moduleColor = dynamicColors, connectivity[,1:2], geneModuleMembership, MMPvalue)
write.csv(geneInfo1,file='geneInfo1_MA.csv',row.names=F)


##########################################################################################

ages.3d7 <- c(2,6,10,10,14,14,18,18,22,22,26,26,30,30,34,34,38,38,42,42,46,46,50,50,54,54)
ages.mal53 <- c(2,2,6,6,10,10,14,14,18,18,22,22,26,26,30,30,34,34,38,38,42,42,46,46,50,50,54,54)

# choose genes wuth kTotal > 150

filter <- read.csv("kTotal150_MA.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)
setwd("C:/Users/xli/Desktop/TimeCouse/TimeCouse_MA_0606/RAWgene/MA_plot")
filterName <- row.names(filter)

# plot gene expressions 

for(i in 1:length(filterName))
	{   
		name <- filterName[i]
		dd3d7=smooth.spline(ages.3d7, raw_3d7[name,])
		ddmal53=smooth.spline(ages.mal53, raw_mal53[name,])
	
		max1 <- max(max(predict(dd3d7,N)$y), max(predict(ddmal53,M)$y))
		max2 <- max(max(raw_3d7[name,]), max(raw_mal53[name,]))
		max <- max(max1,max2)
		
		png(paste(name,".png",sep=""))
		# plot raw point
		plot(ages.3d7, raw_3d7[name,], col =  "#FF0000FF", cex=1.2, xlab="Time(h)",  ylab = rownames(raw_3d7)[name], ylim=range(0:max))
		points(ages.mal53, raw_mal53[name,], col = "#CC00FFFF", cex=1.2)
	
		legend("topright", c("3D7","MAL53"), cex=1, fill=rainbow(2));
		
		# plot predicted lines
		lines(predict(dd3d7,N),col = "#FF0000FF", lwd=2)
		lines(predict(ddmal53,M),col = "#CC00FFFF")
	
		dev.off()
	}

##########################################################################################
# shift
filter<-read.csv("kTotal150_MA.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)
rname <- row.names(filter)
shift_3d7 <- raw_3d7[rname,]
shift_mal53 <- raw_mal53[rname,]

write.table(shift_3d7,"shift_3d7.csv", row.names=rname, sep=",")
write.table(shift_mal47,"shift_mal47.csv", row.names=rname, sep=",")

ages.3d7 <- c(2,6,10,10,14,14,18,18,22,22,26,26,30,30,34,34,38,38,42,42,46,46,50,50,54,54)
ages.mal53 <- c(2,2,6,6,10,10,14,14,18,18,22,22,26,26,30,30,34,34,38,38,42,42,46,46,50,50,54,54)

N <- seq(2,56,1)
M <- seq(2,54,0.125)

sh_mal53 <- shift(shift_3d7, shift_mal53, ages.3d7, ages.mal53, M, N)
boxplot(sh_mal12, use.cols = TRUE, outline = FALSE, ylab ="TimeShift", xlab="Time(h)", col=("gold"), main="MAL12")

abline(h = 0, col = "red", lty = 2)


##########################################################################################

# plot gene expressions after shift

ages2.3d7 <- c(2,6,10,10,14,14,18,18,22,22,26,26,30,30,34,34,38,38,42,42,46,46,50,50,54,54)
ages2.mal53 <- c(2.6495,2.6495,6.7518,6.7518,10.636,10.636,14.493,14.493,17.97391,17.97391,21.92462,21.92462,26.1629,26.1629,30.206,30.206,34.3146,34.3146,38.5465,38.5465,42.4546,42.4546,46.282,46.282,50.2352,50.2352,53.5958,53.5958)

filterName <- row.names(raw_3d7)


for(i in 1:length(filterName))
	{   
		name <- filterName[i]
		dd3d7=smooth.spline(ages2.3d7, raw_3d7[name,])
		ddmal53=smooth.spline(ages2.mal53, raw_mal53[name,])
	
		max1 <- max(max(predict(dd3d7,N)$y), max(predict(ddmal53,M)$y))
		max2 <- max(max(raw_3d7[name,]), max(raw_mal53[name,]))
		max <- max(max1,max2)
		
		png(paste(name,".png",sep=""))
		plot(ages2.3d7, raw_3d7[name,], col =  "#FF0000FF", cex=1.2, xlab="Time(h)",  ylab = rownames(raw_3d7)[name], ylim=range(0:max))
		points(ages2.mal53, raw_mal53[name,], col = "#CC00FFFF", cex=1.2)
	
		
		legend("topright", c("3D7","MAL12","MAL39","MAL47","MAL53"), cex=1, fill=rainbow(5));
		
		lines(predict(dd3d7,N),col = "#FF0000FF", lwd=2)
		lines(predict(ddmal53,M),col = "#CC00FFFF")
	
		dev.off()
	}

##########################################################################################
# limma

name <- row.names(raw_3d7)
geneM <- cbind(raw_3d7,raw_mal53)

targets <- read.csv("targets.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)
row.names(geneM) <- name
use = rowSums(geneM > 5) >=6
countMatrixFiltered = geneM[ use, ]
X <- ns(targets$Time, df = 5)
Group <- factor(targets$Group)
design <- model.matrix(~Group*X)
fit <- lmFit(countMatrixFiltered, design)
fit <- eBayes(fit)
t <- topTable(fit,coef=8:12,number=Inf)
write.table(t,"mal12vsmal47.txt",sep=" ")


##########################################################################################
limma not working very well, too much false positives
will use auc from package MESS to caculate area under the curve (AUC)

P <- seq(2,54,0.125)
name <- row.names(raw_3d7)

geneM <- cbind(raw_3d7,raw_mal53)
row.names(geneM) <- name
use = rowSums(geneM > 5) >=6
countMatrixFiltered = geneM[use, ]
filterName <- row.names(countMatrixFiltered)

A.3d7.mal12 <- matrix(0, ncol= 2, nrow = length(filterName))
for(i in 1:length(filterName))
{
	name <- filterName[i]
	dd = smooth.spline(ages2.3d7, raw_3d7[name,])
	pp = smooth.spline(ages2.mal12, raw_mal12[name,])
	mm <- auc(P, predict(dd,P)$y)
	nn <- auc(P, predict(pp,P)$y)
	p.3d7.mal12[i,1] = mm
	p.3d7.mal12[i,2] = nn
}
write.table(A.3d7.mal12,"A.3d7.mal12", col.names = c("3d7","mal12"), row.names = filterName, sep = ",")
## use fold change to find DEGs

























































