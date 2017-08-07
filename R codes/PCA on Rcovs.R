library(factoextra) 

p = read.csv("Rcovs.csv")
#this is the 10x10 cov matrix obtained from Rcovariance
 
 p$X = NULL
  

myEig = eigen(p)
myEig$values #eigenvalues
myEig$vectors #eigenvectors or Loadings

sdLONG = sqrt(myEig$values)  # calculating singular values from eigenvalues

loadingsLONG = myEig$vectors

loadingsLONG
rownames(loadingsLONG) = c("TFTx","BPP","CRISP","CTL","GF", "OHA","LAAO","PLA2","SVMP","SVSP")
colnames(loadingsLONG) =c("TFTx","BPP","CRISP","CTL","GF", "OHA","LAAO","PLA2","SVMP","SVSP")
#just renaming the rows and cols for plotting


standardize <- function(x) {(x - mean(x))/sd(x)}
X <- apply(mydata, MARGIN=2, FUN=standardize)
# transforming data to zero mean and unit variance (not going to do this since we already have scaled data. I have tried it with this too, nothing changes.
scoresLONG <- X %*% loadingsLONG
# calculating scores from eigenanalysis results


plot(scoresLONG[,1], scoresLONG[,2], xlab="PCA 1", ylab="PCA 2",
       type="n", xlim=c(min(scoresLONG[,1:2]), max(scoresLONG[,1:2])),
       ylim=c(min(scoresLONG[,1:2]), max(scoresLONG[,1:2])))
arrows(0,0,loadingsLONG[,1]*10,loadingsLONG[,2]*10, length=0.1,
         angle=20, col="red")
#This gives the plot I shared.
text(loadingsLONG[,1]*10*1.2,loadingsLONG[,2]*10*1.2,
       rownames(loadingsLONG), col="red", cex=0.7)



biplot(scoresLONG[,1:2], loadingsLONG[,1:2], xlab="PCA 1",
         ylab="PCA 2", cex=0.7)

plot(myEig$values)
plot()


z = prcomp(p)
ggbiplot(z)

