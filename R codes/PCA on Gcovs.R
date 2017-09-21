library(ggplot)
library(reshape2)
library(ggrepel)
library(scatterpie)
library(factoextra)

cov_model<-readRDS("cov_model_scaled1.rds")

summary(cov_model)
cov_model$ginverse


Eig <- eigen(matrix(summary(cov_model)$Gcovariances,10,10))
evectors <- Eig$vectors
rownames(evectors) <- colnames(phydata)[2:11]
print(round(Eig$values/sum(Eig$values) * 100, digits = 2))
print(round(cumsum(Eig$values)/sum(Eig$values) * 100, digits = 2))

scores <- as.matrix(phydata.norm[,2:11]) %*% evectors # calculating scores from eigenanalysis results

correlations <- cor(scores, phydata.norm[,2:11])

quartz(height=7, width=7)
plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",
     type="n", xlim=c(min(scores[,1:2]), max(scores[,1:2])),
     ylim=c(min(scores[,1:2]), max(scores[,1:2])))
arrows(0,0,evectors[,1],evectors[,2], length=0.1,
       angle=20, col="red")
text(evectors[,1]*1.2,evectors[,2]*1.2,
     rownames(evectors), col="red", cex=0.7)
text(scores[,1],scores[,2], rownames(scores), col="blue",
     cex=0.7)

quartz(height=7, width=7)
plot(scores[,2], scores[,3], xlab="PCA 2", ylab="PCA 3",
     type="n", xlim=c(min(scores[,2:3]), max(scores[,2:3])),
     ylim=c(min(scores[,2:3]), max(scores[,2:3])))
arrows(0,0,evectors[,2],evectors[,3], length=0.1,
       angle=20, col="red")
text(evectors[,2]*1.2,evectors[,3]*1.2,
     rownames(evectors), col="red", cex=0.7)
text(scores[,2],scores[,3], rownames(scores), col="blue",
     cex=0.7)

plotData <- phydata[,c("PLA2","SVMP","SVSP","TFTx","family","subfamily","species")]
plotData$PCA1 <- scores[,1]
plotData$PCA2 <- scores[,2]
plotData$shortname <- gsub("([A-Z]).*_([a-z]).*", "\\1\\2", plotData$species)
evedata <- evectors
evedata = as.data.frame(evedata)

evedata
rownames(evedata)
evedata$components = c("TFTx","BPP","CRISP", "CTL", "GF","OHA","LAAO","PLA2","SVMP","SVSP")
evedata$colors = c("grey","pink","brown","lightgoldenrod","sienna","navy","lightslategray","cyan","purple","orange")


ggplot(plotData,aes(PCA1, PCA2,label = shortname))+
  geom_scatterpie(aes(x=PCA1, y=PCA2), data=plotData, cols=c("PLA2","SVMP","SVSP","TFTx"))+
  scale_fill_manual(values=c("cyan","purple","orange","grey"))+
  geom_label_repel(aes(color = family))+
  scale_color_manual(values = c("red","green","blue"))+
  geom_segment(aes(0,0, xend = evedata[,1]*2, yend = evedata[,2]*2), color= evedata$colors , data = evedata,arrow = arrow(), size = 1, show.legend = F ,inherit.aes = F)+
  theme_minimal()




  
  


    
   help("geom_segment")
  help("aes")
traceback()
  


  theme_bw()+geom_label_repel()

g = ggplot(as.data.frame(evedata), aes(evedata[,1],evedata[,2]))
Venom_components<-factor(rownames(evedata))

g + geom_segment(aes(0,0, xend = evedata[,1], yend = evedata[,2], color = Venom_components), arrow = arrow(), size = 0.7)+
    geom_label(aes(label = rownames(evedata)))

    
 scale_fill_discrete()
  geom_segment(mapping = aes(0,0,xend=evedata[,1],yend=evedata[,2]), arrow = arrow(), size =1)

  theme_bw()+geom_label_repel()
  
