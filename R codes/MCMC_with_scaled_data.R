library(MCMCglmmRAM)
library(phytools)
library(ape)



snphylo = read.nexus("MinSnake_tree_u.nex")
plot(snphylo)
is.ultrametric(snphylo)

phydata = read.csv("snakedata2.csv", header = TRUE)
names = phydata$species
colnames(phydata)

#Scaling the data
phydata$species = NULL
phydata = as.matrix(phydata)
phydata = scale(phydata)
phydata

rownames(phydata) = names

phydata = cbind(phydata, species = rownames(phydata))
phydata
row.names(phydata)=NULL
phydata = data.frame(phydata)
phydata

#MCMCglmm
inv.phylo = inverseA(snphylo, nodes = "TIPS", scale =TRUE)

priors1 <- list(G=list(G1=list(V=diag(10), nu=0.002)), R=list(V=diag(10), nu= 0.002)) # standard prior
ptm <- proc.time()
cov_model<-MCMCglmm(cbind(TFTx,BPP,CRISP,CTL,GF,OHA,LAAO,PLA2,SVMP,SVSP) ~ cm, 
                    random = ~us(trait):species,
                    rcov = ~us(trait):units,
                    family=rep("gaussian",10), 
                    ginverse=list(species=inv.phylo$Ainv), 
                    prior=priors1,
                    data=phydata, 
                    nitt=500000, burnin=500,thin=200)
proc.time() - ptm

summary(cov_model)

lambda.PLA2 <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitPLA2:traitPLA2.units'])
posterior.mode(lambda.PLA2)
HPDinterval(lambda.PLA2)


lambda.TFTx <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitTFTx:traitTFTx.units'])
mean(lambda.TFTx)
posterior.mode(lambda.TFTx)
HPDinterval(lambda.TFTx)

lambda.BPP <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitBPP:traitBPP.units'])
mean(lambda.BPP)
posterior.mode(lambda.BPP)
HPDinterval(lambda.BPP)

lambda.CRISP <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitCRISP:traitCRISP.units'])
mean(lambda.CRISP)
posterior.mode(lambda.CRISP)
HPDinterval(lambda.CRISP)

lambda.CTL <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitCTL:traitCTL.units'])
mean(lambda.CTL)
posterior.mode(lambda.CTL)
HPDinterval(lambda.CTL)

lambda.GF <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitGF:traitGF.units'])
mean(lambda.GF)
posterior.mode(lambda.GF)
HPDinterval(lambda.GF)

lambda.OHA <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitOHA:traitOHA.units'])
mean(lambda.OHA)
posterior.mode(lambda.OHA)
HPDinterval(lambda.OHA)

lambda.LAAO <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitLAAO:traitLAAO.units'])
mean(lambda.LAAO)
posterior.mode(lambda.LAAO)
HPDinterval(lambda.LAAO)

lambda.SVMP <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitSVMP:traitSVMP.units'])
mean(lambda.SVMP)
posterior.mode(lambda.SVMP)
HPDinterval(lambda.SVMP)

lambda.SVSP <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitSVSP:traitSVSP.units'])
mean(lambda.SVSP)
posterior.mode(lambda.SVSP)
HPDinterval(lambda.SVSP)

#Obtaining the Rcovs
q=matrix(summary(cov_model)$Rcovariances[,1],10,10)
q

rownames(q) =c("TFTx","BPP","CRISP","CTL","GF", "OHA","LAAO","PLA2","SVMP","SVSP")
colnames(q) =c("TFTx","BPP","CRISP","CTL","GF", "OHA","LAAO","PLA2","SVMP","SVSP")

#saving the Rcovs, and using this file for PCA.
range(q)
write.csv(q, file = "Rcovs.csv")



#this part deals with plotting. Not Imp fpr the PCA
library(reshape2)
library(ggplot2)

cvm<-melt(q)
cvm = data.frame(cvm)
cvm
cvm= cvm[cvm$value!=0,]

scalerange= range(cvm$value)
scalerange
gradientends= scalerange
colorends = c("blue", "white","red")
rescale = scale(cvm$value)
rescale




ggplot(cvm, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=("orangered"), mid = "white", high=("orangered4"), midpoint = 0, space = "Lab", guide = "colourbar") +
  labs(x="Components", y="Components", title="Covariance between venom components") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))







ggplot(cvm, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=rescale)) +
  scale_fill_gradient2(low=("orangered"), mid = "white", high=("orangered4"), 
                       midpoint = 0, space = "Lab", guide = "colourbar") +
  
  labs(x="Components", y="Components", title="Covariance between venom components") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))



help("scale_fill_gradient")

data = read.csv("lambda and CI.csv")
data$CI=NULL

long = data
long


sel = long$lower
seu = long$upper
sz = long$Mean

ggplot(long, aes(x =Component, y= Mean)) +
  geom_point(position=position_dodge(0.1), size=3, shape=21, fill="red")+
  xlab("Component") +
  ylab("Mean")+
  ggtitle("Covariance mean of different venom components with CI(95%)")+
  expand_limits(y=0)+
  scale_y_continuous(breaks = 0:20*0.1)+
  geom_errorbar(aes(ymin=sel, ymax=seu), colour="black", width= .1, position = position_dodge(0.1))

