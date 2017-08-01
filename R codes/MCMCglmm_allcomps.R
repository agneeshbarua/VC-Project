library(phytools)
library(ape)
library(MCMCglmm)


setwd("/Users/agneesh-barua/Dropbox (OIST)/R Training")
list.files()
snphylo = read.nexus("MinSnake_tree_u.nex")
plot(snphylo)
is.ultrametric(snphylo)

phydata = read.csv("snakedata.csv", header = TRUE)
head(phydata)

inv.phylo = inverseA(snphylo, nodes = "TIPS", scale =TRUE)

priors1 <- list(G=list(G1=list(V=1, nu=0.002)), R=list(V=diag(10), nu= 0.002)) # standard prior
cov_model<-MCMCglmm(cbind(TFTx,BPP,CRISP,CTL,GF,OHA,LAAO,PLA2,SVMP,SVSP) ~ cm, 
                    random = ~species,
                    rcov = ~us(trait):units,
                    family=rep("gaussian",10), 
                    ginverse=list(species=inv.phylo$Ainv), 
                    prior=priors1,
                    data=phydata, 
                    nitt=500000, burnin=500,thin=200)


summary(cov_model)
plot(cov_model)
lambda.PLA2 <- cov_model$VCV[,'species']/(cov_model$VCV[,'species']+cov_model$VCV[,'traitPLA2:traitPLA2.units'])
mean(lambda.PLA2)
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

q=matrix(summary(cov_model)$Rcovariances[,1],10,10)


rownames(q) =c("TFTx","BPP","CRISP","CTL","GF", "OHA","LAAO","PLA2","SVMP","SVSP")
colnames(q) =c("TFTx","BPP","CRISP","CTL","GF", "OHA","LAAO","PLA2","SVMP","SVSP")


q

cvm = cov(q)
cvm


library(reshape2)
library(ggplot2)

cvm<-melt(cvm)

cvm= cvm[cvm$value!=0,]

ggplot(cvm, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=("orangered"), mid = "white", high=("orangered4"), midpoint = 0, space = "Lab", guide = "colourbar") +
  labs(x="Components", y="Components", title="Covariance between venom components") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))



