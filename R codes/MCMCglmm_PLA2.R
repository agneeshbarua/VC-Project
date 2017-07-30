library(phytools)
library(ape)
library(MCMCglmm)

snphylo = read.nexus("MinSnake_tree_u")
plot(snphylo)

is.ultrametric(snphylo)


phydata = read.csv("glm trial data.csv", header = TRUE)
phydata



inv.phylo = inverseA(snphylo, nodes = "TIPS", scale =TRUE)
inv.phylo


priors<- list(G=list(G1=list(V=1, nu=0.002)), R=list(V=1, nu= 0.002))


simple_model<-MCMCglmm(ln.PLA2~ln.size, random = ~SNAKE,
                       family="gaussian",ginverse=list(SNAKE=inv.phylo$Ainv),prior=priors,
                       data=phydata,nitt=50000,burnin=500,thin=200)


summary(simple_model)

plot(simple_model)



