library(surface)

dat <- read.csv("./Data csv/New data/surfacedata.csv", row.names = 1)
tree <- read.nexus("MinSnake_tree_u.nex")

tree <-nameNodes(tree)
olist <-convertTreeData(tree,dat)
otree <-olist[[1]]; odata<-olist[[2]]

fwd <- surfaceForward(otree, odata, aic_threshold = 0, exclude = 0, verbose = F, plotaic = F)

k<-length(fwd)
k

fwd[[k]]

fsum <-surfaceSummary(fwd)
names(fsum)

fsum$aics

surfaceTreePlot(tree, fwd[[k]], labelshifts = T, convcol = F)
surfaceTraitPlot(dat, fwd[[k]], convcol = F)

bwd <- surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, only_best = T,verbose = F,plotaic = F)

saveRDS(bwd, "backward_surface.rds")

bsum<- surfaceSummary(bwd)
 kk<-length(bwd)
 kk
 
 bsum$alpha
 bsum$sigma_squared
 bsum$theta
 bsum$n_regimes
 
 surfaceTreePlot(tree, bwd[[kk]], labelshifts = T)
 
 surfaceAICPlot(fwd, bwd)
 par(mfrow =c(1,2),mai=c(0.8,0.8,0.2,0.2))
 surfaceTraitPlot(dat, bwd[[kk]],whattraits = c(3,4))
 surfaceTraitPlot(dat, bwd[[kk]],whattraits = c(1,2))
 surfaceTraitPlot(dat, bwd[[kk]],whattraits = c(2,3)) 
 surfaceTraitPlot(dat, bwd[[kk]],whattraits = c(1,2))

 
 
#Simulating data

#Basic simulation
newsim<- surfaceSimulate(tree, type = "hansen-fit", hansenfit = fwd[[k]]$fit, shifts = fwd[[k]]$savedshifts, sample_optima = T)
par(mfrow =c(1,2),mai=c(0.8,0.8,0.2,0.2))
surfaceTraitPlot(newsim$data, newsim, whattraits = c(1,2), convcol = F)
surfaceTraitPlot(newsim$data, newsim, whattraits = c(3,4), convcol = F)

#Using output of Basic simulation
newout<-runSurface(tree, newsim$data, only_best = T)
newsimres<-newsim$data
newsum<-surfaceSummary(newout$bwd)
newsum$n_regimes
bsum$n_regimes

#Loop simulation
result<-replicate(10, surfaceSimulate(tree, type = "hansen-fit", hansenfit = fwd[[k]]$fit, shifts = fwd[[k]]$savedshifts, sample_optima = T))
result[1,]
resultsim <-as.list.data.frame(result[1,])
resultsim

#Backward analysis with Loop simulation output
sim<-lapply(resultsim, runSurface, tree= read.nexus("MinSnake_tree_u.nex"), only_best=T)

#Obtaining output from sim. #Need to write cleaner code.
surfsum1<-surfaceSummary(sim[[1]]$bwd)
surfsum2<-surfaceSummary(sim[[2]]$bwd)
surfsum3<-surfaceSummary(sim[[3]]$bwd)
surfsum4<-surfaceSummary(sim[[4]]$bwd)
surfsum7<-surfaceSummary(sim[[7]]$bwd)
surfsum1$n_regimes
surfsum2$n_regimes  
surfsum3$n_regimes
surfsum4$n_regimes
surfsum7$n_regimes




