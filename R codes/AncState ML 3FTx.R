library(phytools)


snake.tree = read.nexus("MinSnake tree3")
snake.tree
plotTree(snake.tree)

is.binary.tree(snake.tree)


venom = read.csv("BPP ML.csv", row.names = 1)
venom


venom = as.matrix(venom)[,1]
venom



require(geiger)
name.check(phy = snake.tree,data = venom)



AS = fastAnc(snake.tree, venom, vars = TRUE, CI =TRUE)
AS

range(venom)


obj = contMap(snake.tree,venom,res = 100, plot = TRUE)
plot(obj, legend=FALSE, fsize=0.8, ylim=c(1-0.09*(Ntip(obj$tree)-1),Ntip(obj$tree)),
     mar=c(5.1,0.4,0.4,0.4))
#adding our legend
add.color.bar(50, obj$cols, title="BPP value",
              lims=obj$lims, digits=3,prompt=FALSE, x=0,
              y=1-0.08*(Ntip(obj$tree)-1), lwd=4,fsize=1,subtitle="")






