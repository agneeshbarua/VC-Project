library(ggplot2)
library(devtools)
library(ggbiplot)
library(factoextra)


#PCA short, with Rcovs.csv: The file obtained from Rcovariances in MCMC_with_scaled_data.R file
mydata = read.csv("Rcovs.csv", header = TRUE, row.names = 1,sep = ",")
rownames(mydata)
mydata

mydata.pca = prcomp(mydata,retx = TRUE, center = TRUE, scale. = TRUE)
sd = mydata.pca$sdev
loadings= mydata.pca$rotation
scores = mydata.pca$x

#Simple way of getting summary and ploting using screeplot.
summary(mydata.pca)

eig.val = get_eigenvalue(mydata.pca) #obtaining the eigenvalues.
eig.val

fviz_screeplot(mydata.pca, ncp=11) #Cummilatively, PC 1,2,3 account for 92% of the variation, much better than the previous one using original data. 

fviz_screeplot(mydata.pca, ncp= 11, choice= "eigenvalue") #Any eigen value > 1 can be used as a cutoff point to determine which PCs to retain. Here we see we can potentially retain till PC 3


#biplot
ggbiplot(mydata.pca,obs.scale = 2, var.scale = 2, varname.size = 5, labels = rownames(mydata), labels.size = 3) + 
  scale_color_gradient2(name = '')+ 
  theme_grey()  
#This plot shows that BPP,CRISP,SVMP, CTL are grouped closelt together. Their association is very strong. SVSP and LAAO, are also grouped with them but with a relatively weaker association.
#This could be indicative of a particular stategy involving these 4 components. PLA2,OHA,TFTx.GF, are ore independently related. OHA and PLA2 seem to be antagonists. 





#Results of the PCA 
var = get_pca_var(mydata.pca)
var
var$contrib

# Helper function : 
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- mydata.pca$rotation
sdev <- mydata.pca$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord)

var.cos2 = var.coord^2


#colour is dependent on cos2 value.
fviz_pca_var(mydata.pca, col.var = "contrib") +
  scale_color_gradient2(low = "grey98", mid= "red", high = "blue", midpoint = 55)+
  theme_dark()

#contribution of each variable to the given PC. Or the components driving the envenomation stategy
comp.cos2 = apply(var.cos2, 2, sum)
contri = function(var.cos2,comp.cos2) {var.cos2*100/comp.cos2}
var.contri = t(apply(var.cos2, 1, contri, comp.cos2))
var.contri
range(var$contrib)

#Here we see that BPP, CTL,SVMP,CRISP, are the ones that are the major contributors to variability. Individually they contribute almost the same amount.
fviz_pca_var(mydata.pca, col.var = "contrib" , habillage ="none", label = "var", repel = T, title = "Contribuiton of various components")+
  scale_color_gradient2(low = "white", high = "#FF4500", mid = "lightcoral",midpoint = 9 , limits = c(0, 18))+
  theme_bw()
help("fviz_pca_biplot")

