



library(phytools)
#read the nexua file
read.nexus('Zhang time tree.nex')


#assign variable to the tree
tree = read.nexus('Zhang time tree.nex')

 
#make another tree that will the tree desired
tree2 = tree
write.tree(tree2)

#include the taxa wanted in a new variable
species = c("Agkistrodon_piscivorus", "Anomalepis_mexicanus", "Bitis_gabonica", "Boiga_irregularis", "Bothriechis_lateralis", "Bothriechis_schlegelii", "Bothrops_alternatus", "Bothrops_asper", "Bothrops_atrox", "Bothrops_colombiensis", "Bothrops_insularis", "Bothrops_jararacussu", "Bothrops_jararaca", "Bothrops_neuwiedi", "Bungarus_multicinctus", "Cerrophidion_godmani", "Crotalus_adamanteus", "Crotalus_durissus", "Crotalus_horridus", "Crotalus_simus", "Drysdalia_coronoides", "Echis_carinatus", "Echis_coloratus", "Echis_ocellatus", "Echis_pyramidum", "Gloydius_intermedius", "Hypsiglena_torquata", "Micrurus_carvalho", "Micrurus_corallinus", "Micrurus_lemniscatus", "Micrurus_paraensis", "Micrurus_spixii", "Micrurus_surinamensis", "Micrurus_altirostris", "Micrurus_fulvius", "Naja_atra", "Naja_kaouthia",  "Ovophis_okinavensis", "Pantherophis_guttatus", "Phalotris_mertensi", "Philodryas_olfersii", "Protobothrops_flavoviridis", "Sistrurus_catenatus", "Thamnodynastes_strigatus", "Opheodrys_aestivus","Dispholidus_typus", "Python_regius", "Eublepharis_macularius")

#prune the tree drop.tip func(new tree name, select the labels[ exclude the no matches(match with the species)])
pruned.tree = drop.tip(tree2, tree$tip.label[-na.omit(match(species, tree$tip.label))])


#post production of tee
write.tree(pruned.tree)

plot(pruned.tree)
is.ultrametric(pruned.tree)
is.rooted(pruned.tree)
write.nexus(pruned.tree, file = "MinSnake tree3.nex")
