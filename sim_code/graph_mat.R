# generate neighborhood matrix for CNN

# family structure
m = 3
m2 = 2
ped.struct = data.frame(index=1:27, 
                      index.mother=c(2, 4, 6, rep(27, 4), rep(2, 5), rep(1, 4), rep(4, 5), rep(6, 5), 27), 
                      index.father=c(3, 5, 7, rep(27, 4), rep(3, 5), rep(27, 4), rep(5, 5), rep(7, 5), 27), 
                      gender=c(0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0))
fam.size = nrow(ped.struct)

# neighborhood size
nh.size = 3 + 2*m + 2*m2

# matrix where row i contains neighbors of family member i (self, mother, father, sisters, brother, daughters, sons)
graph.mat =  matrix(rep(NA, fam.size*nh.size), nrow=fam.size, ncol=nh.size)
  
colnames(graph.mat) = c("self", "mom", "dad", paste0("sis",1:m), paste0("bro",1:m), paste0("dau", 1:m2), paste0("son", 1:m2))
rownames(graph.mat) = c("self", "mom", "dad", "mgm", "mgf", "pgm", "pgf", paste0("sis", 1:(m-1)), paste0("bro", 1:m), paste0("dau", 1:m2), paste0("son", 1:m2), paste0("maunt", 1:(m-1)), paste0("muncle", 1:m), paste0("paunt", 1:m), paste0("puncle", 1:(m-1)), "null")
graph.mat[, 1] = ped.struct$index
graph.mat[, 2] = ped.struct$index.mother
graph.mat[, 3] = ped.struct$index.father
  
for (i in 1:(fam.size-1)) {
  # m sisters
  sis = ped.struct$index[which(ped.struct$index.mother %in% graph.mat[i, 2] & graph.mat[i, 2]!=fam.size & ped.struct$index != i & ped.struct$gender==0)]
  null.sis = rep(ped.struct$index[fam.size], m-length(sis))
  graph.mat[i, 4:(4+m-1)] = c(sis[order(sis)], null.sis)
  # m brothers
  bro = ped.struct$index[which(ped.struct$index.mother %in% graph.mat[i, 2] & graph.mat[i, 2]!=fam.size & ped.struct$index != i & ped.struct$gender==1)]
  null.bro = rep(ped.struct$index[fam.size], m-length(bro))
  graph.mat[i, (4+m):(4+2*m-1)] = c(bro[order(bro)], null.bro)
  # m daughters
  dau = ped.struct$index[which( (ped.struct$index.mother == i | ped.struct$index.father == i) & ped.struct$gender==0 )]
  null.dau = rep(ped.struct$index[fam.size], pmax(0, m2-length(dau)))
  graph.mat[i, (4+2*m):(4+2*m+m2-1)] = c(dau[order(dau)], null.dau)[1:m2]
  # m sons
  son = ped.struct$index[which( (ped.struct$index.mother == i | ped.struct$index.father == i) & ped.struct$gender==1 )]
  null.son = rep(ped.struct$index[fam.size], pmax(0, m2-length(son)))
  graph.mat[i, (4+2*m+m2):(4+2*m+2*m2-1)] = c(son[order(son)], null.son)[1:m2]
}
  
graph.mat[fam.size,] = fam.size





graph.mat.py = graph.mat - 1
write.csv(graph.mat.py, "../nn_input/graph_mat.csv", row.names=F)
