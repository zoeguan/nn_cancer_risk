library(dplyr)
library(data.table)

# reference pedigree structure
rel.types.ref = data.frame(relationship=c("Self", "Mother", "Father", "Sister", "Brother", "Daughter", "Son", "Maternal Grandmother", "Maternal Grandfather", "Paternal Grandmother", "Paternal Grandfather", "Maternal Aunt", "Maternal Uncle", "Paternal Aunt", "Paternal Uncle"), 
                           num=c(1, 1, 1, 2, 3, 3, 3, 1, 1, 1, 1, 2, 3, 3, 2),
                           gender=c(0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
                           suffix=c("self", "mom", "dad", "sis", "bro", "dau", "son", "mgm", "mgf", "pgm", "pgf", "maunt", "muncle", "paunt", "puncle"))
rel.types.ref$relationship = as.character(rel.types.ref$relationship)

### function that standardizes and flattens pedigrees
# ped: dataframe of pedigrees in BayesMendel format 
# rel.types: dataframe with same format as rel.types.ref
# fixed.rels: vector relative types that appear exactly once in each pedigree (besides "Self")
standardize_ped = function(ped, rel.types = rel.types.ref, fixed.rels=c("Mother", "Father", "Maternal Grandmother", "Maternal Grandfather", "Paternal Grandmother", "Paternal Grandfather")) {
  
  # shorten names of cancer and cancer age variables 
  ped[, c("AffB", "AgeB", "AffO", "AgeO", "Age")] = ped[, c("AffectedBreast", "AgeBreast", "AffectedOvary", "AgeOvary", "agecur")]
  ped$AgeB[which(ped$AffB==0 | ped$AgeB %in% NA)] = 0
  ped$AffO[which(ped$AffO %in% NA)] = 0
  ped$AgeO[which(ped$AffO==0 | ped$AgeO %in% NA)] = 0
  # add missing indicator
  ped$Miss = 0
  
  features = c("Gender", "AffB", "AgeB", "AffO", "AgeO", "Age", "Miss")
  
  ped = ped[, c("FamID", "ID", features, "relationship")]
  
  ped.std = ped[which(ped$relationship=="Self"), c("FamID", features)]
  
  # get features of fixed relatives
  for (i in 1:length(fixed.rels)) {
    rel.i = ped[which(ped$relationship %in% fixed.rels[i]), ]
    suffix.i = rel.types$suffix[which(rel.types$relationship==fixed.rels[i])]
    rel.i = rel.i[, c("FamID", features)]
    names(rel.i)[which(names(rel.i) %in% features)] = paste0(features, "_", suffix.i)
    ped.std = left_join(ped.std, rel.i, by="FamID")
  }
  
  # get features of other relative types
  other.rels = setdiff(rel.types$relationship, c(fixed.rels, "Self"))
  for (i in 1:length(other.rels)) {
    rel.i = ped[which(ped$relationship %in% other.rels[i]), ]
    num.i = rel.types$num[which(rel.types$relationship==other.rels[i])]
    suffix.i = rel.types$suffix[which(rel.types$relationship==other.rels[i])]
    # count the number of relatives of type i
    rel.counts = data.frame(table(rel.i$FamID))
    rel.counts$FamID = as.numeric(as.character(rel.counts$Var1))
    ped.std = left_join(ped.std, rel.counts[, c("FamID", "Freq")], by="FamID")
    ped.std$Freq[which(ped.std$Freq %in% NA)] = 0
    # get the number of relatives of type i that need to be added to the pedigree
    ped.std$num.add = pmax(0, num.i - ped.std$Freq)
    # "null" relatives who will be added to the pedigree
    null.vals = data.frame(FamID = rep(ped.std$FamID, ped.std$num.add),
                           Gender = rel.i$Gender[1], AffB=0, AgeB=0, AffO=0, AgeO=0, Age=0, Miss=1)
    # actual relatives of type i from families that don't require pruning wrt i
    subset.1 = rbind(rel.i[which(rel.i$FamID %in% ped.std$FamID[which(ped.std$Freq <= num.i)]), names(null.vals)], null.vals)
    # actual relatives of type i from families that require pruning wrt i
    subset.2 = rel.i[which(rel.i$FamID %in% ped.std$FamID[which(ped.std$Freq > num.i)]), names(null.vals)]
    # sample the required number of relatives
    subset.2 = data.frame(subset.2 %>% group_by(FamID) %>% sample_n(size=num.i))
    
    # combine inputs for relatives of type i
    rel.i = rbind(subset.1, subset.2)
    rel.i = rel.i[order(rel.i$FamID), ]
    rel.i$num = rep(paste0(suffix.i, 1:num.i), length(unique(rel.i$FamID)))
    # flatten inputs
    rel.i.wide = data.frame(dcast(data.table(rel.i), FamID ~ num, value.var = features, fun.aggregate=function(x) x[1]))
    rel.i.wide = rel.i.wide[ , c(1, sapply(1:num.i, function(x) which(grepl(x, names(rel.i.wide))))), ]

    ped.std[, c("Freq", "num.add")] = NULL
    ped.std = left_join(ped.std, rel.i.wide, by="FamID")
    
  }
  
  return(ped.std)

}



