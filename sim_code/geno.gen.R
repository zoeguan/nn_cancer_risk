library(dplyr)

### generate BRCA1 and BRCA2 genotypes
geno.gen = function(fam, pbrca1=0.01366243, pbrca2=0.01168798, max.depth=5) {

  fam$BRCA1.a1 = NA
  fam$BRCA1.a2 = NA
  fam$BRCA2.a1 = NA
  fam$BRCA2.a2 = NA
  
  # alleles of founders and people with unknown MotherID/FatherID
  ind.unk.father = which(fam$FatherID==0)
  fam$BRCA1.a1[ind.unk.father] = rbinom(length(ind.unk.father), 1, pbrca1)
  fam$BRCA2.a1[ind.unk.father] = rbinom(length(ind.unk.father), 1, pbrca2)
  ind.unk.mother = which(fam$MotherID==0)
  fam$BRCA1.a2[ind.unk.mother] = rbinom(length(ind.unk.mother), 1, pbrca1)
  fam$BRCA2.a2[ind.unk.mother] = rbinom(length(ind.unk.mother), 1, pbrca2)
  
  # add alleles of mother and father to each row 
  fam = left_join(fam, fam[, c("unique.id", "BRCA1.a1", "BRCA1.a2", "BRCA2.a1", "BRCA2.a2")], by=c("unique.fid" = "unique.id"), suffix=c("", ".f"))
  fam = left_join(fam, fam[, c("unique.id", "BRCA1.a1", "BRCA1.a2", "BRCA2.a1", "BRCA2.a2")], by=c("unique.mid" = "unique.id"), suffix=c("", ".m"))
  
  # sample BRCA1/2 allele from father
  ind.update.brca.f = which(fam$BRCA1.a1 %in% NA & fam$BRCA1.a1.f>=0)
  fam$BRCA1.a1[ind.update.brca.f] = ifelse(sample(0:1, length(ind.update.brca.f), replace=TRUE), fam$BRCA1.a1.f[ind.update.brca.f], fam$BRCA1.a2.f[ind.update.brca.f])
  fam$BRCA2.a1[ind.update.brca.f] = ifelse(sample(0:1, length(ind.update.brca.f), replace=TRUE), fam$BRCA2.a1.f[ind.update.brca.f], fam$BRCA2.a2.f[ind.update.brca.f])
  # sample BRCA1/2 allele from mother 
  ind.update.brca.m = which(fam$BRCA1.a2 %in% NA & fam$BRCA1.a1.m>=0)
  fam$BRCA1.a2[ind.update.brca.m] = ifelse(sample(0:1, length(ind.update.brca.m), replace=TRUE), fam$BRCA1.a1.m[ind.update.brca.m], fam$BRCA1.a2.m[ind.update.brca.m])
  fam$BRCA2.a2[ind.update.brca.m] = ifelse(sample(0:1, length(ind.update.brca.m), replace=TRUE), fam$BRCA2.a1.m[ind.update.brca.m], fam$BRCA2.a2.m[ind.update.brca.m])
  
  # repeat until everyone's alleles have been generated
  depth = 1
  while(length(which(fam$BRCA1.a1 %in% NA | fam$BRCA1.a2 %in% NA))>0 & depth<=max.depth) {
    depth = depth + 1
    fam[, paste0(c("BRCA1.a1", "BRCA1.a2", "BRCA2.a1", "BRCA2.a2"), ".f")] = NULL
    fam[, paste0(c("BRCA1.a1", "BRCA1.a2", "BRCA2.a1", "BRCA2.a2"), ".m")] = NULL
    fam = left_join(fam, fam[, c("unique.id", "BRCA1.a1", "BRCA1.a2", "BRCA2.a1", "BRCA2.a2")], by=c("unique.fid" = "unique.id"), suffix=c("", ".f"))
    fam = left_join(fam, fam[, c("unique.id", "BRCA1.a1", "BRCA1.a2", "BRCA2.a1", "BRCA2.a2")], by=c("unique.mid" = "unique.id"), suffix=c("", ".m"))
    
    # sample BRCA1/2 allele from father
    ind.update.brca.f = which(fam$BRCA1.a1 %in% NA & fam$BRCA1.a1.f>=0)
    fam$BRCA1.a1[ind.update.brca.f] = ifelse(sample(0:1, length(ind.update.brca.f), replace=TRUE), fam$BRCA1.a1.f[ind.update.brca.f], fam$BRCA1.a2.f[ind.update.brca.f])
    fam$BRCA2.a1[ind.update.brca.f] = ifelse(sample(0:1, length(ind.update.brca.f), replace=TRUE), fam$BRCA2.a1.f[ind.update.brca.f], fam$BRCA2.a2.f[ind.update.brca.f])
    # sample BRCA1/2 allele from mother 
    ind.update.brca.m = which(fam$BRCA1.a2 %in% NA & fam$BRCA1.a1.m>=0)
    fam$BRCA1.a2[ind.update.brca.m] = ifelse(sample(0:1, length(ind.update.brca.m), replace=TRUE), fam$BRCA1.a1.m[ind.update.brca.m], fam$BRCA1.a2.m[ind.update.brca.m])
    fam$BRCA2.a2[ind.update.brca.m] = ifelse(sample(0:1, length(ind.update.brca.m), replace=TRUE), fam$BRCA2.a1.m[ind.update.brca.m], fam$BRCA2.a2.m[ind.update.brca.m])
  }
  
  # set BRCA1/2 carrier status to 1 if either allele is mutated
  fam$BRCA1 = fam$BRCA1.a1 | fam$BRCA1.a2
  fam$BRCA2 = fam$BRCA2.a1 | fam$BRCA2.a2
  
  return(fam)
  
}
