
### missing values in test set

load("../nn_input/nn_input.RData")

rel.names = setdiff(gsub("AgeB|_", "" ,names(nn.input)[which(grepl("AgeB", names(nn.input)))]), "")

prop.missing = c(0.05, 0.1, 0.3)

# missing diagnosis ages
set.seed(1)
for (i in prop.missing) {
  print(i)
  nn.input.i = nn.input[1:87353, ]
  for (rel.type in rel.names) {
    ind.aff = which(nn.input[1:87353, paste0("AffB_", rel.type)]==1)
    ind.sample = sample(ind.aff, round(i*length(ind.aff)), replace=F)
    nn.input.i[ind.sample, paste0("AgeB_", rel.type)] = pmin(50, nn.input.i[ind.sample, paste0("Age_", rel.type)])
    
    ind.aff.ov = which(nn.input[1:87353, paste0("AffO_", rel.type)]==1)
    ind.sample = sample(ind.aff.ov, round(i*length(ind.aff.ov)), replace=F)
    nn.input.i[ind.sample, paste0("AgeO_", rel.type)] = pmin(50, nn.input.i[ind.sample, paste0("Age_", rel.type)])
  }
  save(nn.input.i, file=paste0("../nn_input/test_input_missingages", i, ".RData"))
  write.csv(nn.input.i, file=paste0("../nn_input/test_input_missingages", i, ".csv"))
  
  cnn.input = nn.input.i[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input.i)))]
  cnn.input[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
  cnn.input$Miss_null = 1
  cnn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input.i[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]
  write.csv(cnn.input, file=paste0("../nn_input/cnntest_input_missingages", i, ".csv"))
}


# missing relatives
set.seed(1)
for (i in prop.missing) {
  print(i)
  nn.input.i = nn.input[1:87353, ]
  for (rel.type in rel.names) {
    ind.aff = which(nn.input[1:87353, paste0("Miss_", rel.type)]==0)
    ind.sample = sample(ind.aff, round(i*length(ind.aff)), replace=F)
    nn.input.i[ind.sample, paste0(c("AffB_", "AgeB_", "AffO_", "AgeO_", "Age_"), rel.type)] = 0
    nn.input.i[ind.sample, paste0(c("Miss_"), rel.type)] = 1
  }
  save(nn.input.i, file=paste0("../nn_input/test_input_missingrels", i, ".RData"))
  write.csv(nn.input.i, file=paste0("../nn_input/test_input_missingrels", i, ".csv"))
  
  cnn.input = nn.input.i[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input.i)))]
  cnn.input[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
  cnn.input$Miss_null = 1
  cnn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input.i[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]
  write.csv(cnn.input, file=paste0("../nn_input/cnntest_input_missingrels", i, ".csv"))
}


# missing unaffected relatives
set.seed(1)
for (i in prop.missing) {
  print(i)
  nn.input.i = nn.input[1:87353, ]
  for (rel.type in rel.names) {
    ind.aff = which(nn.input[1:87353, paste0("AffB_", rel.type)]==0 & nn.input[1:87353, paste0("AffO_", rel.type)]==0)
    ind.sample = sample(ind.aff, round(i*length(ind.aff)), replace=F)
    nn.input.i[ind.sample, paste0(c("AffB_", "AgeB_", "AffO_", "AgeO_", "Age_"), rel.type)] = 0
    nn.input.i[ind.sample, paste0(c("Miss_"), rel.type)] = 1
  }
  save(nn.input.i, file=paste0("../nn_input/test_input_missingunaff", i, ".RData"))
  write.csv(nn.input.i, file=paste0("../nn_input/test_input_missingunaff", i, ".csv"))
  
  cnn.input = nn.input.i[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input.i)))]
  cnn.input[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
  cnn.input$Miss_null = 1
  cnn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input.i[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]
  write.csv(cnn.input, file=paste0("../nn_input/cnntest_input_missingunaff", i, ".csv"))
}