
prop.missing = c(0.05, 0.1, 0.2, 0.4)

for (suffix in paste0("missingage", prop.missing)) {
  load(paste0("../nn_input/nn_input_", suffix, ".RData"))
  
  cnn.input = nn.input.i[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input.i)))]
  cnn.input[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
  cnn.input$Miss_null = 1
  cnn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input.i[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]
  write.csv(cnn.input, file=paste0("../nn_input/cnn_input_", suffix, ".csv"))
}


for (suffix in paste0("missing", prop.missing)) {
  load(paste0("../nn_input/nn_input_", suffix, ".RData"))

  
  cnn.input = nn.input.i[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input.i)))]
  cnn.input[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
  cnn.input$Miss_null = 1
  cnn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input.i[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]
  write.csv(cnn.input, file=paste0("../nn_input/cnn_input_", suffix, ".csv"))
}
