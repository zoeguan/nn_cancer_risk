load(paste0("../families/simFamAJ_", 3, ".RData"))

# 5 risk scenarios for a fixed family structure
# vary phenotypes of mother, maternal grandmother

# 1 # 0, 0
# 2 # 0, 80b
# 3 # 0 , 60b
# 4 # 50b, 60b
# 5 # 50b and 60o, 60b
# 6 # 50b and 60o, 40b

fam1 = fam[which(fam$FamID==56426), ]

load("../nn_input/nn_input.RData")
nn.input[which(nn.input$FamID %in% c(56426)), which(grepl("maunt", names(nn.input)))]

age.pro = 40

# remove maternal aunts to increase strength of signal from affected mother and maternal grandmother
fam1 = fam1[which(!fam1$ID %in% c(40, 41, 43)), c(1:9, 36:40)]

fam1$agecur[which(fam1$ID==6)] = 70
fam1$agecur[which(fam1$ID==50)] = age.pro
fam1$AffectedBreast = 0
fam1$AgeBreast = fam1$agecur
fam1$AffectedOvary = 0
fam1$AgeOvary = fam1$agecur

### NN inputs
nn.in = nn.input[which(nn.input$FamID %in% c(56426)), ]
nn.in[, which(grepl("AffB", names(nn.input)))] = 0
nn.in[, which(grepl("AffO", names(nn.input)))] = 0
nn.in[, which(grepl("AgeB", names(nn.input)))] = 0
nn.in[, which(grepl("AgeO", names(nn.input)))] = 0
nn.in[, which(grepl("maunt2", names(nn.input)))] = 0
nn.in[, which(grepl("Miss_maunt2", names(nn.input)))] = 1
nn.in$Age_mom = 70
nn.in$Age = age.pro

in1 = nn.input[1:6, ]
for (i in 1:6) {
  in1[i, ] = nn.in
}
in1$AffB_mom = c(0, 0, 0, 1, 1, 1)
in1$AgeB_mom = c(0, 0, 0, 50, 50, 50)
in1$AffO_mom = c(0, 0, 0, 0, 1, 1)
in1$AgeO_mom = c(0, 0, 0, 0, 60, 60)
in1$AffB_mgm = c(0, 1, 1, 1, 1, 1)
in1$AgeB_mgm = c(0, 80, 60, 60, 60, 40)


### BRCAPRO inputs

library(BayesMendel)

fam1$FatherID = as.numeric(fam1$FatherID)
fam1$MotherID = as.numeric(fam1$MotherID)
fam1$AgeBreast[which(fam1$ID==50)] = age.pro
fam1$AgeOvary[which(fam1$ID==50)] = age.pro

bparams = brcaparams(age.by=5, age.to=94)
bparams$comprisk = 0*bparams$comprisk + 1

# version 1
out1 = brcapro(fam1, counselee.id=50, params=bparams)
# version 2
fam2 = fam1
fam2$AffectedBreast[which(fam1$ID==4)] = 1
fam2$AgeBreast[which(fam1$ID==4)] = 80
out2 = brcapro(fam2, counselee.id=50, params=bparams)
# version 3 
fam3 = fam2
fam3$AgeBreast[which(fam1$ID==4)] = 60
out3 = brcapro(fam3, counselee.id=50, params=bparams)
# version 4
fam4 = fam3
fam4$AffectedBreast[which(fam1$ID==6)] = 1
fam4$AgeBreast[which(fam1$ID==6)] = 50
out4 = brcapro(fam4, counselee.id=50, params=bparams)
# version 5
fam5 = fam4
fam5$AffectedOvary[which(fam1$ID==4)] = 1
fam5$AgeOvary[which(fam1$ID==4)] = 60
out5 = brcapro(fam5, counselee.id=50, params=bparams)
# version 6
fam6 = fam5
fam6$AgeBreast[which(fam1$ID==4)] = 40
out6 = brcapro(fam6, counselee.id=50, params=bparams)

in1$brcapro10 = c(out1@predictions$`Breast Ca Risk`[2],
                  out2@predictions$`Breast Ca Risk`[2],
                  out3@predictions$`Breast Ca Risk`[2],
                  out4@predictions$`Breast Ca Risk`[2],
                  out5@predictions$`Breast Ca Risk`[2],
                  out6@predictions$`Breast Ca Risk`[2])

write.csv(in1, file="../nn_input/example_family.csv")


cnn1 = in1[, which(!grepl("FDR|SDR|brcapro|BC.", names(in1)))]
cnn1[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
cnn1$Miss_null = 1
cnn1[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = in1[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]

write.csv(cnn1, file="../nn_input/example_family_cnn.csv")
