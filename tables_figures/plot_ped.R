# generate figure of example pedigree

library(BayesMendel)
library(gridExtra)

data(brca.fam)
fam = brca.fam
fam$Twins = 0
fam$MotherID[which(fam$ID %in% c(13, 5, 25))] = 10
fam$FatherID[which(fam$ID %in% c(13, 5, 25))] = 11
fam = fam[which(fam$ID %in% c(1, 2, 3, 8, 9, 10, 11, 12, 13, 14)), ]
fam[which(fam$ID==14), c("AgeBreast", "AgeOvary")] = 55
myfamily <- new("BayesMendel", family=fam, counselee.id=1)
#plot.BayesMendel(myfamily, cex=0.2)

fam.copy = fam
fam.copy$Death = c(0, 0, 0, 0, 1, 0, 1, 1, 0, 0)
fam.copy$Death = 0
myfamily.death <- new("BayesMendel", family=fam.copy, counselee.id=1)

fam2 = fam
fam2$AffectedBreast = 0
fam2$AffectedOvary = 0
myfamily2 <- new("BayesMendel", family=fam2, counselee.id=1)

ids = expression(paste("H"^"0"))
ped = pedigree(fam$ID, fam$FatherID, fam$MotherID, ifelse(fam$Gender==0, 2, 1), rep(0, nrow(fam)))
#paste0("H", "^", 0:nrow(fam))

#plot.pedigree(ped, id = rep("", length(ped$id)))


plot2 = function(x, plot.status=F) {
  canc = "breast"
  family <- slot(x, "family")
  posterior <- slot(x, "posterior")
  probs <- slot(x, "probs")
  counselee.id <- slot(x, "counselee.id")
  predictions <- slot(x, "predictions")
  if (is.null(family$Death)) {
    status = rep(0, nrow(family))
  }
  if (!is.null(family$Death)) {
    status = family$Death
  }
  else {
    if (!is.null(family$AffectedBreast)) {
      canc <- "breast"
    }
  }
  id <- family$ID
  dadid <- family$FatherID
  momid <- family$MotherID
  sex <- family$Gender
  sex[sex == 0] <- 2
  refs <- 1:length(id)
  get1 <- refs[dadid == 0 & momid > 0]
  get2 <- refs[dadid > 0 & momid == 0]
  rem <- sort(c(get1, get2))
  num <- length(unique(momid[get1])) + length(unique(dadid[get2]))
  if (num > 0) {
    newid <- max(id) + (1:num)
    newid1 <- newid[1:length(unique(momid[get1]))]
    newid2 <- newid[(length(unique(momid[get1])) + 1):(length(newid))]
    id <- c(id, newid)
    sex <- c(sex, rep(1, length(unique(momid[get1]))), rep(2, 
                                                           length(unique(dadid[get2]))))
    dadid <- c(dadid, rep(0, num))
    momid <- c(momid, rep(0, num))
    if (length(get1) > 0) {
      for (i in 1:length(get1)) {
        dadid[id == id[get1[i]]] <- newid1[unique(momid[get1]) == 
                                             momid[get1[i]]]
      }
    }
    if (length(get2) > 0) {
      for (i in 1:length(get2)) {
        momid[id == id[get2[i]]] <- newid2[unique(dadid[get2]) == 
                                             dadid[get2[i]]]
      }
    }
  }
  if (canc == "breast") {
    bc2 <- family$AffectedBreast
    bc2[bc2 <= 1] <- 1
    affected <- cbind(family$AffectedBreast, family$AffectedOvary)
    affected <- affected + 1
    affected[affected == 3] <- 2
    ages.bc1 <- c(family$AgeBreast, rep(1, num))
    ages.oc <- c(family$AgeOvary, rep(1, num))
    ages.bc2 <- c(family$AgeBreastContralateral, rep(1, num))
    bc <- family$AffectedBreast
    oc <- family$AffectedOvary
  }
  if (canc %in% c("breast", "colon")) 
    affected <- rbind(affected, matrix(1, nrow = num, ncol = 2))
  if (canc %in% c("panc", "mela")) 
    affected <- c(affected, rep(1, num))
  status <- c(status, rep(0, num))
  relations <- NULL
  if (!is.null(family$Twins) & sum(family$Twins) > 0) {
    twin.ids <- unique(family$Twins[family$Twins > 0])
    id1 <- id2 <- code <- NULL
    for (ttt in 1:length(twin.ids)) {
      id1 <- c(id1, family$ID[family$Twins == unique(twin.ids)[ttt]][1])
      id2 <- c(id2, family$ID[family$Twins == unique(twin.ids)[ttt]][2])
      code <- c(code, 1)
    }
    relations <- matrix(c(id1, id2, code), byrow = F, ncol = 3)
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, 
                                status, relation = relations)
  }
  if ((!is.null(family$Twins) & sum(family$Twins) == 0) | is.null(family$Twins)) {
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, 
                                status)
  }
  par(xpd = TRUE)
  plt <- kinship2::plot.pedigree(ped.y, id = rep("", length(ped.y$id)))
  xp <- plt$x
  yp <- plt$y
  inc <- (abs(min(yp)) - abs(max(yp)))/12
  xp.inc <- (abs(min(xp)) - abs(max(xp)))/50
  couns <- counselee.id
  xp.couns <- xp[id == couns]
  yp.couns <- yp[id == couns]
  arrows(xp.couns - 0.35, yp.couns - 0.15, xp.couns - 0.1, 
         yp.couns - 0.1, length = 0.1, col="blue")
  if (canc == "breast") {
    xph <- xp[bc == 0 & oc == 0]
    yph <- yp[bc == 0 & oc == 0]
    ageh <- ages.bc1[bc == 0 & oc == 0]
    xph <- xph[ageh > 1]
    yph <- yph[ageh > 1]
    ageh <- ageh[ageh > 1]
    if (length(ageh) > 0) {
      text(x = xph, y = yph - inc, as.character(ageh), 
           cex = 0.75, adj=c(0.5, 2.2))
    }
    if (plot.status) {
      status.char = rep("alive", length(family$Death))
      status.char[which(family$Death==1)] = "dead"
      text(x = xp, y = yp - 2*inc, as.character(status.char), 
           cex = 0.75, adj=c(0.5, 2.2))
    }
    xph <- xp[bc >= 1]
    yph <- yp[bc >= 1]
    ageh <- ages.bc1[bc >= 1]
    ageh <- ageh[ageh > 0]
    xph <- xph[ageh > 1]
    yph <- yph[ageh > 1]
    ageh <- ageh[ageh > 1]
    if (length(ageh) > 0) {
      text(x = xph, y = yph - inc, paste("BC ", as.character(ageh)), 
           cex = 0.75, col = "red", font = 2, adj=c(0.5, 2.2))
    }
    xph <- xp[bc == 2]
    yph <- yp[bc == 2]
    ageh <- ages.bc2[bc == 2]
    ageh <- ageh[ageh > 0]
    xph <- xph[ageh > 1]
    yph <- yph[ageh > 1]
    ageh <- ageh[ageh > 1]
    if (length(ageh) > 0) {
      text(x = xph, y = yph - 1.5 * inc, paste("BC (2) ", 
                                               as.character(ageh)), cex = 0.75, col = "red", 
           font = 2)
    }
    xph <- xp[oc > 0 & bc == 0]
    yph <- yp[oc > 0 & bc == 0]
    ageh <- ages.oc[oc > 0 & bc == 0]
    xph <- xph[ageh > 1]
    ageh <- ageh[ageh > 0]
    yph <- yph[ageh > 1]
    ageh <- ageh[ageh > 1]
    if (length(ageh) > 0) {
      text(x = xph, y = yph - inc, paste("OC ", as.character(ageh)), 
           cex = 0.75, col = "red", font = 2, adj=c(0.5, 2.2))
    }
    xph <- xp[oc > 0 & bc == 1]
    yph <- yp[oc > 0 & bc == 1]
    ageh <- ages.oc[oc > 0 & bc == 1]
    xph <- xph[ageh > 1]
    ageh <- ageh[ageh > 0]
    yph <- yph[ageh > 1]
    ageh <- ageh[ageh > 1]
    if (length(ageh) > 0) {
      text(x = xph, y = yph - 1.5 * inc, paste("OC ", as.character(ageh)), 
           cex = 0.75, col = "red", font = 2, adj=c(0.5, 2.2))
    }
    xph <- xp[oc > 0 & bc == 2]
    yph <- yp[oc > 0 & bc == 2]
    ageh <- ages.oc[oc > 0 & bc == 2]
    xph <- xph[ageh > 1]
    ageh <- ageh[ageh > 0]
    yph <- yph[ageh > 1]
    ageh <- ageh[ageh > 1]
    if (length(ageh) > 0) {
      text(x = xph, y = yph - 2 * inc, paste("OC ", as.character(ageh)), 
           cex = 0.75, col = "red", font = 2, adj=c(0.5, 2.2))
    }
  }
}

plot3 = function(x) {
  canc = "breast"
  family <- slot(x, "family")
  posterior <- slot(x, "posterior")
  probs <- slot(x, "probs")
  counselee.id <- slot(x, "counselee.id")
  predictions <- slot(x, "predictions")
  if (is.null(family$Death)) {
    status = rep(0, nrow(family))
  }
  if (!is.null(family$Death)) {
    status = family$Death
  }
  else {
    if (!is.null(family$AffectedBreast)) {
      canc <- "breast"
    }
  }
  id <- family$ID
  dadid <- family$FatherID
  momid <- family$MotherID
  sex <- family$Gender
  sex[sex == 0] <- 2
  refs <- 1:length(id)
  get1 <- refs[dadid == 0 & momid > 0]
  get2 <- refs[dadid > 0 & momid == 0]
  rem <- sort(c(get1, get2))
  num <- length(unique(momid[get1])) + length(unique(dadid[get2]))
  if (num > 0) {
    newid <- max(id) + (1:num)
    newid1 <- newid[1:length(unique(momid[get1]))]
    newid2 <- newid[(length(unique(momid[get1])) + 1):(length(newid))]
    id <- c(id, newid)
    sex <- c(sex, rep(1, length(unique(momid[get1]))), rep(2, 
                                                           length(unique(dadid[get2]))))
    dadid <- c(dadid, rep(0, num))
    momid <- c(momid, rep(0, num))
    if (length(get1) > 0) {
      for (i in 1:length(get1)) {
        dadid[id == id[get1[i]]] <- newid1[unique(momid[get1]) == 
                                             momid[get1[i]]]
      }
    }
    if (length(get2) > 0) {
      for (i in 1:length(get2)) {
        momid[id == id[get2[i]]] <- newid2[unique(dadid[get2]) == 
                                             dadid[get2[i]]]
      }
    }
  }
  if (canc == "breast") {
    bc2 <- family$AffectedBreast
    bc2[bc2 <= 1] <- 1
    affected <- cbind(family$AffectedBreast, family$AffectedOvary)
    affected <- affected + 1
    affected[affected == 3] <- 2
    ages.bc1 <- c(family$AgeBreast, rep(1, num))
    ages.oc <- c(family$AgeOvary, rep(1, num))
    ages.bc2 <- c(family$AgeBreastContralateral, rep(1, num))
    bc <- family$AffectedBreast
    oc <- family$AffectedOvary
  }
  if (canc %in% c("breast", "colon")) 
    affected <- rbind(affected, matrix(1, nrow = num, ncol = 2))
  if (canc %in% c("panc", "mela")) 
    affected <- c(affected, rep(1, num))
  status <- c(status, rep(0, num))
  relations <- NULL
  if (!is.null(family$Twins) & sum(family$Twins) > 0) {
    twin.ids <- unique(family$Twins[family$Twins > 0])
    id1 <- id2 <- code <- NULL
    for (ttt in 1:length(twin.ids)) {
      id1 <- c(id1, family$ID[family$Twins == unique(twin.ids)[ttt]][1])
      id2 <- c(id2, family$ID[family$Twins == unique(twin.ids)[ttt]][2])
      code <- c(code, 1)
    }
    relations <- matrix(c(id1, id2, code), byrow = F, ncol = 3)
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, 
                                status, relation = relations)
  }
  if ((!is.null(family$Twins) & sum(family$Twins) == 0) | is.null(family$Twins)) {
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, 
                                status)
  }
  par(xpd = TRUE)
  plt <- kinship2::plot.pedigree(ped.y, id = rep("", length(ped.y$id)))
  xp <- plt$x
  yp <- plt$y
  inc <- (abs(min(yp)) - abs(max(yp)))/12
  xp.inc <- (abs(min(xp)) - abs(max(xp)))/50
  couns <- counselee.id
  xp.couns <- xp[id == couns]
  yp.couns <- yp[id == couns]
  arrows(xp.couns - 0.35, yp.couns - 0.15, xp.couns - 0.1, 
         yp.couns - 0.1, length = 0.1, col="blue")
  if (canc == "breast") {
    ##text(x = xp, y = yp - inc, paste0("H", 0:(length(xp)-1)), cex = 0.75, font = 2, adj=c(0.5, 2.2))
    #text(x = xp, y = yp - inc, parse(text=paste("H^", 0:(length(xp)-1))), cex = 0.75, font = 2, adj=c(0.5, 2.2))
    text(x = xp, y = yp - inc, parse(text=paste("italic(a[", 0:(length(xp)-1), "]^{l-1})")), 
         cex = 0.9, font = 2, adj=c(0.5, 1.8))
    
    #expression(paste("H"^"0"))
    
  }
}

fam$agecur = c(50, 70, 75, 90, 92, 88, 90, 74, 72, 73)
fam$AgeBreast[which(fam$AffectedBreast==0)] = fam$agecur[which(fam$AffectedBreast==0)]
fam$AgeOvary[which(fam$AffectedOvary==0)] = fam$agecur[which(fam$AffectedOvary==0)]
myfamily <- new("BayesMendel", family=fam, counselee.id=1)

png('ped_vital_status.png', width=700, height=650, res=120, pointsize=15)
plot2(myfamily.death, T)
legend(legend=c("Breast Cancer (BC)", "Ovarian Cancer (OC)"), fill=c("black", "black"), 
       density=c(NA, 35), bty="n", border=c("black", "black"), cex=0.9, x=2.4, y=2.8)
dev.off()


