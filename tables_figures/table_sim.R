library(xtable)

### performance table for simulations
# results from pred_fam.ipynb

digits = c(1, 2, 2, 3, 2)

res.fcnn = read.csv("../results/res_800000_30_0.csv")
res.fcnn = res.fcnn[, c(1, 5, 2, 3, 4)]
res.fcnn = sapply(1:5, function(x) round(res.fcnn[, x], digits[x]))
res.cnn = read.csv("../results/res_cnn_800000_15.csv")
res.cnn = res.cnn[, c(1, 5, 2, 3, 4)]
res.cnn = sapply(1:5, function(x) round(res.cnn[, x], digits[x]))
res.brca = read.csv("../results/res_brca.csv")
res.brca$corr = 1
res.brca = res.brca[, c(1, 4, 2, 3, 5)]
res.brca = sapply(1:5, function(x) round(res.brca[, x], digits[x]))
res.lr = read.csv("../results/res_lr.csv")
res.lr = res.lr[, c(1, 5, 2, 3, 4)]
res.lr = sapply(1:5, function(x) round(res.lr[, x], digits[x]))

perf.table = data.frame(OE=rep(NA, 4), AUC=NA, BS=NA, corr=NA)
perf.table[1, ] = apply(res.fcnn[, 2:5], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
perf.table[2, ] = apply(res.cnn[, 2:5], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
perf.table[3, ] = apply(res.brca[, 2:5], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))

perf.table[4, ] = apply(res.lr[, 2:5], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
rownames(perf.table) = c("FCNN", "CNN", "BRCAPRO", "LR")

library(xtable)
xtable(perf.table)

### comparisons across bootstrap replicates

oe = t(read.csv("../results/oe_boot.csv"))[-1, ]
auc = t(read.csv("../results/auc_boot.csv"))[-1, ]
brier = t(read.csv("../results/bs_boot.csv"))[-1, ]
corr = t(read.csv("../results/corr_boot.csv"))[-1, ]

boot.table = data.frame(comparison=c("FCNN>CNN", "FCNN>BRCAPRO", "FCNN>LR", "CNN>BRCAPRO", "CNN>LR"), OE=NA, AUC=NA, BS=NA, cor=NA)
boot.table[, "AUC"] = c(sapply(2:4, function(x) length(which(auc[,1] > auc[,x]))), sapply(3:4, function(x) length(which(auc[,2] > auc[,x]))))
boot.table[, "BS"] = c(sapply(2:4, function(x) length(which(brier[,1] < brier[,x]))), sapply(3:4, function(x) length(which(brier[,2] < brier[,x]))))
boot.table[, "OE"] = c(sapply(2:4, function(x) length(which( abs(oe[,1]-1) < abs(oe[,x]-1) ))), sapply(3:4, function(x) length(which(abs(oe[,2]-1) < abs(oe[,x]-1) ))))
boot.table[, "cor"] = c(sapply(2:4, function(x) length(which(corr[,1] > corr[,x]))), sapply(3:4, function(x) length(which(corr[,2] > corr[,x]))))
boot.table[, 2:5] = boot.table[, 2:5]/nrow(oe)


print(xtable(boot.table, digits=c(0, 0, 3, 3, 3, 3)), include.rownames=F)


