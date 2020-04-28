library(xtable)

### performance table for simulations under misreporting
# results from train_lr_noise.py, train_fcnn_noise.py, train_cnn_noise.py, train_fcnn_noise.ipynb

digits = c(1, 2, 2, 3)

res.fcnn = read.csv("../results/res_noise800000_50_0.csv")
res.fcnn = res.fcnn[, c(1, 3, 2, 4)]
res.fcnn = sapply(1:4, function(x) round(res.fcnn[, x], digits[x]))
res.cnn = read.csv("../results/res_cnn_noise800000_15_0.csv")
res.cnn = res.cnn[, c(1, 3, 2, 4)]
res.cnn = sapply(1:4, function(x) round(res.cnn[, x], digits[x]))
res.brca = read.csv("../results/res_brca_noise.csv")
res.brca = res.brca[, c(1, 4, 2, 3)]
res.brca = sapply(1:4, function(x) round(res.brca[, x], digits[x]))
res.lr = read.csv("../results/res_lr_noise.csv")
res.lr = res.lr[, c(1, 4, 2, 3)]
res.lr = sapply(1:4, function(x) round(res.lr[, x], digits[x]))

perf.table = data.frame(OE=rep(NA, 4), AUC=NA, BS=NA)
perf.table[1, ] = apply(res.fcnn[, 2:4], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
perf.table[2, ] = apply(res.cnn[, 2:4], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
perf.table[3, ] = apply(res.brca[, 2:4], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))

perf.table[4, ] = apply(res.lr[, 2:4], 2, function(x) paste0(x[1], " (", x[2], ", ", x[3], ")"))
rownames(perf.table) = c("FCNN", "CNN", "BRCAPRO", "LR")

library(xtable)
xtable(perf.table)




oe = t(read.csv("../results/oe_boot_noise2.csv"))[-1, ]
auc = t(read.csv("../results/auc_boot_noise2.csv"))[-1, ]
brier = t(read.csv("../results/bs_boot_noise2.csv"))[-1, ]



boot.table = data.frame(comparison=c("FCNN>CNN", "FCNN>BRCAPRO", "FCNN>LR", "CNN>BRCAPRO", "CNN>LR"), OE=NA, AUC=NA, BS=NA)

boot.table[, "AUC"] = c(sapply(2:4, function(x) length(which(auc[,1] > auc[,x]))), sapply(3:4, function(x) length(which(auc[,2] > auc[,x]))))
boot.table[, "BS"] = c(sapply(2:4, function(x) length(which(brier[,1] < brier[,x]))), sapply(3:4, function(x) length(which(brier[,2] < brier[,x]))))
boot.table[, "OE"] = c(sapply(2:4, function(x) length(which( abs(oe[,1]-1) < abs(oe[,x]-1) ))), sapply(3:4, function(x) length(which(abs(oe[,2]-1) < abs(oe[,x]-1) ))))

boot.table[, 2:4] = (boot.table[, 2:4])/nrow(oe)

print(xtable(boot.table, digits=c(0, 0, 3, 3, 3)), include.rownames=F)
