

### AUC and correlation plots
# results from train_fcnn.py, train_cnn.py

results = data.frame(size = c(25000, 50000, 100000, 200000, 400000, 800000))
results$log2.size = log2(results$size)
results[, c("auc", "auc_lo", "auc_hi", "cor", "cor_lo", "cor_hi")] = NA
results$Model = "FCNN"
results_c = results
results_c[, c("auc", "auc_lo", "auc_hi", "cor", "cor_lo", "cor_hi")] = NA
results_c$Model = "CNN"

options(scipen=999)

### FCNN results
epochs = c(70, 50, 30, 40, 30, 30)
seed = 0
for (i in 1:length(epochs)) {
  res.i = read.csv(paste0("../results/res_", results$size[i], "_", epochs[i], "_", seed, ".csv"))
  results[i, c("auc", "auc_lo", "auc_hi")] = res.i$auc
  results[i, c("cor", "cor_lo", "cor_hi")] = res.i$corr
}

### CNN results
epochs2 = c(80, 50, 30, 80, 15, 15)
seed2 = 0
for (i in 1:length(epochs)) {
  res.i = read.csv(paste0("../results/res_cnn", results$size[i], "_", epochs2[i], "_", seed2, ".csv"))
  results_c[i, c("auc", "auc_lo", "auc_hi")] = res.i$auc
  results_c[i, c("cor", "cor_lo", "cor_hi")] = res.i$corr
}

brca = read.csv("../results/res_brca.csv")

results_all = rbind(results, results_c)

library(ggplot2)

p1 = ggplot(results_all, aes(x=log2.size, y=auc, color=Model)) + 
  geom_errorbar(aes(ymin=auc_lo, ymax=auc_hi), width=0.16) + xlab(expression("Log"[2]*"(Training Set Size)")) + geom_point(aes(shape=Model), size=1.5, stroke=1) +
  ylab("AUC") + 
  geom_hline(yintercept = brca$auc[1], color="blue") +
  geom_text(aes(16.02, brca$auc[1], 
                label="BRCAPRO AUC", hjust=0.8, vjust=1.5), color="blue",
            size=5) +
  geom_hline(yintercept = brca$auc[2], color="blue", linetype=2) +
  geom_hline(yintercept = brca$auc[3], color="blue", linetype=2) + 
  theme_bw(base_size=20) + scale_shape(solid = FALSE) + theme(axis.title.x = element_blank(),
                                                              axis.text.x = element_blank())  #+ ylim(c(0.5, 0.7))


ggsave("../results/auc_plot2.png", scale=1, dpi=300)

p2 = ggplot(results_all, aes(log2.size, cor, color=Model), show.legend = FALSE) + 
  geom_errorbar(aes(ymin=cor_lo, ymax=cor_hi), width=0.16) +
  geom_point(aes(shape=Model), size=1.5, stroke=1) +
  xlab(expression("Log"[2]*"(Training Set Size)")) +
  ylab(expression(rho)) + 
  #geom_hline(yintercept = 1, color="blue") +
  theme_bw(base_size=20) + scale_shape(solid = FALSE) + theme(legend.text = element_text(color = "white"),
                                                              legend.title = element_text(color = "white"),
                                                              legend.key = element_rect(fill = "white")) + 
  scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) #+ ylim(c(0.4, 1))

ggsave("../results/corr_plot.png", scale=1, dpi=300)

library(cowplot)

p4 = cowplot::plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = c(0.5, 0.5), axis = "l")
ggsave("../results/auc_corr_plots.png", p4, scale=1.2, dpi=300)

#egg::ggarrange(p1, p2, heights = c(0.5, 0.5))



p1b = ggplot(results_all, aes(x=log2.size, y=auc, color=Model)) + 
  geom_errorbar(aes(ymin=auc_lo, ymax=auc_hi), width=0.16) + xlab(expression("Log"[2]*"(Training Set Size)")) + geom_point(aes(shape=Model), size=1.5, stroke=1) +
  ylab("AUC") + 
  geom_hline(yintercept = brca$auc[1], color="blue") +
  geom_text(aes(16.02, brca$auc[1], 
                label="AUC of true model", hjust=0.9, vjust=1.5), color="blue",
            size=5) +
  geom_hline(yintercept = brca$auc[2], color="blue", linetype=2) +
  geom_hline(yintercept = brca$auc[3], color="blue", linetype=2) + 
  theme_bw(base_size=20) + scale_shape(solid = FALSE) + theme(axis.title.x = element_blank(),
                                                              axis.text.x = element_blank())  #+ ylim(c(0.5, 0.7))


p5 = cowplot::plot_grid(p1b, p2, align = "v", ncol = 1, rel_heights = c(0.5, 0.5), axis = "l")
ggsave("../results/auc_corr_plots2.png", p5, scale=1.2, dpi=300)
