
# plot predictions for family from family_interactions.R
# predictions calculated in pred_fam.ipynb

### mom, mgm
# 0, 0
# 0, 80b
# 0 , 60b
# 50b, 60b
# 50b and 60o, 60b
# 50b and 60o, 40b

ind = 1:5

fam.pred = read.csv('../results/fam_pred.csv')

pred = data.frame(scenario=c("A", "B" ,"C", "D", "E"),
                  FCNN=fam.pred$fcnn[ind],
                  CNN=fam.pred$cnn[ind],
                  LR=fam.pred$lr[ind], 
                  BRCAPRO=fam.pred$brcapro[ind])

library(reshape2)
pred2 <- melt(pred, id.vars=c("scenario"))
names(pred2)[2] = "Model"

library(ggplot2)
ggplot(pred2, aes(x=scenario, y=value, colour=Model)) + geom_point(size=4.5, aes(shape=Model), stroke=1.25) +
  xlab("Scenario") +
  ylab("10-year Risk") + theme_bw(base_size=20) + scale_color_manual(values = c("#00BFC4","#F8766D", "#C77CFF", "black")) + scale_shape_manual(values=c(2, 1, 0, 4))


ggsave("family_interactions.png", scale=1, dpi=300)

