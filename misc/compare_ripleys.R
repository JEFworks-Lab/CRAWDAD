library(crawdad)
data(slide)
head(slide)

## run ripley's L
library(spatstat)
## need to convert to these point process objects
pp  <- ppp(x=slide$x, y=slide$y, 
           window = owin(range(slide$x), range(slide$y)), 
           marks=factor(slide$celltypes))
plot(pp)
## now can run for one pair
#test <- Lcross.inhom(pp,'UBCs','Granule')
#test <- Lcross(pp,'UBCs','Granule')
#test <- localLcross.inhom(pp,'UBCs','Granule')
test <- Kcross(pp,'UBCs','Granule')
#test <- Kcross.inhom(pp,'UBCs','Granule')

## "The standard interpretation of the plots of the K function is that
## if the estimated K function curve lies above the theoretical curve, 
## then the pattern is clustered, while if the estimated K lies below the theoretical curve, 
## then the pattern is regular."
plot(test)
names(test)
## plot isotropic-corrected minus theoretical across distances
plot(x = test$r, y = test$iso - test$theo, type='l')

## run on all
results <- do.call(rbind, lapply(levels(pp$marks), function(x) {
  test <- Kcross(pp,'Bergmann',x)
  ## isotropic-corrected minus theoretical
  test$iso - test$theo
}))
rownames(results) <- levels(pp$marks)
colnames(results) <- test$r

## plot
results.melt <- reshape2::melt(results)

library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
ggplot(results.melt, aes(x = Var2, y = value, col=Var1)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(aes(label = Var1),
                  nudge_x = 1,
                  force = 1,
                  box.padding = 1,
                  segment.alpha = .5,
                  data = results.melt %>%
                    group_by(Var1) %>% 
                    filter(Var2 == max(Var2)))

## focus on results of abundant cell-types?
vi <- names(which(table(slide$celltypes) > 20))
results.melt <- reshape2::melt(results[vi,])
ggplot(results.melt, aes(x = Var2, y = value, col=Var1)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(aes(label = Var1),
                  nudge_x = 1,
                  force = 1,
                  box.padding = 1,
                  segment.alpha = .5,
                  data = results.melt %>%
                    group_by(Var1) %>% 
                    filter(Var2 == max(Var2)))

## look at results
plt <- crawdad::vizAllClusters(cells = slide,
                               coms = as.factor(slide$celltypes),
                               ofInterest = c("Bergmann", "Purkinje"),
                               s = 2) 
plt


