breaks = (1/seq.int(min(dat$scale),max(dat$scale),length.out = 5))
breaks = scales::trans_breaks(identity, identity, n = 5)


labels = function(x) round(1/x)
(1/scale)


vizColocDotplot(dat)



# Test --------------------------------------------------------------------

sig_dat <- dat %>%
  dplyr::filter(abs(Z) >= 2) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE))

ggplot(sig_dat, aes(x=reference, y=neighbor, 
                    color=Z, size=scale)) + 
  geom_point() +
  scale_size_continuous(trans = "reverse")
  


# Reorder -----------------------------------------------------------------

## reorder
## get mean Z
mean_dat <- dat %>% 
  group_by(neighbor, scale, reference) %>% 
  summarize(Z = mean(Z))

## calculate sig z scores
sig_dat <- mean_dat %>%
  filter(abs(Z) >= 1.96) %>% 
  group_by(neighbor, reference) %>% 
  filter(scale == min(scale, na.rm = TRUE))
head(sig_dat)

## compare to the original code
## this code will get the zscore for all permutations that present a significant
## result, then they will be displayed. So, if only a single permutation is
## significant, the colocalization will also be.
o_sig_dat <- dat %>%
  filter(abs(Z) >= 1.96) %>% 
  group_by(neighbor, reference) %>% 
  filter(scale == min(scale, na.rm = TRUE))
head(o_sig_dat)
