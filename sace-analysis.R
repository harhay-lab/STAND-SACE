# Load libraries
library(tidyverse)
library(lme4)

##### Replicate STAND data primary outcome analysis #########

stdat=read.csv('/Volumes/stand trial/Data/Analysis/StartingFiles/STAND_SampleFile_211123.csv', header=T)
stdat$time= as.factor(stdat$step)

### Marital status
stdat$zMari_status=ifelse(stdat$zMarital_status=='Partner'|stdat$zMarital_status== 'Married','Married/Partner', ifelse(stdat$zMarital_status=='Divorced'|stdat$zMarital_status== 'Separated','Divorced/Seperated', stdat$zMarital_status))

stdat$zMari_status=relevel(as.factor(stdat$zMari_status), ref='Married/Partner')

st_prim=stdat[stdat$sampleflag_1==1 & is.na(stdat$max_ims_adj_death)==F,]

st_prim$unit_name2 <-
  ifelse(st_prim$unit_name %in% c('PAH CCU WIDENER 3A','PAH ICU SCHIEDT 3',
                                  'PAH NCU WIDENER 3B'),
         'PAH',st_prim$unit_name)

#### Unadjusted analysis using LMM with primary sample

ims_mod_lmm=lmer(max_ims_adj_death~ phase + time + (1|unit_name2), data=st_prim)
summary(ims_mod_lmm)

no.clusters= length(unique(stdat$unit_name))
tstat_lmm=fixef(ims_mod_lmm)[2]/sqrt(vcov(ims_mod_lmm)[2,2])
pval=2 * (pt(-abs(tstat_lmm),df= no.clusters-1))
pval

library(broom.mixed)
tidy(ims_mod_lmm, conf.int=T)

qqnorm(resid(ims_mod_lmm))
randoms<-ranef(ims_mod_lmm, condVar = TRUE)
qq <- attr(ranef(ims_mod_lmm, condVar = TRUE)[[1]], "postVar")

library(redres)
plot_ranef(ims_mod_lmm)


######################################
####### adjusted analysis ###########
ims_adj_lmm=lmer(max_ims_adj_death~ phase+zSex+ time + (1|unit_name2) +AGE_AT_ADMSN+ zrace+zhisp+zMari_status+icu_admit_source+days_since_launch+LAPS2+rass_ind+iculos_strict, data=st_prim)
summary(ims_adj_lmm)

tstat_lmm=fixef(ims_adj_lmm)[2]/sqrt(vcov(ims_adj_lmm)[2,2])
pval=2 * (pt(-abs(tstat_lmm),df= no.clusters-1))
pval


tidy(ims_adj_lmm,confint=T )

qqnorm(resid(ims_adj_lmm))
qqline(resid(ims_adj_lmm))

plot_ranef(ims_adj_lmm)


#################################################################
# Supplemental SACE analysis

# Model (instead use unadjusted comparison as in Chiba Vanderweele 2013)
iculoslogm_adj <-
  lmer(log(iculos_outcome_readm) ~ phase + time + (1|unit_name2) +
         AGE_AT_ADMSN + zrace + zhisp + zSex + zMari_status +
         factor(icu_admit_source) + days_since_launch + LAPS2 + rass_ind,
       data=st_prim)
summary(iculoslogm_adj)

# Unadjusted
mod <- t.test(st_prim$iculos_outcome_readm[st_prim$phase == "UC" &
                                           st_prim$icu_death_readm == 0],
              st_prim$iculos_outcome_readm[st_prim$phase == "IN" &
                                           st_prim$icu_death_readm == 0])

plot(seq(-120, 20, by = 1),
     mod$estimate[2] - mod$estimate[1] + seq(-120, 20, by = 1),
     type = "l", xlab = "Sensitivity parameter",
     ylab = "Difference in avg. length of stay")

lines(seq(-120, 20, by = 1), -mod$conf.int[2] + seq(-120, 20, by = 1),
      type = "l", lty = 2)
lines(seq(-120, 20, by = 1), -mod$conf.int[1] + seq(-120, 20, by = 1),
      type = "l", lty = 2)
abline(h = 0, lty = 3)
lines(-(mod$estimate[2] - mod$estimate[1]))

df <- data.frame(alpha = seq(-120, 20, by = 1),
                 ests = mod$estimate[2] - mod$estimate[1] +
                   seq(-120, 20, by = 1),
                 lower = -mod$conf.int[2] + seq(-120, 20, by = 1),
                 upper = -mod$conf.int[1] + seq(-120, 20, by = 1))

library(ggplot2)
p1 <- ggplot(data = df, aes(x = alpha, y = ests)) +
  geom_line() +
  geom_point(aes(x = 0, y = mod$estimate[2] - mod$estimate[1])) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_segment(y = 0, yend = -80, x = mod$estimate[1] - mod$estimate[2],
               xend = mod$estimate[1] - mod$estimate[2], lty = 3) +
  geom_segment(y = 0, yend = -80, x = mod$conf.int[2],
               xend = mod$conf.int[2], lty = 3) +
  annotate("text", x = -80, y = -50,
           label = paste("\u03b1 =",
                         round(mod$estimate[1] - mod$estimate[2], 1))) +
  annotate("text", x = -36, y = -50,
           label = paste(expression("\u03b1 ="),
                         round(mod$conf.int[2], 1))) +
  xlab("Sensitivity parameter \u03b1") +
  ylab("SACE") +
  theme_bw()

ggsave(filename = "SACE-plot.pdf", plot = p1, device = cairo_pdf)
