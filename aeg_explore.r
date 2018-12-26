library(tidyverse)
library(rms)

df <- read_csv('AEG_XH_2018.csv')

# Var transformation
ndf <- df

ndf$name <- gsub("[0-9（）]","", df$姓名_1)
ndf <- mutate(ndf, sex = ifelse(性别==1, 'M', 'F'))
ndf$sex <- factor(ndf$sex)
ndf$age <- ndf$年龄
ndf$tumsiz <- ndf$COX肿瘤大小
ndf$tumsiz[ndf$COX肿瘤大小 == '#NULL!'] <- NA
ndf$tumsiz <- as.numeric(ndf$tumsiz)
ndf$grade <- factor(ndf$COX分化程度代码, 
                    levels = c(1, 2, 3, 3.00, 4),
                    labels = c(1, 2, 3, 3, 4))
ndf$vasinv <- factor(ndf$COX脉管浸润,
                     levels = c(0, 1),
                     labels = c('No', 'Yes'))
ndf$neuinv <- factor(ndf$神经浸润,
                     levels = c(0, 1),
                     labels = c('No', 'Yes'))
ndf$r0 <- factor(ndf$COX_Rstatus,
                 levels = c(0, 1),
                 labels = c('No', 'Yes'))
ndf$dajcct <- factor(ndf$COX_T分级代码,
                 levels = c(1, 2, 3, 4, 5),
                 labels = c('T1', 'T2', 'T3', 'T4a', 'T4b'))
ndf$pn <- as.numeric(ndf$COX阳性淋巴结数)
ndf$dajccn <- factor(ndf$N分级,
                     levels = c(0, 1, 2, 3, 4),
                     labels = c('N0', 'N1', 'N2', 'N3a', 'N3b'))
ndf$ne <- as.numeric(ndf$淋巴结数目)
## no M1 patient observed
ndf$siewiert <- factor(ndf$肿瘤部位,
                       levels = c(1, 2, 3),
                       labels = c('SI', 'SII', 'SIII'))
ndf$surgprim <- factor(ndf$手术切除范围,
                       levels = c(1, 2),
                       labels = c('PG', 'TG'))
ndf$apprch <- factor(ndf$手术入路,
                     levels = c(1, 2, 3),
                     labels = c('TT', 'TM', 'Combined'))
ndf$srv_mon <- as.numeric(ndf$生存月份)
ndf$stat_rec <- factor(ndf$生存状态,
                       levels = c(0, 1),
                       labels = c('Dead', 'Survived'))
ndf$op_year <- paste0('20', substr(ndf$手术日期, 8,9))
ndf$op_year <- as.numeric(ndf$op_year)
ndf$inhos_drt <- ndf$住院时间
ndf$surg_drt <- as.numeric(ndf$手术时间分钟)
ndf$surg_bld <- as.numeric(ndf$出血)
ndf$po_bld <- factor(ndf$并发症出血,
                     levels = c('00', 0, 1),
                     labels = c('No', 'No', 'Yes'))
ndf$inc_inf <- factor(ndf$切口感染,
                      levels = c(0, 1),
                      labels = c('No', 'Yes'))
ndf$lung_inf <- factor(ndf$肺部感染,
                      levels = c(0, '00', 1),
                      labels = c('No', 'No', 'Yes'))
ndf$inc_inf <- factor(ndf$腹腔感染,
                      levels = c(0, 1),
                      labels = c('No', 'Yes'))
ndf$fstl <- factor(ndf$吻合口瘘,
                   levels = c(0, 1, '1（腹腔感染）', 'N'),
                   labels = c('No', 'Yes','Yes', 'No'))
ndf$pgs <- factor(ndf$胃瘫,
                  levels = c(0, 1),
                  labels = c('No', 'Yes'))
ndf$her2 <- factor(ndf$`HER-2`,
                   levels = c(0, 'cerbB2-', 'N','无扩增','cerbB2+', 'cerbB2++', 'cerbB2+++','Y',1),
                   labels = c('-', '-', '-', '-', '+', '+', '+', '+', '+'))
ndf$chemo <- ifelse(ndf$化疗==0,0,1)
ndf$chemo <- factor(ndf$chemo,
                    levels = c(0, 1),
                    labels = c('No', 'Yes'))
ndf <- select(ndf,
              'name':'chemo', -'name')

## 5 pts with survival month = 0 were perioperative deaths.
## when predicting survival, these pts should be excluded.
ndf <- filter(ndf, srv_mon!=0)

## too few pts with g1 and g4
ndf$grade2c <- ifelse(as.numeric(ndf$grade)<=2, '1-2', '3-4')
ndf$grade2c <- factor(ndf$grade2c)
##
ndf$dajcct2c <- ifelse(ndf$dajcct=='T4b'|ndf$dajcct=='T4a', 4, ndf$dajcct)
ndf$dajcct2c <- factor(ndf$dajcct2c,
                       levels = c(1, 2, 3, 4),
                       labels = c('T1', 'T2', 'T3', 'T4'))
ndf$dajccn2c <- ifelse(ndf$dajccn=='N3b'|ndf$dajccn=='N3a', 4, ndf$dajccn)
ndf$dajccn2c <- factor(ndf$dajccn2c,
                       levels = c(1, 2, 3, 4),
                       labels = c('N0', 'N1', 'N2', 'N3'))

## Event:177
## basic strategy: seperate the cohort by op_year (04-12 and 13-16),
## thus we got 121/235 + 56/99 events.
## we shall use the former to model and the latter to validate.
## by the rule of thumb (Harell), the maximum of d.f should be exceed 121/15 = 8

## Vars included in the model:
## basic: age(2), sex(1), siewiert(2)
## classic: tumsiz, dajcct2c(3), dajccn2c(3), grade2c(3), r0(1)
## explorative: pn/ne, complication, surgprim, apprch, 
## too many NAs: her2(88), chemo(91)

df_train <- filter(ndf, op_year <= 2012)
df_validate <- filter(ndf, op_year > 2012)

units(ndf$srv_mon) <- 'Month'
dd <- datadist(ndf)
options(datadist = 'dd')

S <- Surv(df_train$srv_mon, df_train$stat_rec == 'Dead')
fit_lm_1 <- cph(S ~ rcs(age,3) + sex + siewiert + dajcct2c + dajccn2c + grade2c + r0 + rcs(tumsiz, 3), 
    data = df_train, x=T, y=T, surv=T)
estimates <- survest(fit_lm_1, newdata = as.data.frame(df_validate), times=12)
S_val <- Surv(df_validate$srv_mon, df_validate$stat_rec == 'Dead')
rcorr.cens(x=estimates$surv,S=S_val)
validated=val.surv(fit_lm_1,newdata=as.data.frame(df_validate),S=S_val)
plot(validated)

