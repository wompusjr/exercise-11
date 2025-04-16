#####data loading####
library(tidyverse)
library(skimr)
library(naniar)
library(ggplot2)
library(MuMIn)
m <- read_tsv("https://raw.githubusercontent.com/difiore/ada-datasets/main/Mammal_lifehistories_v2.txt", col_names = TRUE)
skim(m) 
#categorical include: order, family, genus, species 
#numeric include: mass, gestation, newborn, weaning, wean mass, AFR, max. life, litter size, litters/year, refs
#####Step One####
m <- m |> replace_with_na_all(condition = ~.x == -999)

#####Step Two####
m <- m |> select(-"litter size", -"refs")
#also renaming the variables so I don't lose it
names(m)[names(m) == 'mass(g)'] <- 'mass'
names(m)[names(m) == 'gestation(mo)'] <- 'gestation'
names(m)[names(m) == 'newborn(g)'] <- 'newborn'
names(m)[names(m) == 'wean mass(g)'] <- 'wean_mass'
names(m)[names(m) == 'weaning(mo)'] <- 'weaning'
names(m)[names(m) == 'AFR(mo)'] <- 'AFR'
names(m)[names(m) == 'max. life(mo)'] <- 'max_life'
names(m)[names(m) == 'litters/year'] <- 'LPY'

#####Step Three####
m <- m |> 
  mutate(logMass = log(mass),
  logGest = log(gestation),
  logNBmass = log(newborn),
  logWeaning = log(weaning),
  logWeanMass = log(wean_mass),
  logAFR = log(AFR),
  logMaxLife = log(max_life),
  logLPY = log(LPY)
)

#####Step Four####
##age
#mass by gestation
md1 <- lm(logGest ~ logMass, data = m, na.action = na.exclude)
m <- m |> mutate(resGest = residuals(md1))
#mass by weaning
md2 <- lm(logWeaning ~ logMass, data = m, na.action = na.exclude)
m <- m |> mutate(resWeaning = residuals(md2))
#mass by age at first repro
md3 <- lm(logAFR ~ logMass, data = m, na.action = na.exclude)
m <- m |> mutate(resAFR = residuals(md3))
#mass by max lifespan
md4 <- lm(logMaxLife ~ logMass, data = m, na.action = na.exclude)
m <- m |> mutate(resMax_Life = residuals(md4))
##mass
#mass by newborn
md5 <- lm(logNBmass ~ logMass, data = m, na.action = na.exclude)
m <- m |> mutate(resNBmass = residuals(md5))
#mass by wean mass
md6 <- lm(logWeanMass ~ logMass, data = m, na.action = na.exclude)
m <- m |> mutate(resWeanMass = residuals(md6))

####Step Five####
(p1 <- ggplot(data = m, mapping = aes(x=order, y=resMax_Life)) +
  geom_boxplot())
##plot indicates that primates, scandentia, and xenarthra have the highest residual lifespan
(p2 <- ggplot(data = m, mapping = aes(x=order, y=resNBmass)) +
    geom_boxplot())
##plot is a lot more even but macroscelidae looks to be the highest on average
##although other groups have higher tails or outliers I believe the average and inter-quartile range is more telling
(p3 <- ggplot(data = m, mapping = aes(x=order, y=resWeanMass)) +
    geom_boxplot())
##plot suggests that perissodactyla is the beefiest at weaning age
##cetacea is pretty much next but carnivora, primates and lagamorpha have big tails

####Step Six####
m <- m |> drop_na(logGest, logNBmass, logWeaning, logWeanMass, logLPY, logMass)
#testing max.life
m2 <- m |> drop_na(logMaxLife)
maxfull <- lm(data = m2, logMaxLife ~ logGest + logNBmass + logWeaning + logWeanMass + logLPY + logMass, na.action = "na.fail")
(mods1 <- dredge(maxfull)) 
##model with logGest, logLPY, logMass, and logWeaning is best (a.k.a model 24)
##there are 5 models beneath 4 (model 24, 32 (everything but wean mass), 56 (everything but newborn mass), 64 (everything), 52 (exclude mass and newborn mass))
max.avg <- summary(model.avg(mods1, subset = delta <= 4, fit = TRUE))
confint(max.avg)
plot(max.avg)

#testing AFR
m3 <- m |> drop_na(logAFR)
afrfull <- lm(data = m3, logAFR ~ logGest + logNBmass + logWeaning + logWeanMass + logLPY + logMass, na.action = "na.fail")
(mods2 <- dredge(afrfull)) 
##model with logGest, logLPY, logMass, and logWeaning is best (a.k.a model 24)
##there are 7 models beneath delta = 4 
  ##model 24, 
  ##model 42 (excluding logmass and lognewborn)
  ##model 32 (everything but logweaning mass)
  ##model 56 (everything but lognewbonr)
  ##model 60 (everything but logmass)
  ##model 28 (excluding logmass and logweaning mass)
  ##model 64 (everything)
afr.avg <- summary(model.avg(mods2, subset = delta <= 4, fit = TRUE))
confint(afr.avg)
plot(afr.avg)