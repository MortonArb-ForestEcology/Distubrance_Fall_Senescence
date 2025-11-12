library(nlme)
library(MuMIn)  # for R-squared calculations

# Load the main analysis # CR: Don't do this.  Save an output from that script and load it here!
source("hansen_distubence_explore.R")

#filter data for +/- 3 years 
dat.lme <- flossdvi[abs(flossdvi$yrfromdisturb) <= 5 & 
                       !is.na(flossdvi$MidGreendown_DOY) & 
                       !is.na(flossdvi$yrfromdisturb), ]

# Check data structure to make sure there aren't glaring errors 
str(flossdvi[, c("Year", "Label", "MidGreendown_DOY", "dstrbyr", "yrfromdisturb")])
head(dat.lme)

# Check sample sizes to ensure there are an equal number of labels (sample size) and no .Na values 
length(unique(flossdvi$Label))
table(flossdvi$Label)
sum(is.na(flossdvi$MidGreendown_DOY))

# Look at the range of Years_From_Disturbance
range(flossdvi$yrfromdisturb, na.rm = TRUE)

#Model accounting for weather and site as a random effect
lmemod <- lme(MidGreendown_DOY ~ as.factor(yrfromdisturb), random = list(Year = ~1, Label = ~1), data = dat.lme)

#testing for signifigance
summary(lmemod)
r.squaredGLMM(lmemod)
anova(lmemod)
#aaaannnnd there does not seem to be

#model for just weather 'oops all weather'
lmemodwe <- lme(MidGreendown_DOY ~ as.factor(yrfromdisturb),random = ~1|Year, data = dat.lme)

# signif testing
summary(lmemodwe)
r.squaredGLMM(lmemodwe)
anova(lmemodwe)

# again no signifigance 

#model for just site 'oops all sites'
lmemodsit <- lme(MidGreendown_DOY ~ as.factor(yrfromdisturb), random = ~1|Label, data = dat.lme)

summary(lmemodsit)
r.squaredGLMM(lmemodsit)
anova(lmemodsit)

#three strikes
("you're out")



# Test to see if NDVI is lower in the year of disturbance to make sure something isn't looking weird
ggplot(data=dat.lme) +
  geom_line(aes(x=yrfromdisturb, y=NDVI, color=Label)) +
  geom_vline(xintercept=0) +
  theme_bw()


NDVItest <- lme(NDVI ~ relevel(as.factor(yrfromdisturb), "0"), random = list(Year = ~1, Label = ~1), data = dat.lme)
summary(NDVItest)
anova(NDVItest)

NDVItest2 <- lme(NDVI ~ as.factor(yrfromdisturb), random = list(Label = ~1), data = dat.lme)
summary(NDVItest2)
anova(NDVItest2)

# Looking at groupign pre & post disturbance 
dat.lme$yrGroup[dat.lme$yrfromdisturb<0] <- "before"
dat.lme$yrGroup[dat.lme$yrfromdisturb==0] <- "during"
dat.lme$yrGroup[dat.lme$yrfromdisturb>0] <- "after"
dat.lme$yrGroup <- factor(dat.lme$yrGroup, levels=c("before", "during", "after"))
summary(dat.lme)

timeGroupLME <- lme(MidGreendown_DOY ~ yrGroup, random = list(Year=~1, Label = ~1), data = dat.lme)
summary(timeGroupLME)
anova(timeGroupLME)
