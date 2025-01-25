#BS852 Final Project
library(tidyverse)
library(kableExtra)
library(olsrr)
library(car)
library(GGally)
library(survival)
library(survminer)

data0 <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS852 Statistical Methods in Epidemiology/Final Project/framdat4.csv")
view(data0)
#attach(data0)

data <- data0  %>% filter(!is.na(CHOL4)&!is.na(CIGS4)&!is.na(WGT4)&!is.na(FVC4)&!is.na(BMI4)&!is.na(HTN4)&!is.na(CHD)&!is.na(T2D))
attach(data)

data$AGE4_CAT <- ifelse(data$AGE4 >= 34 & data$AGE4 < 42, 1,
                        ifelse(data$AGE4 >= 42 & data$AGE4 < 49, 2,
                               ifelse(data$AGE4 >= 49 & data$AGE4 < 56, 3, 
                                      ifelse(data$AGE4 >= 56 & data$AGE4 <= 69, 4, NA))))

data$SPF4_CAT <- ifelse(data$SPF4 <= median(data$SPF4, na.rm = TRUE), "Low", "High")
data$SPF4_CAT <- as.factor(data$SPF4_CAT)
str(data)
data$SPF4_CAT <- factor(data$SPF4_CAT, levels = c("Low", "High"), labels = c(1, 2))
levels(data$SPF4_CAT)
#Columns with NA data in the original dataset
sum(is.na(CHOL4))
sum(is.na(CIGS4))
sum(is.na(WGT4))
sum(is.na(FVC4))
sum(is.na(BMI4))
sum(is.na(HTN4))
sum(is.na(SMOKE))
sum(is.na(CHD))
sum(is.na(CHD_SURV))
sum(is.na(T2D))
sum(is.na(T2D_SURV))
#removing NA individuals

#from 2000 individuals to 1506


#
mod <- glm(DTH ~ SEX + AGE4 + CHOL4 + CIGS4 + SPF4 + DPF4 + WGT4 + FVC4 + BMI4 + HTN4 + SMOKE, family = "binomial", data = data)
summary(mod)
#exploring collinearity and correlation
ols_vif_tol(mod)
vif(mod)
ols_eigen_cindex(mod)
#ggpairs(data)



#1 Survival analysis
#Variable Selection, backwards, AIC

#Reduced Model
fit0<- coxph(Surv(SURV, DTH)~SPF4, data = data);summary(fit0)
extractAIC(fit0)

#Univariate
coxph(Surv(SURV, DTH)~SPF4 + AGE4, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +SEX, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +CIGS4, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +CHOL4, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +CHD, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +DPF4, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +WGT4, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +BMI4, data = data0)
coxph(Surv(SURV, DTH)~SPF4 +HTN4, data = data0)

#Full model
fit1 <- coxph(Surv(SURV, DTH)~SPF4 + SEX + AGE4 + CHOL4 + CIGS4 + T2D + CHD + DPF4 + WGT4 + FVC4 + BMI4 + HTN4 + SMOKE)
extractAIC(fit1)

#Stepwise Analysis to find coefficients
a <- step(fit1, scope=list(lower=fit0, upper=fit1), direction="backward", data=data)
fit_sum <- summary(a)$coefficients; fit_sum
fit <- coxph(Surv(SURV, DTH)~SPF4 + SEX + AGE4 + CIGS4 +  CHD + FVC4, data = data)
summary(fit)

cox.zph(fit)

fit.final <- coxph(Surv(SURV, DTH)~SPF4 + SEX + AGE4 + CIGS4 + FVC4, data = data)
summary(fit.final)
cox.zph(fit.final)


fit.final2 <- coxph(Surv(SURV, DTH)~SPF4_CAT + SEX + AGE4 + CIGS4 + FVC4, data = data)
summary(fit.final2)
cox.zph(fit.final2)


#2
summary(data$SPF4)

#SEX CONFOUNDER?
fit.test <- coxph(Surv(SURV, DTH)~SPF4 + SEX , data = data); summary(fit.test)
((exp(fit0$coefficients[1]*50) - exp(fit.test$coefficients[1]*50)) / exp(fit0$coefficients[1]*50))*100
((exp(fit0$coefficients[1]) - exp(fit.test$coefficients[1])) / exp(fit0$coefficients[1]))*100
fit.testb <- coxph(Surv(SURV, DTH)~SPF4 + FVC4 + AGE4 + CIGS4 , data = data); summary(fit.testb)
((exp(fit.testb$coefficients[1]*50) - exp(fit.final$coefficients[1]*50)) / exp(fit.testb$coefficients[1]*50))*100

#AGE CONFOUNDER?
fit.test2 <- coxph(Surv(SURV, DTH)~SPF4 + AGE4 , data = data); summary(fit.test2)
((exp(fit0$coefficients[1]*50) - exp(fit.test2$coefficients[1]*50)) / exp(fit0$coefficients[1]*50))*100
((exp(fit0$coefficients[1]) - exp(fit.test2$coefficients[1])) / exp(fit0$coefficients[1]))*100
fit.test2b <- coxph(Surv(SURV, DTH)~SPF4 + SEX + FVC4 + CIGS4 , data = data); summary(fit.test2b)
((exp(fit.test2b$coefficients[1]*50) - exp(fit.final$coefficients[1]*50)) / exp(fit.test2b$coefficients[1]*50))*100

#CIGS COnFOUNDER?
fit.test3 <- coxph(Surv(SURV, DTH)~SPF4 + CIGS4 , data = data); summary(fit.test3)
((exp(fit0$coefficients[1]*50) - exp(fit.test3$coefficients[1]*50)) / exp(fit0$coefficients[1]*50))*100
((exp(fit0$coefficients[1]) - exp(fit.test3$coefficients[1])) / exp(fit0$coefficients[1]))*100
fit.test3b <- coxph(Surv(SURV, DTH)~SPF4 + SEX + AGE4 + FVC4 , data = data); summary(fit.test3b)
((exp(fit.test3b$coefficients[1]*50) - exp(fit.final$coefficients[1]*50)) / exp(fit.test3b$coefficients[1]*50))*100

#FVC4 CONFOUNDER?
fit.test4 <- coxph(Surv(SURV, DTH)~SPF4 + FVC4 , data = data); summary(fit.test4)
((exp(fit0$coefficients[1]*50) - exp(fit.test4$coefficients[1]*50)) / exp(fit0$coefficients[1]*50))*100
((exp(fit0$coefficients[1]) - exp(fit.test4$coefficients[1])) / exp(fit0$coefficients[1]))*100
fit.test4b <- coxph(Surv(SURV, DTH)~SPF4 + SEX + AGE4 + CIGS4 , data = data); summary(fit.test4b)
((exp(fit.test4b$coefficients[1]*50) - exp(fit.final$coefficients[1]*50)) / exp(fit.test4b$coefficients[1]*50))*100


#JOINT CONFOUNDERs
((exp(fit0$coefficients[1]*50) - exp(fit.final$coefficients[1]*50)) / exp(fit0$coefficients[1]*50))*100
((exp(fit0$coefficients[1]) - exp(fit.final$coefficients[1])) / exp(fit0$coefficients[1]))*100



#3
fit3<- coxph(Surv(SURV, DTH)~SPF4+SEX + SPF4*SEX, data = data);summary(fit3)
cox.zph(fit3)

fit3b<- coxph(Surv(SURV, DTH)~SPF4+SEX + SPF4*SEX + AGE4 + FVC4 + CIGS4, data = data);summary(fit3b)
cox.zph(fit3b)

male <- data %>% filter (SEX == 1)
female <- data %>% filter(SEX ==2)
fitm<- coxph(Surv(SURV, DTH)~SPF4 + AGE4 + CIGS4 + FVC4, data = male);summary(fitm)
fitf <- coxph(Surv(SURV, DTH)~SPF4 + AGE4 + CIGS4 + FVC4, data = female);summary(fitf)


#4
fitchd <- coxph(Surv(SURV,DTH)~CHD + SEX + AGE4 + FVC4 + CIGS4 + SPF4, data = data); summary(fitchd)

#Summary Table
summary_table <- data %>%
  summarise(
    `Sex (Male %)` = paste0(sum(SEX == 1), " (", round(mean(SEX == 1) * 100, 1), "%)"),
    `Age (Mean ± SD)` = paste0(round(mean(AGE4, na.rm = TRUE), 1), " ± ", round(sd(AGE4, na.rm = TRUE), 1)),
    `FVC4 (Mean ± SD)` = paste0(round(mean(FVC4, na.rm = TRUE), 1), " ± ", round(sd(FVC4, na.rm = TRUE), 1)),
    `Cigs4 (Median [IQR])` = paste0(median(CIGS4, na.rm = TRUE), " [", 
                                    quantile(CIGS4, 0.25, na.rm = TRUE), ", ", 
                                    quantile(CIGS4, 0.75, na.rm = TRUE), "]"),
    `SPF4 (Mean ± SD)` = paste0(round(mean(SPF4, na.rm = TRUE), 1), " ± ", round(sd(SPF4, na.rm = TRUE), 1)),
    `Deaths (%)` = paste0(sum(DTH == 1), " (", round(mean(DTH == 1) * 100, 1), "%)"),
    `Survival Time (Median [IQR])` = paste0(median(SURV, na.rm = TRUE), " [", 
                                            quantile(SURV, 0.25, na.rm = TRUE), ", ", 
                                            quantile(SURV, 0.75, na.rm = TRUE), "]")
  ) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Summary")

kable(summary_table, 
      caption = "Patient Characteristics",
      align = "l") %>%
  kable_styling(full_width = F)



km_fit <- survfit(Surv(SURV, DTH) ~ AGE4_CAT + SPF4_CAT, data = data)

# Plot the Kaplan-Meier curve stratified by AGE4_CAT and SPF4_CAT
ggsurvplot(km_fit, data = data, pval = T, risk.table = T,      
xlab = "Survival Time (Days)", 
ylab = "Survival Probability", 
legend.title = "Age & SPF4 Group",    
legend.labs = c("34-42, Low SPF4", "34-42, High SPF4", "42-49, Low SPF4", "42-49, High SPF4", "49-56, Low SPF4", "49-56, High SPF4", "56-69, Low SPF4", "56-69, High SPF4"), 
palette = c("lightblue", "#2E9FDF", "lightpink", "#FF33CC", "lightgreen", "green", "red", "darkred")) 

