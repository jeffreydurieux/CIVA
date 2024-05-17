# Wed Nov  8 12:00:03 2023 ------------------------------
# JD
# What: analyse CIVA sim 1 results

library(doBy)


#### simulation 1 #####
setwd("C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation1/")
load(file = 'allres.Rdata')

summary(alldf$ari)
summary(alldf$perc)
summary(alldf$Stuck)
summary(alldf$Atuck)
summary_by(data = alldf, formula = Stuck~R*Q*Nr*Err)
summary_by(data = alldf, formula = Atuck~R*Q*Nr*Err)

summary_by(data = alldf, formula = ari~R, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ari~Q, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ari~Err, FUN = c(mean,sd)) %>% round(digits = 3)

summary_by(data = alldf, formula = Stuck~R, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = Stuck~Q, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = Stuck~Err, FUN = c(mean,sd)) %>% round(digits = 3)

summary_by(data = alldf, formula = Atuck~R, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = Atuck~Q, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = Atuck~Err, FUN = c(mean,sd)) %>% round(digits = 3)



df=summary_by(data = alldf, formula = Stuck~R*Q*Nr*Err)
id1 = which(df$Err==.1)
id2 = which(df$Err==.3)
id3 = which(df$Err==.6)
p1=plot_ly(data = df[id1,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter') 
p2=plot_ly(data = df[id2,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
p3=plot_ly(data = df[id3,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
subplot(p1,p2,p3, shareY = TRUE)


#### simulation 2 #####
setwd("C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation2/")
load(file = 'output_df_Sim2.Rdata')
alldf <- df
summary(alldf$ari)
summary(alldf$STuck);sd(alldf$STuck)
summary(alldf$ATuck);sd(alldf$ATuck)
summary_by(data = alldf, formula = STuck~E*Ti)
summary_by(data = alldf, formula = ATuck~E*Ti)

summary_by(data = alldf, formula = ari~E, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ari~Ti, FUN = c(mean,sd)) %>% round(digits = 3)

summary_by(data = alldf, formula = STuck~E, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = STuck~Ti, FUN = c(mean,sd)) %>% round(digits = 3)

summary_by(data = alldf, formula = ATuck~E, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ATuck~Ti, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ATuck~Ti*E, FUN = c(mean,sd)) %>% round(digits = 3)



df=summary_by(data = alldf, formula = Stuck~R*Q*Nr*Err)
id1 = which(df$Err==.1)
id2 = which(df$Err==.3)
id3 = which(df$Err==.6)
p1=plot_ly(data = df[id1,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter') 
p2=plot_ly(data = df[id2,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
p3=plot_ly(data = df[id3,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
subplot(p1,p2,p3, shareY = TRUE)

#### simulation 3 #####

setwd("C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation3/")
load(file = 'output_df_Sim3.Rdata')
alldf <- df

summary(alldf$ari); sd(alldf$ari)
summary(alldf$STuck);sd(alldf$STuck)
summary(alldf$ATuck);sd(alldf$ATuck)
summary_by(data = alldf, formula = STuck~E*Ti)
summary_by(data = alldf, formula = ATuck~E*Ti)

summary_by(data = alldf, formula = ari~E, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ari~Ti, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ari~Ti*E, FUN = c(mean,sd)) %>% round(digits = 3)

summary_by(data = alldf, formula = STuck~E, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = STuck~Ti, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = STuck~Ti*E, FUN = c(mean,sd)) %>% round(digits = 3)

summary_by(data = alldf, formula = ATuck~E, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ATuck~Ti, FUN = c(mean,sd)) %>% round(digits = 3)
summary_by(data = alldf, formula = ATuck~Ti*E, FUN = c(mean,sd)) %>% round(digits = 3)


df=summary_by(data = alldf, formula = Stuck~R*Q*Nr*Err)
id1 = which(df$Err==.1)
id2 = which(df$Err==.3)
id3 = which(df$Err==.6)
p1=plot_ly(data = df[id1,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter') 
p2=plot_ly(data = df[id2,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
p3=plot_ly(data = df[id3,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
subplot(p1,p2,p3, shareY = TRUE)