# Wed Nov  8 12:00:03 2023 ------------------------------
# JD
# What: analyse CIVA sim 1 results

library(doBy)

setwd("C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation1/")
load(file = 'allres.Rdata')

summary(alldf$ari)
summary(alldf$perc)
summary(alldf$Stuck)
summary(alldf$Atuck)
summary_by(data = alldf, formula = Stuck~R*Q*Nr*Err)
summary_by(data = alldf, formula = Atuck~R*Q*Nr*Err)


df=summary_by(data = alldf, formula = Stuck~R*Q*Nr*Err)
id1 = which(df$Err==.1)
id2 = which(df$Err==.3)
id3 = which(df$Err==.6)
p1=plot_ly(data = df[id1,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter') 
p2=plot_ly(data = df[id2,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
p3=plot_ly(data = df[id3,], x=~jitter(Q, amount = .1), y = ~Stuck.mean,color=~as.factor(R),type = 'scatter')
subplot(p1,p2,p3, shareY = TRUE)
