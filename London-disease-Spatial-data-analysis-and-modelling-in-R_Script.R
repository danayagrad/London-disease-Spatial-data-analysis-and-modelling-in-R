#######################Spatial Statistics Assignment########################

library(spdep)
library(sp)
library(R2OpenBUGS)
library(CARBayes)
library(coda)


###################################  1  ###################################
#1.a

load("london_3.RData")
head(london_3@data)
names(london_3@data)

london3s = london_3[london_3$year == "2003",c("Observed","Expected","JSA","Price","PM25","year")]
head(london3s@data)
NROW(london3s@data)
names(london3s@data)

#1.b
sir = london3s$Observed/london3s$Expected
head(sir)
NROW(sir)
class(sir)


london3s$Sir = sir
names(london3s@data)
head(london3s@data)

#1.3
spplot(london3s,"Sir", colorkye = T, main=list(label="SIR of Respiratory disease admission across London",cex=1))

#1.4
par(mfrow = c(1,3))
# Sir and PM25
plot(london3s$Sir~london3s$PM25, data = as.data.frame(london3s), col = "Purple", xlab = "PM25", ylab = "SIR", main = "Association between SIR and PM25", cex.main=1)
abline(lm(london3s$Sir~london3s$PM25), col = "Purple2", lty= 2, lwd=2)

# Sir and JSA
plot(london3s$Sir~london3s$JSA, data = as.data.frame(london3s), col = "slateblue1", xlab = "JSA", ylab = "SIR", main = "Association between SIR and JSA", cex.main=1)
abline(lm(london3s$Sir~london3s$JSA), col = "slateblue2", lty= 2, lwd=2)

# Sir and Price
plot(london3s$Sir~london3s$Price, data = as.data.frame(london3s), col = "maroon3", xlab = "Price", ylab = "SIR", main = "Association between SIR and Price", cex.main=1)
abline(lm(london3s$Sir~london3s$Price), col = "maroon4", lty= 2, lwd=2)


#############################2. Poisson Regression##################################
#2.a 
## Assign data open bugs.

Data = list(Y = london3s$Observed, E = london3s$Expected, PM25 = as.numeric(scale(london3s$PM25)), JSA = as.numeric(scale(london3s$JSA)), Price= as.numeric(scale( london3s$Price)), N = length(london3s$Expected))

Inits = list(list(beta0 = 0, beta1 =0, beta2 = 0, beta3 =0))



## run MCMC sample using bugs function.
preg.bugs = bugs(data = Data,
                 model.file = "poisson-normal_draft2.txt",
                 parameters.to.save = c("beta0", "beta1", "beta2", "beta3"),
                 inits = rep(Inits,4),
                 n.iter = 5000,
                 n.burnin = 2500,
                 n.thin = 1,
                 n.chains = 4)


## trace plot for beta 0-3.
par(mfrow = c(2,2))
plot(preg.bugs$sims.list$beta0, type ="l", ylab = "Samples beta0")
plot(preg.bugs$sims.list$beta1, type ="l", ylab = "Samples beta1")
plot(preg.bugs$sims.list$beta2, type ="l", ylab = "Samples beta2")
plot(preg.bugs$sims.list$beta3, type ="l", ylab = "Samples beta3")


## Posterior median and credible intervals.
print(preg.bugs, digits.summary = 3)
preg.bugs$summary[2,1:7]


###2b: Use Gelman-Rubin diagostic (Rhat) and geweke diagnostic plots to assess convergence.

#### Gelman-rubin, checking RHAT
par(mfrow = c(1,1))
preg.bugs$summary[,8]
hist(preg.bugs$summary[,8], main = " Distribution of Rhat", ylab = " Frequency", xlab = "Rhat",col= c("slategray1", "lightblue", "lavender"))
### it is between 1 -1.2 --> converge.

### geweke convergence

geweke.plot(as.mcmc(as.data.frame(preg.bugs$sims.list$beta0)))
geweke.plot(as.mcmc(as.data.frame(preg.bugs$sims.list$beta1)))
geweke.plot(as.mcmc(as.data.frame(preg.bugs$sims.list$beta2)))
geweke.plot(as.mcmc(as.data.frame(preg.bugs$sims.list$beta3)))

###2c: Calculate pearson residuals and check 2 model assumptions
# Overdispersion

preg.bugs$summary[,1]
mu =london3s$Expected*exp(-2.926854e-01 + 3.598773e-02*as.numeric(scale(london3s$PM25)) +1.784286e-01*as.numeric(scale(london3s$JSA)) -3.472947e-02*as.numeric(scale( london3s$Price)))
pearson = (london3s$Observed - mu)/sqrt(mu)

## check overdispersion - plot fitted values and residual

par(mfrow = c(1,1))
plot(pearson ~ mu, xlab = "Fitted values", ylab = "Pearson residual", main = "Pearson residuals vs mu", col = "mediumorchid4")
var(pearson)

## check independence - moran's i statistics
Wnb <- nb2listw(poly2nb(london3s))
moran =moran.test(pearson, Wnb)
moran$estimate[1]
moran$p.value

########################3 Poisson Regression with CAR #########################
## 3.a  
#proudce the adjacency matrix, and various pieces of associated info Bugs needs
W = nb2mat(poly2nb(london3s), style = "B")
inds = lapply(1:nrow(W), function(i) which (W[i, ] ==1))
Adj = Reduce("c", inds)
Num.Adj = rowSums(W)
SumNumNeigh = sum(Num.Adj)

# combine all of the data, and constatns into a single list
DataCar= list(observed = london3s$Observed, expected = london3s$Expected, N = nrow(london3s), PM25 = as.numeric(scale(london3s$PM25)), JSA = as.numeric(scale(london3s$JSA)), Price= as.numeric(scale( london3s$Price)), Adj = Adj, Num = Num.Adj, SumNumNeigh = SumNumNeigh)
InitsCar = list(list(beta0 =0,beta1 =0, beta2 = 0, beta3 =0))

pcar.bugs = bugs(data=DataCar, 
                 model.file = "poisson-CAR_draft2.txt", 
                 parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "phi"), 
                 inits = rep(InitsCar,5), 
                 n.iter = 10000, 
                 n.burnin = 5000, 
                 n.thin = 1, 
                 n.chains = 5)

par(mfrow = c(2,2))
plot(pcar.bugs$sims.list$beta0, type ="l", ylab = "Samples beta0")
plot(pcar.bugs$sims.list$beta1, type ="l", ylab = "Samples beta1")
plot(pcar.bugs$sims.list$beta2, type ="l", ylab = "Samples beta2")
plot(pcar.bugs$sims.list$beta3, type ="l", ylab = "Samples beta3")
plot(pcar.bugs$sims.list$phi[,1], type ="l", ylab = "Samples phi")
plot(pcar.bugs$sims.list$phi[,50], type ="l", ylab = "Samples phi")
plot(pcar.bugs$sims.list$phi[,100], type ="l", ylab = "Samples phi")
plot(pcar.bugs$sims.list$phi[,150], type ="l", ylab = "Samples phi")
plot(pcar.bugs$sims.list$phi[,300], type ="l", ylab = "Samples phi")
plot(pcar.bugs$sims.list$phi[,600], type ="l", ylab = "Samples phi")

print(pcar.bugs, digits.summary = 3)
pcar.bugs$summary[2,1:7]

## 3.b - Use Gelman-Rubin diagostic (Rhat) and geweke diagnostic plots to assess convergence.

#### Gelman-rubin, checking RHAT
par(mfrow = c(1,1))
hist(pcar.bugs$summary[,8],main = " Distribution of Rhat", ylab = " Frequency", xlab = "Rhat",col= c("slategray1", "turquoise1","lightblue", "lavender", "lightcyan1", "cyan"))
### it is between 1 -1.2 --> converge.

### geweke convergence
par(mfrow = c(2,2))
geweke.plot(as.mcmc(as.data.frame(pcar.bugs$sims.list$beta0)))
geweke.plot(as.mcmc(as.data.frame(pcar.bugs$sims.list$beta1)))
geweke.plot(as.mcmc(as.data.frame(pcar.bugs$sims.list$beta2)))
geweke.plot(as.mcmc(as.data.frame(pcar.bugs$sims.list$beta3)))
geweke.plot(as.mcmc(as.data.frame(pcar.bugs$sims.list$phi[,c(1,50,150,300)])))


###3c: Calculate pearson residuals and check 2 model assumption 
pcar.bugs$summary[,1]
phi.mean = apply(pcar.bugs$sims.list$phi, 2, mean)
phi.mean

#########
mu.car = london3s$Expected*exp(-3.492363e-01 -2.529970e-02*as.numeric(scale(london3s$PM25)) +2.052273e-01*as.numeric(scale(london3s$JSA))-7.061127e-02*as.numeric(scale( london3s$Price))+phi.mean)
pearson.car = (london3s$Observed - mu.car)/sqrt(mu.car)
var(pearson.car)

## check overdispersion - plot fitted values and residual
par(mfrow = c(1,1))
plot(pearson.car ~ mu.car, xlab = "Fitted values", ylab = "Pearson residual", main = "Pearson residuals vs mu", col = "hotpink3")
var(pearson.car)

## check independence - moran's i statistics
Wnb.car <- nb2listw(poly2nb(london3s))
moran.car = moran.test(pearson.car, Wnb.car)
moran.car
moran.car$estimate[1]

#############################4 Summary table ################################
names(preg.bugs)
names(pcar.bugs)


Model = c("Poisson Regression", "Poisson Car")
DIC = c(preg.bugs$DIC,pcar.bugs$DIC) 
MoranI = c(moran$estimate[1],moran.car$estimate[1])
p_value = c(moran$p.value,moran.car$p.value)
PS_Var = c(var(pearson), var(pearson.car))



PM25_mean_exp = c(exp(preg.bugs$summary[2,1]),exp(pcar.bugs$summary[2,1]))
PM25_95Low_exp = c(exp(preg.bugs$summary[2,3]),exp(pcar.bugs$summary[2,3]))
PM25_95_HI_exp = c(exp(preg.bugs$summary[2,7]),exp(pcar.bugs$summary[2,7]))
PM25_mean = c(preg.bugs$summary[2,1],pcar.bugs$summary[2,1])
PM25_95Low = c(preg.bugs$summary[2,3],pcar.bugs$summary[2,3])
PM25_95_HI = c(preg.bugs$summary[2,7],pcar.bugs$summary[2,7])



summary.table.converted = data.frame(Model, DIC, MoranI, p_value, PS_Var, PM25_mean_exp, PM25_95Low_exp, PM25_95_HI_exp)
summary.table.converted

summary.table = data.frame(Model, DIC, MoranI, p_value, PS_Var, PM25_mean, PM25_95Low, PM25_95_HI)
summary.table

PM25_mean
PM25_95_HI
