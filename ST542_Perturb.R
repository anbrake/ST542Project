#### Load Libraries ####
library(rjags)
library(ggplot2)
library(gridExtra)
library(grid)

#### Generate Data ####
set.seed(10) #set seed for repeatability
x = matrix(rnorm(100*3),100,3); epsi = rnorm(100) #generate x data and noise

betas = sample(-10:10, 4, replace = TRUE)
Y = cbind(1, x) %*% betas + epsi #generate the Y vector
dta = data.frame(y=Y, x=x)

#### Functions ####
R2 = function(Y, X, beta, pert= 0){
  beta.new = beta + pert
  Y.hat = X %*% beta.new
  r2 = Y - Y.hat
  sum(r2^2)
}

#### Beta Estimates ####
# OLS #
fit = lm(y ~ x, data = dta)
OLS.bhat = coef(fit)

# Bayes #
X = cbind(1,x)
p=ncol(X)
n=nrow(X)

modelData = list(X=as.matrix(X), Y=as.vector(Y), n=nrow(X), p=ncol(X))
model_string = textConnection("model {
  #Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] = inprod(X[i, ], beta[])
  }
  
  #Priors
  for(j in 1:p){
    beta[j] ~ dnorm(0, 1/100)
  }
  tau ~ dgamma(0.1, 0.1)
  
  sigma = 1/sqrt(tau)
}")
inits = list(beta = rep(0, p), tau=1)
jags_model = jags.model(model_string, data=modelData, inits=inits, n.chains = 2, quiet = T)
update(jags_model, 1000, progress.bar="none")

BayesianBetas = coda.samples(jags_model, variable.names = c("beta", "sigma"), n.iter = 1000, progress.bar="none")
Bayes.bhat = colMeans(BayesianBetas[[1]])[1:4]

#### Perturb Set-Up ####
total.perturb = seq(-4,4,length = 20)

#### Equal Perturbations ####
eq.pert = total.perturb / 4
eqRSS.OLS = NA
eqRSS.Bayes = NA

#Calculate Resid
for(p in 1:20){
eqRSS.OLS[p]   = R2(Y,cbind(1,x), OLS.bhat, pert=eq.pert[p])
eqRSS.Bayes[p] = R2(Y,cbind(1,x), Bayes.bhat, pert=eq.pert[p])
}



#### Unequal Perturbations ####
un.pert = cbind(total.perturb/2, total.perturb/3, total.perturb/4, total.perturb/5)
tot.unPtb = rowSums(un.pert)

unRSS.OLS = NA
unRSS.Bayes = NA

#Calculate Resid
for(p in 1:20){
  unRSS.OLS[p]   = R2(Y,cbind(1,x), OLS.bhat, pert=un.pert[p,])
  unRSS.Bayes[p] = R2(Y,cbind(1,x), Bayes.bhat, pert=un.pert[p,])
}

#### Random Perturbations ####
set.seed(1)
rand.pertPos = matrix(runif(10*4,0,3),10,4)
set.seed(2)
rand.pertNeg = matrix(runif(10*4,-3,0),10,4)
rand.pert = rbind(rand.pertNeg, 0,rand.pertPos)
tot.randPtb = rowSums(rand.pert)


randRSS.OLS = NA
randRSS.Bayes = NA

#Calculate Resid
for(p in 1:21){
  randRSS.OLS[p]   = R2(Y,cbind(1,x), OLS.bhat, pert=rand.pert[p,])
  randRSS.Bayes[p] = R2(Y,cbind(1,x), Bayes.bhat, pert=rand.pert[p,])
}


#### Combine Data into Dataframes ####
eq = data.frame(RSS.o = eqRSS.OLS, RSS.b = eqRSS.Bayes, ptb = total.perturb)
un = data.frame(RSS.o = unRSS.OLS, RSS.b = unRSS.Bayes, ptb = tot.unPtb)
rand = data.frame(RSS.o = randRSS.OLS, RSS.b = randRSS.Bayes, ptb = tot.randPtb)
#### Make Plots ####
par(mfrow=c(1,2))
plot(total.perturb, eqRSS.OLS, xlim = c(-5,5))
plot(total.perturb, eqRSS.Bayes)

plot(tot.randPtb, randRSS.OLS)
grid()
plot(tot.randPtb,randRSS.Bayes)
grid()

#### Equal PTB PLots ####
eq.OLS.plt = ggplot(data = eq, aes(x=ptb, y=RSS.o))+
  geom_smooth(size = 0.5)+
  geom_point(size = 2)+
  ggtitle('OLS Estimates')+
  ylab("RSS")+
  xlab("Perturbations")

eq.Bayes.plt = ggplot(data = eq, aes(x=ptb, y=RSS.b))+
  geom_smooth(size = 0.5)+
  geom_point(size = 2)+
  ggtitle('Bayes Estimates')+
  ylab("RSS")+
  xlab("Perturbations")

#Combine plots
grid.arrange(eq.OLS.plt, eq.Bayes.plt, ncol=2,
    top=textGrob('Equal Perturbation RSS Plots',gp=gpar(fontsize=20,font=3)))
#### Unequal PTB PLots ####
un.OLS.plt = ggplot(data = un, aes(x=ptb, y=RSS.o))+
  geom_smooth(size = 0.5)+
  geom_point(size = 2)+
  ggtitle('OLS Estimates')+
  ylab("RSS")+
  xlab("Perturbations")+
  expand_limits(x=c(-6,6))

un.Bayes.plt = ggplot(data = un, aes(x=ptb, y=RSS.b))+
  geom_smooth(size = 0.5)+
  geom_point(size = 2)+
  ggtitle('Bayes Estimates')+
  ylab("RSS")+
  xlab("Perturbations")+
  expand_limits(x=c(-6,6))

#Combine 
grid.arrange(un.OLS.plt, un.Bayes.plt, ncol=2,
  top=textGrob('Unequal Perturbation RSS Plots',gp=gpar(fontsize=20,font=3)))

#### Random Unif PTB PLots ####
rand.OLS.plt = ggplot(data = rand, aes(x=ptb, y=RSS.o))+
  geom_smooth(se=F,size = 0.5)+
  geom_point(size = 2)+
  ggtitle('OLS Estimates')+
  ylab("RSS")+
  xlab("Perturbations")

rand.Bayes.plt = ggplot(data = rand, aes(x=ptb, y=RSS.b))+
  geom_smooth(se=F,size = 0.5)+
  geom_point(size = 2)+
  ggtitle('Bayes Estimates')+
  ylab("RSS")+
  xlab("Perturbations")
#Combine 
grid.arrange(rand.OLS.plt, rand.Bayes.plt, ncol=2,
top=textGrob('Random Perturbation RSS Plots',gp=gpar(fontsize=20,font=3)))
