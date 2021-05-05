library(ggplot2)
library(dplyr)
library(ggpubr)

n=10
X = seq(-0.9, 1.1,len=n)[-10]
#Decided to go with no random noise
#set.seed(623); Y = -X^2 + rnorm(n,sd=0.1)

#Non-random 
Y = -X^2

# Problem that Y is negative, need to make positive and retain shape
Y.beta = Y + -min(Y)
plot(X, Y.beta)

Y = 5 * (Y.beta / sum(Y.beta))
plot(X,Y)

# #Drop 10th point to force even range about 0:
# X = X[-n]
# Y = Y[-n]
# Y.beta = Y.beta[-n]

# #technically some invented dist 
# test.fit = nls(Y~exp(-((X-mu)/sigma)^2), start=list(mu=0, sigma=0.5), trace = T,
#           control = list(maxiter = 100))
# plot(X,Y)
# lines(X, predict(test.fit), lty=2)

### Real X ###
##Works with better starting values
fit.normal = nls(Y~dnorm(X, mean, sd), start=c(mean=0, sd=0.25), trace=F,
          control = list(maxiter=100))
fit.cauchy = nls(Y~dcauchy(X, location, scale), start=c(location=0, scale=0.1), trace = F,
          control = list(maxiter=100))
fit.logistic = nls(Y~dlogis(X, location, scale), start=c(location=0, scale=0.1), trace = F,
          control = list(maxiter=1000))
lines(X, predict(fit.normal), col=1)
lines(X, predict(fit.cauchy), col=2)
lines(X, predict(fit.logistic), col=3)

#### Positive X #### 
X.pos = X+1 #Forces X positive
fit.gamma = nls(Y~dgamma(X.pos, shape, rate), start=c(shape=5, rate=5), trace = F,
                control = list(maxiter = 1000))
fit.exponential = nls(Y~dexp(X.pos, rate), start=c(rate=1), trace = F,
                control = list(maxiter=1000))
#Had to set initial DF=1
fit.chisq = nls(Y~dchisq(X.pos, df), start=c(df=2), trace=F, 
                control = list(maxiter = 1000))
#Don't think F is flexible enough to fit
#fit.f = nls(Y~df(X.pos, df1, df2), start=c(df1=10, df2=1), trace = T, control = list(maxiter=1000))
fit.logNormal = nls(Y~dlnorm(X.pos, meanlog, sdlog), start=c(meanlog=0, sdlog=1/5), trace = F,
                control = list(maxiter=100))
fit.weibull = nls(Y~dweibull(X.pos, shape, scale), start=c(shape=3, scale=1), trace = F,
                control=list(maxiter=100))

plot(X.pos,Y)
x.pos = seq(0,2,len=100)
lines(X.pos, predict(fit.gamma), col="red")
lines(X.pos, predict(fit.exponential), col="blue")
lines(X.pos, predict(fit.chisq), col="green")
lines(X.pos, predict(fit.normal), col="black")
lines(X.pos, predict(fit.logNormal), col="pink")
lines(X.pos, predict(fit.weibull), col="orange")

#### Make X e[0,1] #### 
X.0t1 = (X+1) / 2 
Y.beta = Y.beta + 1/2
#any var obtained from these fits should be multiplied by 4
#to get us on the same scale as the original X data 
fit.beta = nls(Y.beta~dbeta(X.0t1, shape1, shape2), start=c(shape1=1.7, shape2=1.7), trace=T,
               control=list(maxiter=1000))
plot(X.0t1, Y.beta, ylim=c(0,2))
lines(X.0t1, predict(fit.beta))




### Plotting ###
X.long = seq(-0.9, 1.1, len=100)[-c(91:100)]
Y.long = -X.long^2
Y.long = Y.long + -min(Y.long)

### Real X ###
realPlot = data.frame(
  X.long = X.long,
  Y.long = Y.long,
  normal = dnorm(X.long, mean=coef(fit.normal)[1], sd=coef(fit.normal)[2]),
  cauchy = dcauchy(X.long, location=coef(fit.cauchy)[1], scale=coef(fit.cauchy)[2]),
  logistic = dlogis(X.long, location=coef(fit.logistic)[1], scale=coef(fit.logistic)[2])
)
frame1 = ggplot(mapping=aes(x=X.long, y=Y.long))+
  geom_point(mapping=aes(x=X, y=Y))+
  geom_line(mapping=aes(y=realPlot$normal, color="Normal"))+
  geom_line(mapping=aes(y=realPlot$cauchy, color="Cauchy"))+
  geom_line(mapping=aes(y=realPlot$logistic, color="Logistic"))+
  labs(x="X", y="Y", color="Distributions")+
  ggtitle("Real-Valued Support distributions")

### Positive X ###
X.pos.long = X.long+1 
positivePlot = data.frame(
  X.pos.long = X.pos.long,
  Y.long = Y.long,
  gamma = dgamma(X.pos.long, shape=coef(fit.gamma)[1], rate=coef(fit.gamma)[2]),
  exponential = dexp(X.pos.long, rate=coef(fit.exponential)),
  chisq = dchisq(X.pos.long, df=coef(fit.chisq)),
  logNormal = dlnorm(X.pos.long, meanlog=coef(fit.logNormal)[1], sdlog=coef(fit.logNormal)[2]),
  weibull = dweibull(X.pos.long, shape=coef(fit.weibull)[1], scale=coef(fit.weibull)[2])
)
frame2 = ggplot(mapping=aes(x=X.pos.long, y=Y.long))+
  geom_point(mapping=aes(x=X.pos, y=Y))+
  geom_line(aes(y=positivePlot$gamma, color="Gamma"))+
  geom_line(aes(y=positivePlot$exponential, color="Exponential"))+
  geom_line(aes(y=positivePlot$chisq, color="Chi Squared"))+
  geom_line(aes(y=positivePlot$logNormal, color="Log-Normal"))+
  geom_line(aes(y=positivePlot$weibull, color="Weibull"))+
  labs(x="X", y="Y", color="Distributions")+
  ggtitle("Positive-Valued Support Distributions")

#### X e[0,1] #### 
X.long01 = (X.long+1) / 2 
Y.long.beta = Y.long + 1/2
zeroOnePlot = data.frame(
  X = X.long01,
  Y = Y.long.beta, 
  beta = dbeta(X.long01, shape1=coef(fit.beta)[1], shape2=coef(fit.beta)[2])
)
frame3 = ggplot(mapping=aes(x=X.long01, y=Y.long.beta))+
  geom_point(mapping=aes(x=X.0t1, y=Y.beta))+
  geom_line(mapping=aes(y=zeroOnePlot$beta, color="Beta"))+
  labs(x="X", y="Y", color="Distributions")+
  ggtitle("[0, 1] Support Distributions")

#Display plots separately
frame1
frame2
frame3

#Display together
ggarrange(frame1, frame2, frame3, nrow=3, ncol=1)

