#############################
## Monte carlo integration ##
#############################

## MC approximation of a uniform (0,1) distribution

nsim=500
u=runif(nsim,0,1)
mc=mean(u)
mc
mc=c(); v=c(); upper=c(); lower=c()
for(j in 1:nsim){mc[j]=mean(u[1:j])
  v[j]=(j^{-1})*var(u[1:j])
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])}
mc
plot(mc,type="l")

# Using ggplot2 to obtain good graphs
library(ggplot2)
values=c(mc,upper,lower)
type=c(rep("mc",nsim),rep("upper",nsim),rep("lower",nsim))
iter=rep(seq(1:nsim),3)
data=data.frame(val=values, tp=type,itr=iter)
data
ggplot(data, aes(itr,val))+geom_point(aes(colour=factor(tp)))

### Exercise 1. Approximate P(0<X<3), X~cauchy(0,1) using MC approximation.
### Result: X1~N(0,1), and X2~N(0,1) and X1 and X2 are independent
### then X1/X2 ~ Cauchy (0,1)
### Use ordinary MC approximation sum_{j=1}^{nsim} 1[0<X_i<3]/nsim

nsim=10000
x1=normal_sim(nsim)
x2=normal_sim(nsim)
x=x1/x2
hist(x)
length(x)
indicator=function(x){
  y=ifelse((x>0 & x<3),1,0)
  return(y)}
indicator(4)
indicator(2)
mc=c(); v=c(); upper=c(); lower=c()
for(j in 1:nsim){ mc[j]=mean(indicator(x[1:j]))
  v[j]= (j^{-1})*var(indicator(x[1:j]))
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])}

plot(mc,type="l")
values=c(mc,upper,lower)
type=c(rep("mc",nsim),rep("upper",nsim),rep("lower",nsim))
iter=rep(seq(1:nsim),3)
data=data.frame(val=values, tp=type,itr=iter)
#data
ggplot(data, aes(itr,val))+geom_point(aes(colour=factor(tp)))

## Exercise 2. Do the same question using importance sampling and
## using a uniform (0,3) instrument.

cdensity=function(x){y=1/(pi*(1+x^2));return(y)}
cdensity(1)
nsim=500
x=runif(nsim,0,3)
mc=mean(3*cdensity(x))
mc
mc=c(); v=c(); upper=c();lower=c()
for(j in 1:nsim){
  mc[j]=mean(3*cdensity(x[1:j]))
  v[j]=(j^{-1})*var(3*cdensity(x[1:j]))
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])}
plot(mc,type="l")