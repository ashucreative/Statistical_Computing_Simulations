#########################################################
# Marsaglia's polar method to simulate std normal R.V's #
#########################################################

# Nsim should be much smaller than 10^5

Nsim = 10^3

normal_sim = function(Nsim){
x1 = c(); x2 = c(); its = 0

for(j in 1:10^5){ # j is generating two unifrom r.v.'s at a time
u1 = runif(1, -1, 1); u2 = runif(1,-1,1)
s = u1^2 + u2^2
if (s<1){x1[j] = sqrt(-2*log(s))*(u1/sqrt(s))
    x2[j] = sqrt(-2*log(s))*(u2/sqrt(s)) 
    its = its + 1}
if (its >= (Nsim/2)){break}}

x1r = x1[is.na(x1)!= T]; x2r =x2[is.na(x2)!=T]
x=c(x1r, x2r)
return(x)
}

x = normal_sim(400)
hist(x)

z = mu+sigma*x
hist(z)

######################################
## Multivariate normal distribution ##
######################################

Nsim = 10^3; p =10; mu=rep(0,p)
rho = 0.5; Sig = rho^toeplitz(0:(p-1))

multi_normal = function(Nsim, mu,Sig){
  p = dim(Sig)[1]
  x = matrix(normal_sim(Nsim*p), nrow = Nsim, ncol = p)
  Gamma=chol(Sig)
  z = mu + t(Gamma)%*%t(x) # Here x is nxp matrix
  return(z)}

z = multi_normal(3,rep(1,p), Sig = Sig)  
z
dim(z)

##################################################################
## Truncated normal from shifted exponential using AR algorithm ##
##################################################################

d=1
f = function(x){exp(-x^2/2)} # Here we dont need to know the normalizing constant.
f(1)

# Now we need to know the instrumental density
g = function(x){lambda*exp(-lambda*(x-d))} #density of shifted exponential
lambda = 2; g(1)
lambda = (d+sqrt((d^2)+4))/2
lambda

M = (lambda^(-1))*(exp(lambda^2/2-lambda*d))
M

bound = function(f, g, M, x){f(x)/(M*g(x))}
bound(f, g, M, 1.5)

####################################
## Simulating from the instrument ##
####################################

shifted_exp = function(Nsim, lambda, d){ #this is the realization of shifted exponential
  u = runif(Nsim)
  e = -log(u)/lambda
  se = e+d
  return(se)}
shifted_exp(10, lambda = 2, d=1)

Nsim = 10^3
y = c(); its = 0
for(j in 1: 10^6){
  #step1
  x = shifted_exp(1, lambda = lambda, d = d) # Simulate from instrument
  u = runif(1)
  ## step 2
  if(u<=bound(f, g, M, x)){y[j]=x; its = its+1}
  if (its == Nsim){break}}  

yr = y[is.na(y)!=T]
hist(yr)


