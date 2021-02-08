# In all exercises, I will be using simulation of unifor rvs as given

n = 100
x = runif(n,0,1) # simulate 100 realization of U(0,1)
x

################################
# Implementing inverse method  #
################################

# Simulate exponential
# 1. generate U(0,1)
#2. x = -(1/lambda)*log(U)

n = 10^4
U = runif(n,0,1)
E = -log(U) #Exp(1) distribution
hist(U)
hist(E)

# Bernoulli(p) using the inverse method
#1. generate U(0,1)
#2. x = =0 if U<1-p and x = 1 if U> 1-p

p = 0.25
U = runif(n,0,1)
X = ifelse(U<1-p, 0,1)
X

mean(X==1)

###############################
#### Transformation Method ####
###############################

#simulate gamma (alpha, beta) alpha --> integer
n = 10^4; alpha = 2; beta = 1
U = matrix(runif(n*alpha, 0 , 1), alpha, n) # to simulate 1 gamma, you need alpha uniforms
E = -log(U) #exp(1) # Negative log transformation
G = beta*colSums(E)
hist(G)

### Make a function to simulate gamma's
gamma_sim = function(Nsim, alpha,beta){
  if(floor(alpha) != alpha) stop("please supply integer alpha")
  U = matrix(runif(Nsim*alpha,0,1), alpha, Nsim)
  E = -log(U) #exp(1)
  G = beta*colSums(E)
  return(G)
}
gamma_sim (10,1.5,2)

hist(gamma_sim(Nsim = 1000, alpha = 1, beta = 0.2))

#simulate a beta distribution using a transformation of uniforms

##Y = sum_{i = 1}^{a}E_i/sum_{i = 1}^{a+b}E_i ~ beta(a,b)
Nsim = 10^4; a = 1; b = 1
beta_sim = function(Nsim,a,b){
  U = matrix(runif(Nsim*(a+b)), (a+b), Nsim)
E = -log(U)
A = matrix(E[(1:a),],nrow = a)
numerator = colSums(A)
denominator = colSums(E)
B = numerator/denominator
return(B)
}
beta_sim(Nsim = 10, a = 1, b=1)  


######################
#### AR algorithm ####
######################

# Simulate gamma (alpha,beta) from gamma(a,b)
alpha = 2.5; beta = 2

#instrumental density: gamma(a,b)
a = floor(alpha); b = a*beta/alpha
a;b

M = (beta^alpha)*((b)^(-a))*(((alpha-a)/(beta-b))^(alpha-a))*exp(a-alpha)
M

#compute the ratio f(x)/Mg(x) for any x

bound = function(M,x){  # Here we have stored alpha, beta thats why not adding values in fun definition
  z = dgamma(x,shape = alpha, scale = (1/beta))/(M*dgamma(x, shape = a, scale = (1/b)))
  return (z)
  }

bound (M, x = 2)
Nsim = 10
ysim= c()

######################
#implementation of AR#
######################

#Step 1: simulate from instrument and a U(0,1)
while(length(ysim)<Nsim){
x = gamma_sim(1,alpha = a, beta = b); u = runif(1)

#Step 2: accept if u<=bound(M,x)
y = ifelse(u<=bound(M=M, x=x), x, NA)
if(is.na(y)==F){ysim = c(ysim,y)}
}

y; ysim
