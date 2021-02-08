## Approximate the value of Sqrt(2) using Newton Raphson update.

## Solve x2-2=0
g=function(x){x^2-2}
gprime=function(x){2*x}

# Initialize X value, tolerance, maximum iteration, and error
xold=1;tol=1e-5
maxits=3000;its=1;err=1

while(err>tol & its<maxits){
  xnew=xold-(g(xold)/gprime(xold))
  err=abs(xnew-xold);its=its+1;xold=xnew}
for(j in 1:maxits){
  xnew=xold-(g(xold)/gprime(xold))
  err=abs(xnew-xold);its=its+1;xold=xnew
  if(err<tol){break}}
xnew

# Exercise 2: Newton Raphson for logistic regression
## Generate data
set.seed(123)
n=300; x=rnorm(n); b=c(0.25,0.75)
X=cbind(1,x)
invlogit=function(X,b){p=exp(X%*%b)/(1+exp(X%*%b))
return(p)}
p=invlogit(X,b)
y=c(); for(i in 1:n){y[i]=rbinom(1,1,p[i])}
y

## Now we have a dataset (x,y) where X-->predictors, y-->binary responses
## use (x,y) to estimate the b coefficients for any matrix H: H^-1=solve(H)

p=dim(X)[2] ## Measure the dimension of the model
old_hb=rep(0,p) ## Arbitrary starting point
n=length(y)
V=matrix(0,n,n)
tol=1e-5
maxits=5000
its=0; err=10
while(err>tol & its<maxits){
  old_hp=invlogit(X,old_hb)
  diag(V)=old_hp*(1-old_hp)
  H=-t(X)%*%V%*%X; Hinv=solve(H)
  S=t(X)%*%(y-old_hp)
  new_hb=old_hb-Hinv%*%S
  err=max(abs(new_hb-old_hb))
  its=its+1
  old_hb=new_hb
}
new_hb

for(j in 1:100){
  old_hp=invlogit(X,old_hb)
  diag(V)=old_hp*(1-old_hp)
  H=-t(X)%*%V%*%X; Hinv=solve(H)
  S=t(X)%*%(y-old_hp)
  new_hb=old_hb-Hinv%*%S
  err=max(abs(new_hb-old_hb))
  its=its+1
  old_hb=new_hb
}

new_hb
