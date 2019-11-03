##Bernoulli-beta example
#Estimation
a0 = 4
b0 = 4
n = 20
xbar = 0.4

posterior <- function(x) { 
  return( pbeta(x,a0+n*xbar,n-n*xbar+b0) )
}

prior <- function(x) { 
  return( pbeta(x,a0,b0) )
}
##################################3
#Hypothesis testing
#Discritized version
null = 0.5

RB <- function(x,eps) { 
  a = abs(posterior(x+eps)-posterior(x-eps))
  b = abs(prior(x+eps)-prior(x-eps))
  return( a/b  )
}
m = 1000
x = seq(0.001,0.999,1/m)
eps = 0.05
plot(x,RB(x,eps))

#Direction
a = RB(null,eps) #1.421018

#Strength
solver3 <- function(eps) {
  ans = NULL
  for (i in 1:length(x)) {
    if (RB(x[i],eps)>a) {
      ans = c(ans,x[i])
    }
    else {
    }
  }
  return(ans)
}
result = solver3(eps)

total = sum(RB(x,eps)) #Total area
greater = sum(RB(result,eps)) #Greater than
str = (total-greater)/total #0.3541251
###############################################
##Estimation
#Solve the argsup(RB)
solver <- function(eps) {
  x = seq(0.001,0.9999,1/m)
  ans = NULL
  for (i in 1:length(x)) {
      ans = c(ans,RB(x[i],eps))
  }
  result = x[which.max(ans)]
  
  return(result)
}
solver(eps) #0.4

#Solve for the plausible region
solver2 <- function(eps) {
  x = seq(0.001,0.9999,1/m)
  ans = NULL
  for (i in 1:length(x)) {
    if (RB(x[i],eps)>1) {
      ans = c(ans,x[i])
    }
    else {
    }
  }
  return(ans)
}
result = solver2(eps=10^-6)
RB(result[1])
RB(result[length(result)])
#The plausible region is
result[1] #0.273
result[length(result)] #0.536
result[length(result)]-result[1]  #0.263
#########################################
##Bias for Hypothesis testing
RB2 <- function(c0,t,n,a,b) { #c = lambda
  return(dbinom(t,n,c0)/mt(t,n,a,b))
}

rb2 = NULL
n = 20
t = seq(0,n,1)
for (i in 1:length(t)) {
  rb2 = c(rb2, RB2(c0=0.45,t[i],n=20,a=4,b=4))
}
plot(t,rb2)

#bias against
bias1 = sum(rb2[rb2<=1])/sum(rb2) #0.3269325

#bias in favour of theta = (0.45,0.55)
bias2 = sum(rb2[rb2>1])/sum(rb2) #0.7626085
bias3 = sum(rb2[rb2>1])/sum(rb2) #0.7626085
####################################################################
#Prior-data conflict 
mt <- function(t,n,a,b) { #Prior predictive distribution
  result = choose(n,t)*gamma(a+b)*gamma(t+a)*gamma(n-t+b)/gamma(a)/gamma(b)/gamma(n+a+b)
  return(result)
}
#Example
prob = NULL
n = 10
t = seq(0,n,1)
for (i in 1:length(t)) {
  prob = c(prob, mt(t[i],10,5,20))
}
plot(t,prob)

#p-value prior-data conflict (t0 = 9)
mt0 = mt(9,n,5,20) #0.0001090536
pvalue = sum(prob[prob<=mt0]) #0.0001166874