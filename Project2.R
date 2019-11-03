library(ggplot2)
data =  read.table('Data.txt', header = FALSE)
dat = data$V2 #Height

##1. Problem definition, ingredients, method of inference
##Plot the data (heights)
ggplot(data, aes(x=V1, y=V2)) + 
  geom_point(size=3, shape=16) +
  xlab("Time (in sec)") + 
  ylab("Beer foam height (in cm)") + 
  ggtitle("Decay of Beer Foam after pouring into glass")
################################################################
##2. (Sampling ) Model checking
#1). Data generating/measurement process - Decay
#- Exponential distribution is infinite divisible!

#2). Using the one-sample Kolmogorov-Smirnov test to assess if the data can be modelled by 
#an Exponential distribution
rate = 1/mean(dat) #Maximum likelihood estimates
ks.test(dat, "pexp", rate) #p-value = 0.1842

#3). Using graphical assessment to see if the data can be modelled by an Exponential distribution
x = rexp(n=100,rate=rate) #Exponential distribution, simulate 100 data

ggplot() +
  geom_density(aes(dat, fill='Data'),alpha=0.8) +
  geom_density(aes(x, fill='Model'),alpha=0.8) +
  xlab("Beer foam height (in cm)") + 
  ylab("Density") + 
  ggtitle("Kernel Density Plot - Data vs Model
  Exponential Decay of Beer Foam") +
  theme(legend.position = 'right') +
  scale_color_manual(values = c('Data' = 'orange', 'Model' = 'blue'))

#Plot the density function of the Exponential distribution
x = seq(1,100,1)
plot(x,dexp(x,rate=rate))
plot(seq(1,n,1),dat)
################################################################
##3. Prior elicitation
#1). Conjugate prior - Gamma(a,b)

#2). Hyperparameter selection (a,b) - Using past data(belief) of Half-life of beer foam decay
#Determine the value of (a,b) 
f <- function(rate0,rate1) {
  Fx = 0
  b = 0.01
  while (Fx <= 0.99) {
    a =  b*rate0 + 1 #match the mode
    Fx = pgamma(rate1,shape=a,rate=b)
    b = b + 0.01
  }
  param = c(a,b-0.01)
  return(param)
}

rate0 = log(2)/4 #0.1732868 = mode
rate1 = log(2)/2 #0.3465736 = extreme (1 in 100/99% percentile)
result = f(rate0,rate1) #a = 11.22392 b = 59.00000
a = result[1]
b = result[2]

#Check
#mode = (result[1]-1)/result[2] #0.1732868
#extreme = pgamma(rate1,result[1],result[2]) #0.9900032
#################################################################################
##4. Prior-data conflict
#1) Visual check (is posterior-prior difference surprising?)
#Plot the prior and posterior densities
n = length(dat) #sample size
tx = sum(dat) #sufficient stats

x = seq(0.001,1,1/1500)
fx = dgamma(x,shape=a,rate=b)
fy = dgamma(x,shape=a+n,rate=b+tx)

#Just the prior
ggplot() +
  geom_point(aes(x=x, y=fx),size=1, shape=16,color='Red') +
  xlab("Rate (lambda)") + 
  ylab("Density") + 
  ggtitle("Selected Prior distribution") 

#Prior vs Posterior
ggplot() +
  geom_point(aes(x=x, y=fx, color='Prior'),size=1, shape=16,alpha=0.5) +
  geom_point(aes(x=x, y=fy, color='Posterior'),size=1, shape=16) +
  xlab("Rate (lambda)") + 
  ylab("Density") + 
  ggtitle("Prior vs Posterior densities") +
  theme(legend.position = 'right') +
  scale_color_manual(values = c('Prior' = 'Red', 'Posterior' = 'Blue'))

#2) Use the prior-predictive distribution
mt <- function(t,n,a,b) { #Prior predictive distribution
  result = (b^a)*(t^(n-1))*gamma(a+n)/gamma(a)/gamma(n)/((t+b)^(a+n))
  return(result)
}

#Plot the prior-predictive distribution for different values of t
t = seq(0,max(dat)*n,0.1)
mT = mt(t,n,a,b)
t0 = sum(dat) #observed t0 

ggplot() +
  geom_point(aes(x=t, y=mT),size=1, shape=16,colour='Orange') +
  xlab("t") + 
  ylab("Density") + 
  ggtitle("Prior-predictive distribution") +
  geom_vline(xintercept = t0,color = "blue", size = 1, alpha = 0.5)

#p-value prior-data conflict
mt0 = mt(tx,n,a,b) #observed mT(t0) = 0.004491333
prob = NULL
for (i in 1:length(t)) {
  prob = c(prob, mt(t[i],n,a,b))
}
#plot(t,prob)
pvalue = sum(prob[prob<=mt0])/sum(prob) #0.1499355
################################################################################
##5. Inference
#Relative belief (discritized version)
posterior <- function(x) { #Gamma(a+n,b+tx)
  return( pgamma(x,shape=a+n,rate=b+tx) )
}

prior <- function(x) { #Gamma(a,b)
  return( pbeta(x,a,b) )
}

RB <- function(x,eps) { 
  a = abs(posterior(x+eps/2)-posterior(x-eps/2))
  b = abs(prior(x+eps/2)-prior(x-eps/2))
  return( a/b  )
}

#Plot the relative belief ratio RB(lambda|x)
m = 5000
x = seq(0.001,0.499,1/m)
eps = 0.1 #0.01,0.02,0.03,0.04,0.06.0.07,0.08,0.09.0.1
rb = RB(x,eps)
ggplot() +
  geom_point(aes(x=x, y=rb),size=1, shape=16,colour='Purple') +
  xlab("Rate (lambda)") + 
  ylab("RB") + 
  ggtitle("Plot of relative belief ratio")

##1) Estimation
#Solve the argsup(RB) [using discritized version, if cts version, estimate = ML estimate!]
solver <- function(eps) {
  ans = NULL
  for (i in 1:length(x)) {
    ans = c(ans,RB(x[i],eps))
  }
  result = x[which.max(ans)]
  
  return(result)
}
solver(eps) #0.1244

#Solve for the plausible region
solver2 <- function(eps) {
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
result = solver2(eps)

#The plausible region is (0.084,0.169)
lower = result[1] #0.084
upper = result[length(result)] #0.169
length = upper - lower #0.085 (reasonable uncertainty)
#####################################################
##2) Hypothesis testing
##HT1: H0: lambda = prior mode (0.173)
eps = 0.01
#Direction
null1 = rate0
rb1 = RB(null1,eps) #0.9313572 (against H0)
rb1

#Strength
solver3 <- function(eps,rb0) {
  ans = NULL
  for (i in 1:length(x)) {
    if (RB(x[i],eps)>rb0) {
      ans = c(ans,x[i])
    }
    else {
    }
  }
  return(ans)
}
result = solver3(eps,rb1)
total = sum(RB(x,eps)) #Total area
greater = sum(RB(result,eps)) #Greater than
str = (total-greater)/total #Strength = 0.3088631
str

##HT2: H0: lambda = prior extreme (0.347)
#Direction
null2 = rate1
rb2 = RB(null2,eps) #0.001979244 (against H0)

#Strength
result = solver3(eps,rb2)
total = sum(RB(x,eps)) #Total area
greater = sum(RB(result,eps)) #Greater than
str = (total-greater)/total #Strength = 0.0003155821
#####################################################33
#########################################
##Bias for Hypothesis testing #1
RB0 <- function(c0,t,n,a,b) { #c = lambda
  mt0 = dgamma(t,shape=n,rate=c0)
  mt = mt(t,n,a,b)
  return(mt0/mt)
}

rb2 = NULL
x0 = max(dat) #Initial height = 17.4
t = seq(0.001,n*x0,n*x0/1000)

bias1 <- function(c0) { #against
  for (i in 1:length(t)) {
    rb2 = c(rb2, RB0(c0,t=t[i],n,a,b))
  }  
  result = sum(rb2[rb2<=1])/sum(rb2)
  return(result)
}

bias2 <- function(c0) { #in favour
  for (i in 1:length(t)) {
    rb2 = c(rb2, RB0(c0,t=t[i],n,a,b))
  }  
  result = sum(rb2[rb2>1])/sum(rb2)
  return(result)
}

#bias against
b1 = bias1(null1) #0.3535029

#bias in favour (take the sup)
v1 = null1-0.01
v2 = null1+0.01
v = seq(v1,v2,(v2-v1)/100)

BIAS2 = NULL
for (i in 1:length(v)) {
  BIAS2 = c(BIAS2, bias2(v[i]))  
}
plot(BIAS2)
max(BIAS2)

