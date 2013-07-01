#################### README
# ML + LSDA Assignment 2 Question 3: Supernova Classification
# Zi Chong Kao, 4/18/2013
# 
# Start by running commands from the main section
# Then run all the functions in main
# Some of the functions have variables that are hard-coded from the main
#
# Runtime ~ 3 minutes
#
#################### FUNCTIONS

# from http://www.stat.cmu.edu/~cshalizi/350/lectures/26/lecture-26.pdf
# w0: initial guess for minimum w
# tolerance: minimum change in y values between loops before termination
my.newton = function(f,f.prime,f.prime2,w0,tolerance=1e-10,max.iter=200) {
  w = w0
  old.f = f(w)
  iterations = 0
  made.changes = TRUE
  while(made.changes & (iterations < max.iter)) {
    iterations <- iterations +1
    made.changes <- FALSE
    new.w = w - solve(f.prime2(w)) %*% f.prime(w)
    new.f = f(new.w)
    relative.change = abs(new.f - old.f)/old.f -1
    made.changes = (relative.change > tolerance)
    w = new.w
    old.f = new.f
  }
  if (made.changes) {
    warning("Newton's method terminated before convergence")
  }
  return(list(minimum=w,value=f(w),deriv=f.prime(w),deriv2=f.prime2(w),
  iterations=iterations,converged=!made.changes))
}

# splits data 60% 40% at random
test.index = function(dataset) {
  nrows = dim(dataset)[1]
  ntest = round(nrows/4)
  testindex = sample(1:nrows,ntest)
  return(testindex)
} 

is.error = function(data,b_hat,true_val){
  fitted = data %*% b_hat
  # fitted more than zero means logistics variable more than 1/2
  more_than_zero = fitted>0
  return(more_than_zero!=true_val)
}
  
se = function(vector){
  return(sd(vector)/sqrt(length(vector)))
}

logistics = function(X,beta){
  exponent = X %*% beta
  return(as.vector(1/(1+exp(-exponent))))
}

likelihood = function(beta){
  sum = 0
  for (i in (1:length(Y))){
    Xb = X[i,] %*% beta 
    sum = sum + ((Y[i] * Xb) - log(1+exp(Xb)))
  }
  return(-sum)
}

gradient = function(beta){
  p = logistics(X,beta)
  g = -t(X) %*% (Y - p) + 2*lambda*beta
  #print(g)
  return(g)
}

hessian = function(beta){
  I = diag(ncol(X))
  p = logistics(X,beta)
  W = diag((p)*(1-p),ncol=length(p))
  #print(lambda*I)
  H = t(X) %*% W %*% X + 2*lambda*I
  #print(H)
  return(H)
}

########################### ERROR BAR
# from: http://monkeysuncle.stanford.edu/?p=485

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
  stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


########################### MAIN 

freqdata <- read.table("C:/Users/kao/Google Drive/Desktop/supernovae.dct.dat",header=TRUE, sep=" ")
dataset=freqdata

# set lambda values as first column of main_storage
main_storage = cbind((0:50)/10-2.5,rep(0,50),rep(0,50),rep(0,50),rep(0,50))

# loop through lambda values
for (i in 1:length(storage[,1])) {
lambda = exp(storage[i,1])

error_storage_train = rep(1:50)
error_storage_test  = rep(1:50)

# for each lambda split data 50 times
for (j in 1:50) {
test_index = test.index(dataset)
test_data  = dataset[test_index,]
train_data = dataset[-test_index,]

Y = train_data[,1] 
X = data.matrix(cbind(I = rep(1,length(train_data[,1])),train_data[,-1]))
Y_test = test_data[,1]
X_test = data.matrix(cbind(I = rep(1,length(test_data[,1])),test_data[,-1]))

# initial = glm(y~.,family=binomial(logit), data=train_data)
# init_coef = data.matrix(coef(initial))
# the glm algorithm wasn't converging. 
# since min_b were eventually in the E-2 to E -7 range, I set them to 0 first.
init_coef = rep(0,101)

# now define logistics, likelihood, gradient and hessian functions

final = my.newton(likelihood,gradient,hessian,init_coef)
min_b = as.numeric(final$minimum)

train_err_rec = is.error(X,min_b,Y)
test_err_rec  = is.error(X_test,min_b,Y_test)
error_storage_train[j] = mean(train_err_rec)
error_storage_test[j] = mean(test_err_rec)
}

# recover mean and se of errors of each lambda
main_storage[i,2] = mean(error_storage_train)
main_storage[i,3] = se(error_storage_train)
main_storage[i,4] = mean(error_storage_test)
main_storage[i,5] = se(error_storage_test)

}

######################### PLOTTING

# training set
plot(main_storage[,1],main_storage[,2],main="Freq Domain",ylab="Errors",xlab="log(lambda)",type="o",col="blue3")
error.bar(main_storage[,1],main_storage[,2],main_storage[,3],col="darkblue")

# testing set
lines(main_storage[,1],main_storage[,4],col="darkorange2",type="o")
error.bar(main_storage[,1],main_storage[,4],main_storage[,5],col="darkorange3")

