################# Lasso vs Ridge Constraint ############

model <- "ridge"
nx <- 200; ny <- 200; x0 <- -10; x1 <- 10; y0 <- -10; y1 <- 10
a <- 5; b <- 2.5; rho <- -.5
x <- seq(x0, x1, len = nx)
y <- seq(y0, y1, len = ny)
z <- outer((x-a)^2,(y-b)^2,"+")-2*rho*outer((x-a),(y-b),"*"); z <- (-.5*(z)/(1-rho))
if (model=="lasso") z2 <- outer(abs(x),abs(y),"+") else z2 <- outer(x^2/3,y^2/3,"+")+1.25
rgl::persp3d(x=x,y=y,z=z,col=2)
rgl::persp3d(x=x,y=y,z=z2-2.5,add=T,col=4);
rgl::persp3d(x=x,y=y,z=rep(0,nx*ny),col=3,add=T)
for (i in 1:6*1)
  rgl::persp3d(x=x,y=y,z=z+i,col=2,add=T)

################# Functional Data Analysis #############
t <- 1:100/100;# t <- -50:50; 
p <- 11
X <- rep(1,length(t))
for (i in 1:p) X <- cbind(X,t^i)
kappa(t(X)%*%X,exact = T)

XtX.invr <- solve(t(X)%*%X + 0.0001*diag(11))
ts.plot(X,col=1:(p+1))

# ts.plot(X,col=1:(p+1),log="y")

X <- splines::bs(t, degree = 3,intercept = T, df=p)
ts.plot(X, col=1:p, type="t")
image(X)
image(t(X)%*%X)
image(solve(t(X)%*%X))

#####################################################

S <- X%*%solve(t(X)%*%X)%*%t(X)
sum(diag(S))
