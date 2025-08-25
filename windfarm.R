
## fit multiple dataset using CAR models

library(MASS)
library(fields)
library(akima)
library(maps)
library(maptools)
library(spdep)
library(classInt)
library(RColorBrewer)
library(SpatialEpi)
library(R2WinBUGS)
library(CARBayes)

df1 <- read.csv('result_sample1.csv')
df2 <- read.csv('result_sample2.csv')

Power1 <- 0.001*df1$Power
xx <- df1$x
yy <- df1$y
EFL1 <- 0.001*df1$EFL
FL1 <- 0.001*df1$FL

EFL2 <- 0.001*df2$EFL
FL2 <- 0.001*df2$FL


Coords <- cbind(xx,yy)

# data.knn.nb <- knearneigh(Coords,k=8)
data.nb <- cell2nb(5,5,type="queen")

data.WB.weights <- nb2WB(data.nb)
adj.data <- data.WB.weights$adj
num.data <- data.WB.weights$num
weights.data <- data.WB.weights$weights

N <- length(EFL1)

M <- as.integer(2)

rep.data <- rep(1:N,num.data)

W <- matrix(0,N,N)

for(i in 1:N){
  W[i,adj.data[rep.data==i]] <- rep(1,num.data[i])  
}

input_data <- cbind(EFL1,EFL2)

inits.data <- function(){
  list(beta0=runif(1,6000,10000),beta1=runif(1,0,1000), invtau2=0.01, invsigma2=0.01,phi=c(rep(0,N)))
}

V <- cbind(rep(12,N),rep(10,N))

### Here we organize the data to pass to WinBUGS in a list 
### with all the names used in the model that we are running in WinBUGS
data.car <- list("Y"=input_data,"V"=V,"N"=N,"M"=M,"adj"=adj.data,
                  "weights"=weights.data,"num"=num.data)

mod.data.car <- bugs(data=data.car,inits=inits.data,
                 model.file="C:/Users/mingdyou/PhD/Computing/Gibbs/farm.bug.txt",
                 parameters.to.save=c("beta0","beta1","invtau2","tau2","invsigma2","sigma2","phi","mu","res"),n.chains=1,n.iter=50000,n.sims=25000,
                 bugs.directory="C:/Users/mingdyou/PhD/Computing/winbugs14",debug=FALSE)

plotvar <- EFL1

nclr <- 11
plotclr <- brewer.pal(nclr,"Blues")

# This command divides the variable to plot in nclr classes with fixed breaks. The break points are specified in
# the fixedBreaks. Other methods can be chose. See ?classIntervals
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))

colcode <- findColours(class,plotclr)

length.x <- 25
length.y <- 25

x.grid <- seq(0,1800,length=length.x)
y.grid <- seq(0,1400,length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)

# This first command does not really plot anything except for the legend
image.plot(x.grid,y.grid,z.grid,xlab="x (m)",ylab="y (m)",zlim=z.lim,
           col=plotclr,breaks=class$brks,main="EFL (kNm): mean wind speed = 10m/s")
points(xx,yy,col=colcode,pch=19,cex=2.5)

