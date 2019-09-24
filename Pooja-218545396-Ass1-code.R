#setwd("~/Multivariate and categorical data analysis/assignment1")

#Load the provided data 
#---------------------------------------------------------------
the.data <- as.matrix(read.csv("GBRData.csv", header = TRUE, sep =","))

#Generate a random set of 250 rows
#---------------------------------------------------------------
my.data <- the.data [sample(1:556,250),c(1:4)]

#Write to file
#---------------------------------------------------------------
write.table(my.data,"name-StudentID-GBRMyData.txt") 


#========================================================================================================================
#1.1) Generate histograms
#========================================================================================================================

#histograms for 'HeronIsland-Salinity'
#-------------------------------------
hist(my.data[,1], main="Histogram for HeronIsland-Salinity", xlab="Salinity in practical salinity units", col="green") 

#histograms for 'HeronIsland-Water Temperature'
#----------------------------------------------
hist(my.data[,2], main="Histogram for HeronIsland-WaterTemperature", xlab="WaterTemperature in degree celcius", col="blue")

#========================================================================================================================
#1.2) Parallel Box plot using the two variables; 'HeronIsland-WaterTemperature' and the 'LadyMusgrave-Water Temperature'. 
#========================================================================================================================
WaterTemperature <- c( "LadyMusgrave","HeronIsland") 
boxplot(my.data[,4],my.data[,2], names=WaterTemperature, horizontal = TRUE, main="WaterTemperatures of HeronIsland and LadyMusgrave", xlab="Temperature in degrees celcius", col="Red") 

#5 number summary measures for HeronIsland-WaterTemperature
#--------------------------------------------------------------
summary(my.data[,2])

#5 number summary measures for LadyMusgrave-Water Temperature
#--------------------------------------------------------------
summary(my.data[,4])

#========================================================================================================================
#1.3) Scatter plot of 'HeronIsland-Water Temperature' and 'LadyMusgrave-Water Temperature' for 150 points
#========================================================================================================================
HeronIsland <- my.data[1:150,2]
LadyMusgrave <- my.data[1:150,4]

# Plot with main and axis titles
# Change point shape (pch = 19) and remove frame.
#------------------------------------------------
plot(HeronIsland, LadyMusgrave, main = "Scatter plot of water temperatures between HeronIsland and LadyMusgrave",
     xlab = "HeronIsland - Water Temperature in degrees celcius", ylab = "LadyMusgrave - Water Temperature in degrees celcius",
     pch = 19)
abline(lm(LadyMusgrave ~ HeronIsland, data = mtcars), col = "blue")

#Fit Linear regression model
#------------------------------------------------
lm.wt <- lm(LadyMusgrave ~ HeronIsland)
summary(lm.wt)

#Calculate correlation coefficient
#------------------------------------------------
cor(LadyMusgrave,HeronIsland)

#Calculate coefficient of determination
#------------------------------------------------
summary(lm.wt)$r.squared 

#========================================================================================================================
#5.4) prior, likelihood and the posterior distributions of snowfall received at Bright city
#========================================================================================================================

colors <- c("black", "blue", "red")
labels <- c("prior (mean=10, var=1)", "likelihood (x10=15, var=25)", "posterior")
#prior
mean=10; sd=sqrt(1)
x <- seq(-10,10,length=200)*sd + mean
hx <- dnorm(x,mean,sd)
plot(x, hx, type="n", xlab="", ylab="", ylim=c(0, 0.8),xlim=c(-1, 20), main="Bayesian estimation", axes=TRUE)
lines(x, hx, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors[1])
#liklihood
mean1=15; sd1=5
hx <- dnorm(x,mean1,sd1)
lines(x, hx,lwd=2, col=colors[2])
legend("topleft", inset=.005,
 labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
#posterior
mean2=11.4285; sd2=sqrt(0.7142)
hx <- dnorm(x,mean2,sd2)
lines(x, hx,lwd=2, col=colors[3])
legend("topleft", inset=.005,
 labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
 
#========================================================================================================================
#6.1) K-Means clustering
#========================================================================================================================
 
zz<-read.table("SITEdata2019.txt")
zz<-as.matrix(zz)

# Scatter Plot of the provided data
#------------------------------------------------
plot(zz[,1], zz[,2], main = "Scatter plot of data", pch = 19) 

# K-mean clustering with 5 classes
#------------------------------------------------
km.out =kmeans (zz, 5, nstart =20)
km.out$cluster
plot(zz, col =(km.out$cluster) , main="K-Means Clustering Results with K=5", xlab ="", ylab="", pch =20, cex =2)

# Vary K and calculate Total within sum of squares
#------------------------------------------------
wss <- function(k) {
  kmeans(zz, k, nstart = 10 )$tot.withinss
}

m=20
mymat <- matrix(NA, nrow=19,ncol =2)

for (k in 2:20){
   x<- wss(k)
  mymat[k-1,]<- c(k,x)
 } 
 
plot(mymat[,1], mymat[,2], type="b", pch = 19, main="Plot of K versus TOTWSS", xlab="Number of clusters (K)", ylab="Total within sum of squares")

#========================================================================================================================
# 6.2) Spectral clustering
#========================================================================================================================
 
zz<-read.table("SITEdata2019.txt")
zz<-as.matrix(zz)

#compute similarity matrix
#------------------------------------------------
maxMinScale <- function(x){(x-min(x))/(max(x)-min(x))}
zz<-maxMinScale(zz)
plot(zz)

#Compute affinity matrix
#------------------------------------------------
dZZ<-as.matrix(dist(zz)) # compute Euclidean distance between data points
cParam =1 # parameter of similarity function
S<-exp(-dZZ/cParam) #compute similarity matrix
S


AffMat<-function(S,k) #S-distance matrix and k-no of neighbours
{
 AM <- matrix(0,nrow=nrow(S),ncol=ncol(S))
 for(i in 1:nrow(S)){
 d <- sort(S[i,],decreasing=TRUE)
 for (t in 1:ncol(S))
 {
 if (S[i,t] < d[k])
 {
 AM[i,t]<-0
 AM[t,i]<-0
 }
 else
 {
 AM[i,t] <- S[i,t]
 AM[t,i] <- AM[i,t]
 }
 }
 }
 AM
}
A<-AffMat(S,11)
A

#Plot the Graph
#------------------------------------------------
library(shape)
library(diagram)
B<-A
diag(B) <- 0 # to avoid self loop in the graph.
B
pp <- plotmat(B, curve = 0, lwd = 1, box.lwd = 2, cex.txt = 0.8, box.type = "circle", box.prop = 0.1, arr.width=0, arr.pos = 0.5, shadow.size = 0, main = "Graph: connented components")
 
#Compute degree of Affinity matrix
#------------------------------------------------
D <- diag(apply(A, 1, sum)) # sum rows
D

#compute graph laplasian matrix (un-normalised)
#------------------------------------------------
L <- D - A
L

#find eigenvalues and eigenvectors
#------------------------------------------------
eigL<-eigen(L)
eigL
plot(eigL$values) 

#smallest eigenvalues of L
#------------------------------------------------
k<-5
Z<- eigL$vectors[,(ncol(eigL$vectors)-k+1):ncol(eigL$vectors)]
#plot data using the two eigenvectors
plot(Z)

# Perform k means clustering
#------------------------------------------------
library(stats)
km <- kmeans(Z, centers=k, nstart=20)
plot(Z, col=km$cluster)
plot(zz, col=km$cluster, main="Spectral Clustering with K=5")

#========================================================================================================================
# 7) Data Analysis of LMWT
#========================================================================================================================

#Extract data of LMWTdata
#------------------------------------------------
the.data <- as.matrix(read.csv("GBRData.csv", header =TRUE, sep = ","))
LMWTdata <- the.data[,4]

#Time series plot of LMWTdata
#------------------------------------------------
plot.ts(LMWTdata,main="Time series plot of LadyMusgrave-Water Temperature")

#Histogram of LMWTdata
#------------------------------------------------
hist(LMWTdata)

#MLE of LMWTdata
#------------------------------------------------
library(MASS)
library(mixtools)
fit1<-fitdistr(LMWTdata,"normal")
fit1

#Plot of density function distribution
#------------------------------------------------
x1 <- LMWTdata
y1 <- dnorm(LMWTdata,mean=24.98047914, sd=2.14893739)
plot(x1,y1, type="l", lwd=1, ylim=c(0, 0.2), col="Red", main="Plot of Density function distribution")

#Plot of density function distribution
#------------------------------------------------
dataSynth<-LMWTdata
#Histogram
hist(dataSynth)
#Gaussian mixture
mixmdl = normalmixEM(dataSynth) # default k=2 components
mixmdl
summary(mixmdl)
mixmdl$lambda # the mixing coefficients
mixmdl$mu
mixmdl$sigma

#Plot of density function distribution
#------------------------------------------------
plot(mixmdl,which=2)
lines(density(dataSynth), lty=2, lwd=2)

#Plot of the log likelihood values
#------------------------------------------------
plot(mixmdl$all.loglik, main="Plot of log likelihood values")