library(tseries)
library(tfarima)
library(astsa)
library(TTR)
library(TSA)
library(MTS)
library(lmtest)
library(psych)
library(forecast)
library(MLmetrics)

#data input
data_ts <- read.csv(file.choose())
str(data_ts)
data_ts <- as.ts(data_ts[1:505,], frequency=12)

head(data_ts)

#train test split
training <- data_ts[-c(1:400),]
test <- data_ts[c(1:400),]

y <- as.ts(test[,2])
x <- as.ts(test[,3])

par(mfrow = c(2, 1))

plot(y)
plot(x)

MAPE_x <- training[,3]
MAPE_y <- training[,2]

#using augmented dickey fuller to check for variance stationarity, differencing if not stationary
adf.test(x)

x_diff <- diff(x,1)
adf.test(x_diff)


##using arima to check for initial model
par(mfrow = c(2,1))
acf(x_diff)
pacf(x_diff)

#the determined potential models is as follows. Eliminate using Coef significance, AIC, and MAPE
#(p,d,q)
#model x ((5-0),1,(2-0))
#model y ((3-0),1,(3-0))

model_x1 <- arima(x, order=c(5,1,2))
model_x2 <- arima(x, order=c(5,1,1))
model_x3 <- arima(x, order=c(5,1,0))

model_x4 <- arima(x, order=c(4,1,2))
model_x5 <- arima(x, order=c(4,1,1))
model_x6 <- arima(x, order=c(4,1,0))

model_x7 <- arima(x, order=c(3,1,2))
model_x8 <- arima(x, order=c(3,1,1))
model_x9 <- arima(x, order=c(3,1,0))

model_x10 <- arima(x, order=c(2,1,2))
model_x11 <- arima(x, order=c(2,1,1))
model_x12 <- arima(x, order=c(2,1,0))

model_x13 <- arima(x, order=c(1,1,2))
model_x14 <- arima(x, order=c(1,1,1))
model_x15 <- arima(x, order=c(1,1,0))

model_x16 <- arima(x, order=c(0,1,2))
model_x17 <- arima(x, order=c(0,1,1))
model_x18 <- arima(x, order=c(0,1,0))

coeftest(model_x1)
coeftest(model_x2) #*
coeftest(model_x3)
coeftest(model_x4)
coeftest(model_x5)
coeftest(model_x6)
coeftest(model_x7)
coeftest(model_x8)
coeftest(model_x9)
coeftest(model_x10)
coeftest(model_x11)
coeftest(model_x12)
coeftest(model_x13)
coeftest(model_x14) #*
coeftest(model_x15)
coeftest(model_x16)
coeftest(model_x17)
coeftest(model_x18)


AIC(model_x2) #*
AIC(model_x14) #*


MAPE(predict(model_x2, n.ahead = 105)$pred, MAPE_x)
MAPE(predict(model_x14, n.ahead = 105)$pred, MAPE_x) #

#model_x14 is chosen

#checking residual assumptions for model. both white noise and normality
residx14 <- residuals(model_x14)
n2=length(residx14)
mean2=mean(residx14)
sd2=sd(residx14)
res2=rnorm(n2,mean2,sd2)
ks.test(residx14,res2)

Box.test(residx14,lag=1,type=c("Ljung-Box"))
Box.test(residx14^2,lag=1,type=c("Ljung-Box")))

#cross correlation function to find b,r,s of transfer function
umx <- um(lapply(x, as.numeric), ar = 1, i = 1, ma = 1)
umy <- fit(umx,y)

a <- residuals(umx, y, method = "cond")
b <- residuals(umx, x, method = "cond")

pccf(a, b, um.x = umx, um.y = NULL, lag.max = 16)
CCF =ccf(a,b, ylab="CCF", main="Cross-correlations after prewhitening")

#Potential transfer function brs(0,1-2,1-5). After each potential model is tested using accuracy tests in another machine, it is found that, using accuracy tests such as MAPE, brs(0,1,5) is best

tfx <- tfest(y, x, delay = 0, p = 1, q = 5, um.x = umx, um.y = umy)
tfx

#model coefficients
print(tfx$theta)
print(tfx$phi)

tfmy <- tfarima::tfm(y, inputs=tfx, noise=um(ar=4))
printLagpol(tfmy$noise$phi)
printLagpol(tfmy$noise$theta)

#testing model accuracy against test set
p <- predict.tfm(tfmy, n.ahead = 105)
p
accuracy1 <- forecast::accuracy(p$z[401:505], MAPE_y)
accuracy1

#residual test of final model
Box.test(residuals(tfmy), type = "Ljung-Box")
Box.test(residuals(tfmy)^2, type = "Ljung-Box")
ks.test(residuals(tfmy), "pnorm", mean(residuals(tfmy)), sd(residuals(tfmy)))
