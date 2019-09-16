library(forecast)

monthly.df = read.csv('data/monthly_in_situ_co2_mlo.csv', skip=58, sep=",", header=FALSE)
monthly.ts = ts(data=monthly.df$V5[2:737], frequency = 12, start=c(1958,3), end = c(2019,6))

for (i in 1:736) {
  if(monthly.ts[i]==-99.99) monthly.ts[i]= monthly.ts[i-1]+(monthly.ts[i-1]-monthly.ts[i-2])
}


monthly.dc = decompose(monthly.ts, type = "multiplicative")

plot(monthly.dc)

acf(monthly.dc$random, na.action = na.pass, lag.max = 36)

diff1 = diff(monthly.ts, lag=1)
diff1 = diff(diff1, lag=12)

plot(diff1)

acf(diff1, na.action = na.pass, lag.max = 36)
pacf(diff1, na.action = na.pass, lag.max = 36)

model = arima(diff1, order = c(0,0,1), seasonal = c(0,0,1), include.mean = FALSE)
summary(model)
model$aicc


auto.arima(diff1, d = 0, D = 0)
auto.arima(diff1, d = 0, D = 0, max.p=0, max.q=1, max.P=0, max.Q=1, ic = "aicc")
