library(forecast)
library(EnvStats)

# Load data
monthly.df = read.csv('data/monthly_in_situ_co2_mlo.csv', skip=58, sep=",", header=FALSE)
monthly.ts = ts(data=monthly.df$V5[2:737], frequency = 12, start=c(1958,3), end = c(2019,6))

# Interpolate missing data, simple linear interpolation
for (i in 1:736) {
  if(monthly.ts[i]==-99.99) monthly.ts[i]= monthly.ts[i-1]+(monthly.ts[i-1]-monthly.ts[i-2])
}

# Decompose using moving average for trend and average months for seasonal component
monthly.dc = decompose(monthly.ts, type = "additive")

# Plot series with decomposition, acf and pacf of random part
plot(monthly.dc)
acf(monthly.dc$random, na.action = na.pass, lag.max = 36)
pacf(monthly.dc$random, na.action = na.pass, lag.max = 36)


# Transform and difference series, plot result and print mean, plot acf/pacf
model1 = BoxCox(monthly.ts, lambda = 0)
diff1 = diff(model1, lag=1, differences = 1)
plot(diff1)
diff1 = diff(diff1, lag=12, differences = 1)
plot(diff1)
print(mean(diff1))
acf(diff1, na.action = na.pass, lag.max = 36)
pacf(diff1, na.action = na.pass, lag.max = 36)

# Fit ARIMA-model, print summary, plot residuals with qq and acf/pacf
model = arima(diff1, order = c(2,0,2), seasonal = c(2,0,0), include.mean = TRUE)
summary(model)
plot(model$residuals)
qqPlot(c(model$residuals))
acf(model$residuals)
pacf(model$residuals)

# Fit auto-arima to find "best" model, with and wothout restrains
auto.arima(diff1, d = 0, D = 0, ic="aicc")


