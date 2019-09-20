library(forecast)
library(EnvStats)
library(tseries)

plot_path=paste0(getwd(),"/plot/")


# Load data
monthly.df = read.csv('data/monthly_in_situ_co2_mlo.csv', skip=58, sep=",", header=FALSE)
monthly.ts = ts(data=monthly.df$V5[2:737], frequency = 12, start=c(1958,3), end = c(2019,6))

# Interpolate missing data, simple linear interpolation
for (i in 1:736) {
  if(monthly.ts[i]==-99.99) monthly.ts[i]= monthly.ts[i-1]+(monthly.ts[i-1]-monthly.ts[i-2])
}

# Decompose using moving average for trend and average months for seasonal component
filter = c(1/2,rep(1,10),1/2)/12
print(filter)
monthly.dc = decompose(monthly.ts, type = "multiplicative")#, filter = filter)



# Plot series with decomposition, acf and pacf of random part
pdf(paste0(plot_path,"decomposition.pdf"), width=8, height=8)
plot(monthly.dc)
dev.off()

pdf(paste0(plot_path,"acf_decomposed.pdf"), width=8, height=5)
acf(monthly.dc$random, na.action = na.pass, lag.max = 36)
dev.off()



# Transform and difference series, plot result and print mean, plot acf/pacf
model1 = monthly.ts
model1 = BoxCox(model1, lambda = 0)
model1 = diff(model1, lag=1, differences = 1)

pdf(paste0(plot_path,"simple_differenced.pdf"), width=8, height=5)
plot(model1)
dev.off()

model1 = diff(model1, lag=12, differences = 1)
print(mean(model1))

adf.test(model1)
acf(model1, na.action = na.pass, lag.max = 36)
pacf(model1, na.action = na.pass, lag.max = 36)

# Fit ARIMA-model, print summary, plot residuals with qq and acf/pacf
model = arima(model1, order = c(0,0,1), seasonal = c(0,0,1), include.mean = FALSE)
summary(model)
plot(model$residuals)
qqPlot(c(model$residuals))
acf(model$residuals)
pacf(model$residuals)

# Fit auto-arima to find "best" model, with and wothout restrains
auto.arima(model1, d = 0, D = 0, ic="aicc")


