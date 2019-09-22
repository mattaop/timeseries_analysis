library(forecast)
library(EnvStats)
library(tseries)

plot_path=paste0(getwd(),"/plot/")


# Load data
monthly.df = read.csv('data/monthly_in_situ_co2_mlo.csv', skip=58, sep=",", header=FALSE)
print(monthly.df)
monthly.ts = ts(data=monthly.df$V5[2:737], frequency = 12, start=c(1958,3), end = c(2019,6))

# Interpolate missing data, simple linear interpolation
for (i in 1:736) {
  if(monthly.ts[i]==-99.99) monthly.ts[i]= monthly.ts[i-1]+(monthly.ts[i-1]-monthly.ts[i-2])
}

# Decompose using moving average for trend and average months for seasonal component
filter = c(1/2,rep(1,10),1/2)/12
print(filter)
monthly.dc = decompose(monthly.ts, type = "multiplicative", filter = filter)



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
model = arima(model1, order = c(2,0,1), seasonal = c(0,0,1), include.mean = FALSE)
summary(model)
plot(model$residuals)
qqPlot(c(model$residuals))
acf(model$residuals)
pacf(model$residuals)
# futurVal <- forecast.Arima(model, h=10, level=c(99.5))
# plot.forecast(futurVal)
# Fit auto-arima to find "best" model, with and wothout restrains

# auto.arima(diff1, d = 0, D = 0, ic="aicc")
print(model$coef)
print(model$var.coef)

bootstrap <- function(data, model){
  bootstrap_residuals <- sample(model$residuals, size=736+13, replace=TRUE)
  bootstrap_index <- sample(1:(length(data)-1), 1) # Get index for a random sample from the data, but not the last
  x <- vector(length=736+2)
  x[1] <- data[bootstrap_index] # Get a random x from the data
  x[2] <- data[bootstrap_index+1]
  for(i in 1:736){
    j <- i + 2 # Shift timeseries to be indexed from -1 to 736
    k <- i + 13 # Shift residual array to be indexed from -12 to 736
    x[j] <- (model$coef['ar1']*x[j-1]
             +model$coef['ar2']*x[j-2]
             +bootstrap_residuals[k]
             +model$coef['ma1']*bootstrap_residuals[k-1]
             +model$coef['sma1']*bootstrap_residuals[k-12]
             +model$coef['ma1']*model$coef['sma1']*bootstrap_residuals[k-13]
             )
  }
  generated_timeseries <- ts(data=tail(x, -2), frequency = 12, start=c(1958,3), end = c(2019,6))
  return(generated_timeseries)
}
simulate_sequence <- function(data, model){
  #simulated_data <- simulate(model, nsim=736)
  simulated_data  <- bootstrap(data, model)
  simulated_model <- arima(simulated_data, order = c(2,0,1), seasonal = c(0,0,1), include.mean = FALSE)
  beta <- simulated_model$coef
  return (beta)
}
beta <- simulate_sequence(model1, model)

sample_parameters <- function(B, diff1, model){
  sampled_beta <- matrix(nrow=4, ncol=B) # Matrix to store sampled beta
  for (i in 1:B){
    # Return betas for random generated sequences
    beta_hat <- simulate_sequence(diff1, model) 
    sampled_beta[,i] <- beta_hat # Store betas
  }
  observed_mean <- rowSums(sampled_beta)/B # Compute observed mean for sampled betas for LS
  observed_bias <- model$coef - observed_mean # Compute observed bias for sampled betas for LA
  observed_variance <- c(var(sampled_beta[1,]), var(sampled_beta[2,]),
                         var(sampled_beta[3,])
                         ,var(sampled_beta[4,])
                         ) # Compute 
  # observed variance for sampled betas
  return(list(beta = sampled_beta, mean = observed_mean, bias = observed_bias, variance = observed_variance))
}
B = 1000
bootstrap_samples = sample_parameters(B, model1, model)
# print(bootstrap_samples$beta)
print(bootstrap_samples$mean)
print(bootstrap_samples$bias)
print(bootstrap_samples$variance)
hist(bootstrap_samples$beta[1,], freq = F, breaks=40, main="Histogram of simulated samples", 
     xlab="Ar1", ylab="Denisty")
abline(v=c(model$coef[1], bootstrap_samples$mean[1], quantile(bootstrap_samples$beta[1,], c(0.025, 0.975))), col=c("red", "yellow", "blue", "blue"))
hist(bootstrap_samples$beta[2,], freq = F, breaks=40, main="Histogram of simulated samples", 
     xlab="Ar2", ylab="Denisty")
abline(v=c(model$coef[2], bootstrap_samples$mean[2], quantile(bootstrap_samples$beta[2,], c(0.025, 0.975))), col=c("red", "yellow", "blue", "blue"))
hist(bootstrap_samples$beta[3,], freq = F, breaks=40, main="Histogram of simulated samples", 
     xlab="Ma1", ylab="Denisty")
abline(v=c(model$coef[3], bootstrap_samples$mean[3], quantile(bootstrap_samples$beta[3,], c(0.025, 0.975))), col=c("red", "yellow", "blue", "blue"))
hist(bootstrap_samples$beta[4,], freq = F, breaks=40, main="Histogram of simulated samples", 
     xlab="Sma1", ylab="Denisty")
abline(v=c(model$coef[4], bootstrap_samples$mean[4], quantile(bootstrap_samples$beta[4,], c(0.025, 0.975))), col=c("red", "yellow", "blue", "blue"))
