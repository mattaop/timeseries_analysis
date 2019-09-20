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
model = arima(diff1, order = c(0,0,1), seasonal = c(0,0,1), include.mean = FALSE)
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
observed_residuals <- model$residuals

simulate_arima_values <- function(parameters, residuals){
  sampled_residuals <- sample(residuals, size=736+13, replace=TRUE)
  z <- vector(length=736)
  for(i in 14:(736+13)){
    # index <- sample(14:736, size=1) # Sample T random residuals from 
    z[i-13] <- (sampled_residuals[i]
               +parameters[1]*sampled_residuals[i-1]
               +parameters[2]*sampled_residuals[i-12]
               +parameters[1]*sampled_residuals[2]*residuals[i-13])
  }
  return(z)
}
simulate_sequence <- function(data, model, observed_residuals){
  simulated_data <- simulate(model, nsim=736)
  simulated_model <- arima(simulated_data, order = c(0,0,1), seasonal = c(0,0,1), include.mean = FALSE) # Fit betas for LS to the sequence
  # simulated_data  <- simulate_arima_values(model$coef, observed_residuals)
  # simulated_model <- arima(ts(data=simulated_data, frequency = 12), order = c(0,0,1), seasonal = c(0,0,1), include.mean = FALSE) # Fit betas for LS to the sequence
  beta <- simulated_model$coef
  return (beta)
}
beta <- simulate_sequence(diff1, model, observed_residuals)
print(beta)

sample_parameters <- function(B, diff1, model, observed_residuals){
  sampled_beta <- matrix(nrow=2, ncol=B) # Matrix to store sampled beta
  for (i in 1:B){
    # Return betas for random generated sequences
    beta_hat <- simulate_sequence(diff1, model, observed_residuals) 
    sampled_beta[,i] <- beta_hat # Store betas
  }
  observed_mean <- rowSums(sampled_beta)/B # Compute observed mean for sampled betas for LS
  observed_bias <- model$coef - observed_mean # Compute observed bias for sampled betas for LA
  observed_variance <- c(var(sampled_beta[1,]), var(sampled_beta[2,])) # Compute 
  # observed variance for sampled betas
  return(list(beta = sampled_beta, mean = observed_mean, bias = observed_bias, variance = observed_variance))
}
B = 1000
bootstrap_samples = sample_parameters(B, diff1, model, observed_residuals)
# print(bootstrap_samples$beta)
print(bootstrap_samples$mean)
print(bootstrap_samples$bias)
print(bootstrap_samples$variance)
hist(bootstrap_samples$beta[1,], freq = F, breaks=40, main="Histogram of bootstrap samples", 
     xlab="Sampled beta_1 values", ylab="Denisty", sub="Figure 1: Histogram of bootstrap samples of ma1 parameter.")
abline(v=c(model$coef[1], bootstrap_samples$mean[1], quantile(bootstrap_samples$beta[1,], c(0.025, 0.975))), col=c("red", "yellow", "blue", "blue"))
hist(bootstrap_samples$beta[2,], freq = F, breaks=40, main="Histogram of bootstrap samples", 
     xlab="Sampled beta_2 values", ylab="Denisty", sub = "Figure 2: Histogram of bootstrap samples of sma1 parameter.")
abline(v=c(model$coef[2], bootstrap_samples$mean[2], quantile(bootstrap_samples$beta[2,], c(0.025, 0.975))), col=c("red", "yellow", "blue", "blue"))
