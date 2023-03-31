
# Good workflow for project:
# Apply methods of Ch. 2 (e.g. regression) to remove the overall trend
# Apply methods of Ch. 4 to remove the seasonality
# Then fit an ARIMA model to what is left


# Load libraries
library(data.table)
library(MASS)
library(astsa)


################################################################################
# SET-UP/LOAD DATA
################################################################################

# Import data
df <- read.csv("cookies.csv")


# Overview of data (maybe sub-index df so we just have the data through 2019 for simplicity)
#df <- df[1:480, ]
tail(df)


# Get the column we care about ("Value")
cookies <- df[5]
nrow(cookies)
summary(cookies)


# Create time series object
xt_total <- ts(cookies, start = c(1980,1), deltat = 1/12)

# Modify date in order to use as.Date() function
df[,4] <- paste(df[,4], "01", sep = " ")    # Just designate everything as the 1st of the month
t_total <- as.Date(df[,4], format = "%Y %b %d")
n_total <- length(t_total)

# Train-test split the series
xt <- ts(cookies[1:468, ], start = c(1980,1), deltat = 1/12)
t <- t_total[t_total < as.Date("01/01/2019", "%m/%d/20%y")]
n <- length(t)


# Create time plot (of training data)
plot(t, xt, type = "l", xlab = "year", ylab = "cookie price per lb. ($)")

# Show/plot composition of the series
plot(decompose(xt))

# Estimate the preliminary sample ACF and PACF
# The sample autocovariance function is defined as 
# $$\gamma^{\hat}(h) = n^{-1}\sum_{t=1}^{n-h}{(x_{t+h}-\bar{x})(x_t-\bar{x})}$$
# with $\gamma^{\hat}(-h) = \gamma^{\hat}(h)$ for $h = 0, 1, \hdots, n-1$.
round(acf(xt, lag.max = NULL, plot = TRUE)$acf[-1], 3)
round(pacf(xt, lag.max = NULL, plot = TRUE)$acf[-1], 3)


################################################################################
# REMOVE (LINEAR) TREND
################################################################################

# Use regression to remove the mean (trends)
#m <- lm(xt ~ t)
#summary(m)
# Use this instead:
m <- lsfit(x = julian(t, origin = t[1]), y = xt)
plot(t, xt, xlab = "year", ylab = "price of cookies", type = "l")
lines(t, m$coef[1] + m$coef[2] * julian(t, origin = t[1]), col = "blue", lwd = 2) 


# Make sure we only need to remove LINEAR relationship (look for signal in noise?)
detrended_xt <- m$residuals 
#plot(t, detrended_xt, xlab = "year", ylab = "residuals", type = "l") 


################################################################################
# REMOVE SEASONALITY
################################################################################

# Get rid of seasonal trend/cycles
spectrum(detrended_xt, log = "no", spans = c(12, 12))    # nonparametric
spectrum(detrended_xt, method = "ar", add = T, col = "blue", lty = 2)    # parametric
legend("topright", cex = 0.7, c("nonpar", "par"), col = c("black", "blue"), lty = 1:2) 


# Nothing clear from spectrum, however we may suspect an annual cycle of length 12 using domain knowledge/intuition
# Some years resemble each other more closely than others, still some relationship there...
months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
plot(t[421:432], detrended_xt[421:432], ylim = c(-0.25, 0.25), type = "c", xlab = "month", ylab = "residual")    # 2015
points(t[421:432], detrended_xt[421:432], pch = months, col = 1:4)
plot(t[457:468], detrended_xt[457:468], ylim = c(-0.25, 0.25), type = "c",xlab = "month", ylab = "residual")    # 2018
points(t[457:468], detrended_xt[457:468], pch = months, col = 1:4)

detrended_matrix <- matrix(detrended_xt, byrow = T, ncol = 12) 
monthly_means <- colMeans(detrended_matrix) 
plot(1:12, monthly_means, xlab = "month", ylab = "mean price residual") 


## Now remove this annual cycle and show the residuals. 
final_xt <- detrended_xt 
for(i in 1:n){
  j <- i %% 12 
  if(j == 0) j <- 12
  final_xt[i] <- detrended_xt[i] - monthly_means[j]
} 
plot(t, final_xt, xlab = "year", ylab = "residuals", type = "l") 
#plot(t, detrended_xt - final_xt, type = "l", xlab = "year", ylab = "monthly seasonal component") 


################################################################################
# FIT ARMA MODELS
################################################################################

# Grid search for ARMA parameters: through all 10x10 settings of ARMA(p, q) and calculate the AICs
aics = matrix(100, nrow = 10, ncol = 10) 
for(i in 1:10){
  for(j in 1:10){
    aics[i,j] = sarima(final_xt, p = i-1, q = j-1, d = 0)$AIC 
  }
} 

# Some models may not converge --> also try to fill in the AICs matrix randomly (using class notes)
for(i in 1:1000){ 
  i <- sample(10)[1]; j <- sample(10)[1] 
  if(aics[i,j] > 99){
    aics[i,j] <- sarima(final_xt, p = i-1, q = j-1, d = 0)$AIC 
  }
}


# Investigate the (top 5) best fit model(s)
which(aics == min(aics))
aics[8, 4] 
which(aics <= sort(aics)[5])
#plot(c(0,18),c(3.96,4.01), type="n", xlab="p+q", ylab="AIC") 
#for(i in 1:10){
#for(j in 1:10){
#points(i+j-2, aics[i,j]) 
#}
#}

# Our optimal ARMA model
a <- sarima(final_xt, p = 7, q = 3, d = 0)
# AIC of our selected model
a$AIC
b <- a$fit$coef
# Print ARMA coefficients
b


# Analyze the fit of our model and calculate RMSE
residuals <- rep(0, n) 
preds <- final_xt
for(i in 8:n){
  preds[i] <- sum(b[1:5] * final_xt[(i-1):(i-7)]) + sum(b[8:10] * residuals[(i-1):(i-3)]) + b[11] 
  residuals[i] <- final_xt[i] - preds[i] 
} 
plot(t, preds, type = "l", col = "green", ylab = "value", xlab = "year") 
lines(t, final_xt)
legend("topleft", cex = 0.6, c("Actual", "Predicted"), lty = 1, col = c("black", "green"))
# Find RMSE on training set (dollar units)
sqrt(mean(residuals^2)) 


################################################################################
# MAKE PREDICTIONS (FORECASTING)
################################################################################

## Create final, detrended, deseasoned full_xt: all the data inc. the testing portion, minus linear trend and cycles
full_xt <- xt_total - (m$coefficients[1] + m$coefficients[2] * julian(t_total, origin = t_total[1]) + monthly_means)


# Forecast using our model
preds <- full_xt 
residuals <- rep(0, n_total) 
for(i in 8:n_total){
  preds[i] <- sum(b[1:7] * full_xt[(i-1):(i-7)]) + sum(b[8:10] * residuals[(i-1):(i-3)]) + b[11] 
  residuals[i] <- full_xt[i] - preds[i] 
}

preds_total <- (m$coefficients[1] + m$coefficients[2] * julian(t_total, origin = t_total[1]) + monthly_means + preds) 
plot(t_total, preds_total, col = "green", type = "l", xlab = "year", ylab = "cookie price") 
lines(t_total, xt_total) 
legend("topleft", cex = 0.55, c("Actual", "Predicted"), lty = 1, col = c("black", "green"))


## Look closely at the testing data only
plot(t_total[(n+1):n_total], xt_total[(n+1):n_total], type = "l", xlab = "year", ylab = "cookie price") 
lines(t_total[(n+1):n_total], preds_total[(n+1):n_total], col = "green") 
legend("topleft", cex = 0.55, c("Actual", "Predicted"), lty = 1, col = c("black", "green"))
## Test RMSE of 1-day horizon forecasts
sqrt(mean(residuals[(n+1):n_total]^2))


xt_total[(n+1):n_total]
preds_total[(n+1):n_total]


