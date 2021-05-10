

rm(list = ls())

## Load libraries
library(MASS)


## Empirical distribution of correlation coefficient r (for true correlation = theta and sample size = n)

# Define model
n <- 10
mu <- c(0, 0)
cor <- 0
Sigma <- matrix(data = c(1, cor, cor, 1), nrow = 2)

# Monte Carlo simulation from model
R <- 1e4
rs <- numeric(length = R)
for (i in 1:R) {
  x <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  rs[i] <- cor(x[, 1], x[, 2])
  # rs[i] <- atanh(cor(x[, 1], x[, 2]))
}
hist(rs, 
     breaks = 200, 
     freq = FALSE,
     main = paste0("Cor = ", cor, ", n = ", n), xlim = c(-1, 1))
# qqnorm(rs)


# Theoretical distribution of sample correlation coefficient r (for true correlation = cor and sample size = n)
r <- seq(from = -1, to = 1, by = 0.001)
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}

plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0.5, n = 10"),
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))


## Generate pdf with theoretical distribution of sample correlation coefficient
pdf(file = paste0("cor_pdf.pdf"), 
    width = 21/2.54, height = 16/2.54, pointsize = 16)

par(mfrow = c(3, 2))
par(mar=c(4, 4, 3, 0.5))


# n = 10, cor = 0
n <- 10
cor <- 0
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}
plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0,  n = 10"),# cex.main = 0.9, #font.main = 2, 
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))

# n = 100, cor = 0
n <- 100
cor <- 0
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}
plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0,  n = 100"),
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))

# n = 10, cor = 0.5
n <- 10
cor <- 0.5
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}
plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0.5,  n = 10"),
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))

# n = 100, cor = 0.5
n <- 100
cor <- 0.5
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}
plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0.5,  n = 100"),
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))

# n = 10, cor = 0.9
n <- 10
cor <- 0.9
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}
plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0.9,  n = 10"),
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))

# n = 100, cor = 0.9
n <- 100
cor <- 0.9
fr <- numeric(length = length(r))
for (i in 1:length(r)) {
  f <- function(w) {
    return( 1 / (cosh(w) - cor*r[i])^(n-1) )
  }
  fr[i] <- (n-2)*(1-cor^2)^((n-1)/2) *(1-r[i]^2)^((n-4)/2)/pi * integrate(f, lower = 0, upper = Inf)$value
}
plot(r, fr, type = "l", col = "blue", 
     main = expression(rho * " = 0.9, n =  100"),
     xlab = "r", ylab = "Density",
     xlim = c(-1, 1))


dev.off()








