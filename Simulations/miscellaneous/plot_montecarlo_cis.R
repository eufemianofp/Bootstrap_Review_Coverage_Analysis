
## Remove previous objects in the R session
rm(list = ls())

set.seed(3)  # symmetric noncoverage
set.seed(5)  # asymmetric noncoverage

nsim <- 100
mu <- 0
sd <- 1
n <- 100

tt <- matrix(NA, nrow = nsim, ncol = 2)
for (i in 1:nsim) {
  x <- rnorm(n, mean = mu, sd = sd)
  tt[i, ] <- t.test(x, conf.level = 0.90)$conf.int
}

yy <- matrix(data = rep(seq(1/(nsim+1), 1-1/(nsim+1), 1/(nsim+1)), each = 2), 
             nrow = nsim, ncol = 2, byrow = TRUE)

contains <- tt[, 1] < mu & tt[, 2] > mu
contains[contains == TRUE] <- "blue"
contains[contains == FALSE] <- "red"


# pdf(file = "symmetric_cis2.pdf", width = 15/2.54, height = 15/2.54)
# pdf(file = "symmetric_cis2.pdf")

plot.new()
plot.window(xlim = c(-3*sd/sqrt(n), 3/sqrt(n)), ylim = c(0, 1))
par(mar=c(2, 1, 1, 1), oma=c(0, 0, 0, 0), cex = 2.5)
axis(1, at=0, labels=expression(theta)) # set x axis values
abline(v = mu, col = "black", lty = "dashed")
segments(x0 = tt[, 1], y0 = yy[, 1], x1 = tt[, 2], y1 = yy[, 2], col = contains, lwd = 1.2)

dev.off()
