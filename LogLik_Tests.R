lincomb <- function(x, theta) {
  x %*% t(theta)
}

logf <- function(x, theta) {
  1 / (1 + exp(-lincomb(x, theta)))
}

likf <- function(x, y, t) {
  temp <- lincomb(x, t);
  loglik <- 0
  for (i in seq(1, nrow(temp))) {
    print(paste("Temp[i] = ", temp[i], " Y[i] = ", y[i]))
    z <- logf(x, t)[i]
    print((y[i]*log(z) + ((1 - y[i])*log(1 - z))))
    loglik <- loglik + ((y[i]*log(z) + ((1 - y[i])*log(1 - z))))
  }
  return(loglik) 
}

t1 <- matrix(c(1, 1), nrow = 1, ncol = 2, byrow = TRUE)
x1 <- matrix(c(1,1, 0,0), nrow = 2, ncol = 2, byrow = TRUE)
y1 <- matrix(c(1, 0), nrow = 2, ncol = 1, byrow = TRUE)

lincomb(x1, t1)
logf(x1, t1)
likf(x1, y1, t1)

t2 <- matrix(c(0.25, 7.5), nrow = 1, ncol = 2, byrow = TRUE)
x2 <- matrix(c(1.1,0.2, -6.9,4.20), nrow = 2, ncol = 2, byrow = TRUE)
y2 <- matrix(c(0, 1), nrow = 2, ncol = 1, byrow = TRUE)

lincomb(x2, t2)
logf(x2, t2)
likf(x2, y2, t2)

t3 <- matrix(c(1, -1), nrow = 1, ncol = 2, byrow = TRUE)
x3 <- matrix(c(1,23, 1,22, 1,21), nrow = 3, ncol = 2, byrow = TRUE)
y3 <- matrix(c(1, 0, 1), nrow = 3, ncol = 1, byrow = TRUE)

lincomb(x3, t3)
logf(x3, t3)
likf(x3, y3, t3)

likf_d <- function(x, y, t) {
  temp <- lincomb(x, t);
  ld <- c(0, 0)
  for (j in seq(1, length(ld))) {
    for (i in seq(1, nrow(x))) {
      q <- logf(x, t)[i]
      print((y[i] - q) * x[i, j])
      ld[j] <- ld[j] + ((y[i] - q) * x[i, j])
    }
  }
  return(ld)
}

likf_d(x1, y1, t1)
likf_d(x2, y2, t2)
likf_d(x3, y3, t3)
