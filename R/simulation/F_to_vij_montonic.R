f <- function(x, eps = 0.5, c = 2) {
  pnorm(-eps - c * x) + pnorm(-eps + c * x)
}
x <- seq(-3, 3, length = 100)
x2 <- x^2
plot(x2, f(x, eps = 0.1), type = "l")
