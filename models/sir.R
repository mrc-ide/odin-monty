deriv(S) <- -beta * S * I / N
deriv(I) <- beta * S * I / N - sigma * I
deriv(R) <- sigma * I

initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0

N <- parameter(1e6)
I0 <- parameter(1)
beta <- parameter(4)
sigma <- parameter(2)
