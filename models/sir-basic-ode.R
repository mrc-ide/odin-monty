deriv(S) <- -beta * S * I / N
deriv(I) <- beta * S * I / N - gamma * I
deriv(R) <- gamma * I

initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0

N <- parameter(1000)
I0 <- parameter(10)
beta <- parameter(0.2)
gamma <- parameter(0.1)
