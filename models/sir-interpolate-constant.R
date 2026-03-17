update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(incidence) <- incidence + n_SI

beta <- interpolate(beta_time, beta_value, "constant")
beta_time <- parameter()
beta_value <- parameter()
dim(beta_time, beta_value) <- parameter(rank = 1)

p_SI <- 1 - exp(-beta * I / N * dt)
p_IR <- 1 - exp(-gamma * dt)
n_SI <- Binomial(S, p_SI)
n_IR <- Binomial(I, p_IR)

initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0
initial(incidence, zero_every = 1) <- 0

N <- parameter(1000)
I0 <- parameter(10)
gamma <- parameter(0.1)
