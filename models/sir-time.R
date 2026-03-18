update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(incidence) <- incidence + n_SI

seed <- if (time == seed_time) seed_size else 0

p_SI <- 1 - exp(-beta * I / N * dt)
p_IR <- 1 - exp(-gamma * dt)
n_SI <- min(seed + Binomial(S, p_SI), S)
n_IR <- Binomial(I, p_IR)

initial(S) <- N
initial(I) <- 0
initial(R) <- 0
initial(incidence, zero_every = 1) <- 0

N <- parameter(1000)
beta <- parameter(0.2)
gamma <- parameter(0.1)
seed_time <- parameter()
seed_size <- parameter()
