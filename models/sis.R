update(S) <- S - n_SI + n_IS
update(I) <- I + n_SI - n_IS
update(incidence) <- incidence + n_SI

initial(S) <- N - I0
initial(I) <- I0
initial(incidence, zero_every = 1) <- 0

schools <- interpolate(schools_time, schools_open, "constant")
schools_time <- parameter()
schools_open <- parameter()
dim(schools_time, schools_open) <- parameter(rank = 1)

beta <- ((1 - schools) * (1 - schools_modifier) + schools) * beta0

p_SI <- 1 - exp(-beta * I / N * dt)
p_IS <- 1 - exp(-gamma * dt)
n_SI <- Binomial(S, p_SI)
n_IS <- Binomial(I, p_IS)

N <- parameter(1000)
I0 <- parameter(10)
beta0 <- parameter(0.2)
gamma <- parameter(0.1)
schools_modifier <- parameter(0.6)

cases <- data()
cases ~ Poisson(incidence)
