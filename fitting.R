
set.seed(42)

## Load packages

library(odin2)
library(dust2)
library(monty)



## The data
# Download https://github.com/mrc-ide/odin-monty/blob/main/data/incidence.csv
# and ensure you read from the correct file path
data <- read.csv("data/incidence.csv")

plot(data, pch = 19, col = "red")



## Model with likelihood

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  cases <- data()
  cases ~ Poisson(incidence)
})



## Calculating likelihood
data <- dust_filter_data(data, time = "time")
filter <- dust_filter_create(sir, data = data, time_start = 0,
                             n_particles = 200, dt = 0.25)
dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))
dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))
dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))



## Filtered trajectories

dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2),
                    save_trajectories = TRUE)
y <- dust_likelihood_last_trajectories(filter)
y <- dust_unpack_state(filter, y)
matplot(data$time, t(y$incidence), type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence")
points(data, pch = 19, col = "red")



## Parameter packers

packer <- monty_packer(c("beta", "gamma"))
packer

packer$unpack(c(0.2, 0.1))

packer$pack(c(beta = 0.2, gamma = 0.1))

packer <- monty_packer(c("beta", "gamma"), fixed = list(I0 = 5))
packer$unpack(c(0.2, 0.1))

packer <- monty_packer(c("beta", "gamma"),
                       array = list(alpha = 3, delta = c(2, 2)),
                       fixed = list(eta = c(1, 3, 5)))
packer
packer$unpack(c(0.2, 0.1, 0.31, 0.32, 0.33, 0.15, 0.16, 0.17, 0.18))



## The monty DSL

prior <- monty_dsl({
  beta ~ Exponential(mean = 0.5)
  gamma ~ Exponential(mean = 0.3)
})
prior

prior$parameters

prior$density(c(0.2, 0.1))
dexp(0.2, 1 / 0.5, log = TRUE) + dexp(0.1, 1 / 0.3, log = TRUE)

rng <- monty_rng_create()
prior$direct_sample(rng)

prior$gradient(c(0.2, 0.1))


m <- monty_dsl({
  a ~ Normal(0, 1)
  b ~ Normal(a, 1)
})
m

m <- monty_dsl({
  b ~ Normal(a, 1)
  a ~ Normal(0, 1)
})

m <- monty_dsl({
  mu <- 10
  sd <- 2
  a ~ Normal(mu, sd)
})

m <- monty_dsl({
  a ~ Normal(0, 1)
  b ~ Normal(0, 1)
  mu <- (a + b) / 2
  c ~ Normal(mu, 1)
})

m <- monty_dsl({
  alpha ~ Normal(178, 20)
  beta ~ Normal(0, 10)
  sigma ~ Uniform(0, 50)
})

fixed <- list(alpha_mean = 170, alpha_sd = 20,
              beta_mean = 0, beta_sd = 10,
              sigma_max = 50)
m <- monty_dsl({
  alpha ~ Normal(alpha_mean, alpha_sd)
  beta ~ Normal(beta_mean, beta_sd)
  sigma ~ Uniform(0, sigma_max)
}, fixed = fixed)


fixed = list(gamma_mean = c(0.3, 0.25, 0.2))
m <- monty_dsl({
  alpha ~ Exponential(mean = 0.5)
  beta[, ] ~ Normal(0, 1)
  gamma[] ~ Exponential(mean = gamma_mean[i])
  dim(beta) <- c(2, 2)
  dim(gamma, gamma_mean) <- 3
}, fixed = fixed, gradient = FALSE)
m

fixed = list(gamma_mean = c(0.3, 0.25, 0.2),
             n_beta_1 = 2, n_beta_2 = 2, n_gamma = 3)
m <- monty_dsl({
  alpha ~ Exponential(mean = 0.5)
  beta[, ] ~ Normal(0, 1)
  gamma[] ~ Exponential(mean = gamma_mean[i])
  dim(beta) <- c(n_beta_1, n_beta_2)
  dim(gamma, gamma_mean) <- n_gamma
}, fixed = fixed, gradient = FALSE)
m

fixed = list(n_beta = 5)
m <- monty_dsl({
  beta[1] ~ Exponential(mean = 0.3)
  beta[2:5] ~ Exponential(mean = beta[i - 1])
  dim(beta) <- n_beta
}, fixed = fixed, gradient = FALSE)



## From a dust filter to a monty model

filter

packer <- monty_packer(c("beta", "gamma"))
likelihood <- dust_likelihood_monty(filter, packer)
likelihood



## Posterior from likelihood and prior

posterior <- likelihood + prior
posterior



## Create a sampler

vcv <- diag(2) * 0.2
vcv

sampler <- monty_sampler_random_walk(vcv)
sampler



## Let's sample!

samples <- monty_sample(posterior, sampler, 1000, n_chains = 3)
samples



## The result: diagnostics

samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)



## The results: parameters

bayesplot::mcmc_scatter(samples_df)



## The result: traceplots

bayesplot::mcmc_trace(samples_df)



## The result: density over time

matplot(drop(samples$density), type = "l", lty = 1)
matplot(drop(samples$density[-(1:100), ]), type = "l", lty = 1)



## Better mixing

vcv <- matrix(c(0.01, 0.005, 0.005, 0.005), 2, 2)
sampler <- monty_sampler_random_walk(vcv)
samples <- monty_sample(posterior, sampler, 2000, initial = samples,
                        n_chains = 4)
matplot(samples$density, type = "l", lty = 1)



## Better mixing: the results

samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)

bayesplot::mcmc_scatter(samples_df)

bayesplot::mcmc_trace(samples_df)



## Trajectories

likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE)
posterior <- likelihood + prior
samples <- monty_sample(posterior, sampler, 1000, n_chains = 4)

trajectories <- dust_unpack_state(filter,
                                  samples$observations$trajectories)
matplot(data$time, trajectories$incidence[, , 1], type = "l", lty = 1,
        col = "#00000044", xlab = "Time", ylab = "Infection incidence")
points(data, pch = 19, col = "red")

dim(samples$observations$trajectories)



## Saving a subset of trajectories

likelihood <- dust_likelihood_monty(filter, packer, 
                                    save_trajectories = c("I", "incidence"))
posterior <- likelihood + prior
samples2 <- monty_sample(posterior, sampler, 100, initial = samples)
dim(samples2$observations$trajectories)



## Thinning

samples <- monty_samples_thin(samples,
                              burnin = 500,
                              thinning_factor = 2)



## Fitting in deterministic mode

unfilter <- dust_unfilter_create(sir, data = data, time_start = 0, dt = 0.25)

dust_likelihood_run(unfilter, list(beta = 0.4, gamma = 0.2))
dust_likelihood_run(unfilter, list(beta = 0.4, gamma = 0.2))

likelihood <- dust_likelihood_monty(unfilter, packer, save_trajectories = TRUE)
posterior <- likelihood + prior
samples_det <- monty_sample(posterior, sampler, 1000, n_chains = 4)
samples_det <- monty_samples_thin(samples_det,
                                  burnin = 500,
                                  thinning_factor = 2)



## Stochastic v deterministic comparison

y <- dust2::dust_unpack_state(filter, samples$observations$trajectories)
incidence <- array(y$incidence, c(20, 1000))
matplot(data$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infection incidence", ylim = c(0, 75),
        main ="Stochastic fit")
points(data, pch = 19, col = "red")

y <- dust2::dust_unpack_state(filter, samples_det$observations$trajectories)
incidence <- array(y$incidence, c(20, 1000))
matplot(data$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infection incidence", ylim = c(0, 75),
        main ="Deterministic fit")
points(data, pch = 19, col = "red")

pars_stochastic <- array(samples$pars, c(2, 500))
pars_deterministic <- array(samples_det$pars, c(2, 500))
plot(pars_stochastic[1, ], pars_stochastic[2, ], ylab = "gamma", xlab = "beta",
     pch = 19, col = "blue")
points(pars_deterministic[1, ], pars_deterministic[2, ], pch = 19, col = "red")
legend("bottomright", c("stochastic fit", "deterministic fit"), pch = c(19, 19), 
       col = c("blue", "red"))



## Fitting a time-varying parameter
# Download https://github.com/mrc-ide/odin-monty/blob/main/data/incidence2.csv
# and ensure you read from the correct file path
data <- read.csv("data/incidence2.csv")
head(data)

plot(data, pch = 19, col = "red")

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  beta_t <- interpolate(beta_time, beta, "constant")
  beta <- parameter()
  beta_time <- parameter()
  dim(beta, beta_time) <- parameter(rank = 1)
  
  p_SI <- 1 - exp(-beta_t * I / N * dt)
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
  
  cases <- data()
  cases ~ Poisson(incidence)
})

packer <- monty_packer("gamma", array = list(beta = 2),
                       fixed = list(beta_time = c(0, 20)))

data <- dust_filter_data(data, time = "time")
filter <- dust_filter_create(sir, data = data, time_start = 0,
                             n_particles = 200, dt = 0.25)
likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE)

prior <- monty_dsl({
  beta[] ~ Exponential(mean = 0.3)
  gamma ~ Exponential(mean = 0.1)
  dim(beta) <- n_beta
}, fixed = list(n_beta = 2), gradient = FALSE)
posterior <- likelihood + prior
posterior

vcv <- diag(c(0.005, 0.01, 0.01))
sampler <- monty_sampler_random_walk(vcv)
samples <- monty_sample(posterior, sampler, 1000, initial = c(0.1, 0.5, 0.5), n_chains = 4)
samples <- monty_samples_thin(samples, thinning_factor = 2, burnin = 500)

y <- dust_unpack_state(filter, samples$observations$trajectories)
incidence <- array(y$incidence, c(50, 1000))
matplot(data$time, incidence, type = "l", lty = 1,
        col = "#00000044", xlab = "Time", ylab = "Infection incidence")
points(data, pch = 19, col = "red")



## Fitting to vector/array data
# Download https://github.com/mrc-ide/odin-monty/blob/main/data/incidence-age.csv
# and ensure you read from the correct file path
data_age <- read.csv("data/incidence-age.csv")
head(data_age)

sir_age <- odin({
  # Equations for transitions between compartments by age group
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  update(incidence[]) <- incidence[i] + n_SI[i]
  
  # Individual probabilities of transition:
  p_SI[] <- 1 - exp(-lambda[i] * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  
  # Calculate force of infection
  
  # age-structured contact matrix: m[i, j] is mean number
  # of contacts an individual in group i has with an
  # individual in group j per time unit
  m <- parameter()
  
  # here s_ij[i, j] gives the mean number of contacts an
  # individual in group i will have with the currently
  # infectious individuals of group j
  s_ij[, ] <- m[i, j] * I[j]
  
  # lambda[i] is the total force of infection on an
  # individual in group i
  lambda[] <- beta * sum(s_ij[i, ])
  
  # Draws from binomial distributions for numbers
  # changing between compartments:
  n_SI[] <- Binomial(S[i], p_SI[i])
  n_IR[] <- Binomial(I[i], p_IR)
  
  initial(S[]) <- S0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- 0
  initial(incidence[], zero_every = 1) <- 0
  
  S0 <- parameter()
  I0 <- parameter()
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  n_age <- parameter()
  dim(S, S0, n_SI, p_SI, incidence) <- n_age
  dim(I, I0, n_IR) <- n_age
  dim(R) <- n_age
  dim(m, s_ij) <- c(n_age, n_age)
  dim(lambda) <- n_age
  
  
  ## Likelihood
  cases <- data()
  dim(cases) <- n_age
  cases[] ~ Poisson(incidence[i])
})


## Formatting array data in the data frame

data_age$cases <- I(asplit(data_age[, c("cases_children", "cases_adult")], 1))
data_age$cases_children <- NULL
data_age$cases_adult <- NULL
head(data_age)


## Fitting the age-stratified model

packer <- monty_packer(c("beta", "gamma"),
                       fixed = list(S0 = c(990, 1000),
                                    I0 = c(10, 0),
                                    m = matrix(c(1.8, 0.4, 0.4, 1.2) / 2000, 2, 2),
                                    n_age = 2))


data_age <- dust_filter_data(data_age, time = "time")
filter <- dust_filter_create(sir_age, data = data_age, time_start = 0,
                             n_particles = 200, dt = 0.25)
likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE)

prior <- monty_dsl({
  beta ~ Exponential(mean = 0.3)
  gamma ~ Exponential(mean = 0.1)
})
posterior <- likelihood + prior

vcv <- matrix(c(0.01, 0.005, 0.005, 0.005), 2, 2)
sampler <- monty_sampler_random_walk(vcv)
samples <- monty_sample(posterior, sampler, 500, 
                        initial = c(0.3, 0.1), n_chains = 4)


y <- dust_unpack_state(filter, samples$observations$trajectories)
incidence <- array(y$incidence, c(2, 100, 2000))

matplot(data_age$time, incidence[1, , ], type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence (children)")
cases_children <- vapply(data_age$cases, "[[", numeric(1), "cases_children")
points(data_age$time, cases_children, pch = 19, col = "red")

matplot(data_age$time, incidence[2, , ], type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence (adults)")
cases_adult <- vapply(data_age$cases, "[[", numeric(1), "cases_adult")
points(data_age$time, cases_adult, pch = 19, col = "red")



## Projections and counterfactuals
# Download https://github.com/mrc-ide/odin-monty/blob/main/data/schools.csv
# and ensure you read from the correct file path
data <- read.csv("data/schools.csv")
plot(data, pch = 19, col = "red")

sis <- odin({
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
})

schools_time <- c(0, 50, 60, 120, 130, 170, 180)
schools_open <- c(1,  0,  1,   0,   1,   0,   1)


## Fitting to the SIS model

packer <- monty_packer(c("beta0", "gamma", "schools_modifier"),
                       fixed = list(schools_time = schools_time,
                                    schools_open = schools_open))

data <- dust_filter_data(data, time = "time")
filter <- dust_filter_create(sis, time_start = 0, dt = 0.25,
                             data = data, n_particles = 200)

prior <- monty_dsl({
  beta0 ~ Exponential(mean = 0.3)
  gamma ~ Exponential(mean = 0.1)
  schools_modifier ~ Uniform(0, 1)
})

vcv <- diag(c(2e-4, 2e-4, 4e-4))
sampler <- monty_sampler_random_walk(vcv)


likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE,
                                    save_state = TRUE, save_snapshots = 60)

posterior <- likelihood + prior

samples <- monty_sample(posterior, sampler, 500, initial = c(0.3, 0.1, 0.5),
                        n_chains = 4)
samples <- monty_samples_thin(samples, burnin = 100, thinning_factor = 8)



## Fit to data

y <- dust_unpack_state(filter, samples$observations$trajectories)
incidence <- array(y$incidence, c(150, 200))
matplot(data$time, incidence, type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence")
points(data, pch = 19, col = "red")



## Running projection using the end state

state <- array(samples$observations$state, c(3, 200))
pars <- array(samples$pars, c(3, 200))
pars <- lapply(seq_len(200), function(i) packer$unpack(pars[, i]))

sys <- dust_system_create(sis, pars, n_groups = length(pars), dt = 1)

dust_system_set_state(sys, state)
t <- seq(150, 200)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

matplot(data$time, incidence, type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence", xlim = c(0, 200))
matlines(t, t(y$incidence), col = "blue")
points(data, pch = 19, col = "red")



## Running counterfactual using the snapshot

snapshot <- array(samples$observations$snapshots, c(3, 200))
pars <- array(samples$pars, c(3, 200))
f <- function(i) {
  p <- packer$unpack(pars[, i])
  p$schools_time <- c(0, 50, 130, 170, 180)
  p$schools_open <- c(1, 0, 1, 0, 1)
  p
}
pars <- lapply(seq_len(200), f)
sys <- dust_system_create(sis, pars, n_groups = length(pars), dt = 1)

dust_system_set_state(sys, snapshot)
t <- seq(60, 150)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

matplot(data$time, incidence, type = "l", col = "#00000044", lty = 1,
        xlab = "Time", ylab = "Incidence")
matlines(t, t(y$incidence), col = "blue")
points(data, pch = 19, col = "red")
