
set.seed(42)

## Load packages

library(odin2)
library(dust2)
library(monty)



## The data

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

packer <- monty_packer(array = c(beta = 3, gamma = 3))
packer
packer$unpack(c(0.2, 0.21, 0.22, 0.1, 0.11, 0.12))



## Priors

prior <- monty_dsl({
  beta ~ Exponential(mean = 0.5)
  gamma ~ Exponential(mean = 0.3)
})
prior

monty_model_density(prior, c(0.2, 0.1))

dexp(0.2, 1 / 0.5, log = TRUE) + dexp(0.1, 1 / 0.3, log = TRUE)



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
        xlab = "Time", ylab = "Infection incidence", ylim = c(0, 75))
points(data, pch = 19, col = "red")

y <- dust2::dust_unpack_state(filter, samples_det$observations$trajectories)
incidence <- array(y$incidence, c(20, 1000))
matplot(data$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infection incidence", ylim = c(0, 75))
points(data, pch = 19, col = "red")

pars_stochastic <- array(samples$pars, c(2, 500))
pars_deterministic <- array(samples_det$pars, c(2, 500))
plot(pars_stochastic[1, ], pars_stochastic[2, ], ylab = "gamma", xlab = "beta",
     pch = 19, col = "blue")
points(pars_deterministic[1, ], pars_deterministic[2, ], pch = 19, col = "red")
legend("bottomright", c("stochastic fit", "deterministic fit"), pch = c(19, 19), 
       col = c("blue", "red"))



## Projections and counterfactuals

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

filter <- dust_filter_create(sis, time_start = 0, dt = 1,
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
