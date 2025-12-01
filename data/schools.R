#!/usr/bin/env Rscript
root <- here::here()
dest <- file.path(root, "data", "schools.csv")

library(dust2)
library(odin2)

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
}, debug = TRUE, quiet = TRUE)

set.seed(1)
schools_time <- c(0, 50, 60, 120, 130, 170, 180)
schools_open <- c(1,  0,  1,   0,   1,   0,   1)
schools_modifier <- 0.6

pars <- list(schools_time = schools_time, schools_open = schools_open,
             schools_modifier = 0.6, beta0 = 0.2, gamma = 0.1)
sys <- dust_system_create(sis, pars, dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 150)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

incidence <- y$incidence
dat <- data.frame(time = t, cases = rpois(length(incidence), incidence))[-1, ]
write.csv(dat, dest, row.names = FALSE)
