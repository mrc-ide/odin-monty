#!/usr/bin/env Rscript
root <- here::here()
dest <- file.path(root, "data", "incidence.csv")

library(dust2)
library(odin2)

sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI

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
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
}, debug = TRUE, quiet = TRUE)

set.seed(1)
pars <- list(beta = 1, gamma = 0.6)
sys <- dust_system_create(sir(), pars, dt = 0.25)
dust_system_set_state_initial(sys)
t <- seq(0, 20)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)

incidence <- y$incidence
dat <- data.frame(time = t, cases = rpois(length(incidence), incidence))[-1, ]
write.csv(dat, dest, row.names = FALSE)
