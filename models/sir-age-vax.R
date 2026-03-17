# Equations for transitions between compartments by age group
update(S[, ]) <- new_S[i, j]
update(I[, ]) <- I[i, j] + n_SI[i, j] - n_IR[i, j]
update(R[, ]) <- R[i, j] + n_IR[i, j]
update(incidence) <- incidence + sum(n_SI)

# Individual probabilities of transition:
p_SI[, ] <- 1 - exp(-rel_susceptibility[j] * lambda[i] * dt) # S to I
p_IR <- 1 - exp(-gamma * dt) # I to R
p_vax[, ] <- 1 - exp(-eta[i, j] * dt)

# Force of infection
m <- parameter() # age-structured contact matrix
s_ij[, ] <- m[i, j] * sum(I[j, ])
lambda[] <- beta * sum(s_ij[i, ])

# Draws from binomial distributions for numbers changing between
# compartments:
n_SI[, ] <- Binomial(S[i, j], p_SI[i, j])
n_IR[, ] <- Binomial(I[i, j], p_IR)

# Nested binomial draw for vaccination in S
# Assume you cannot move vaccine class and get infected in same step
n_S_vax[, ] <- Binomial(S[i, j] - n_SI[i, j], p_vax[i, j])
new_S[, 1] <- S[i, j] - n_SI[i, j] - n_S_vax[i, j] + n_S_vax[i, n_vax]
new_S[, 2:n_vax] <- S[i, j] - n_SI[i, j] - n_S_vax[i, j] + n_S_vax[i, j - 1]

initial(S[, ]) <- S0[i, j]
initial(I[, ]) <- I0[i, j]
initial(R[, ]) <- 0
initial(incidence, zero_every = 1) <- 0

# User defined parameters - default in parentheses:
S0 <- parameter()
I0 <- parameter()
beta <- parameter(0.0165)
gamma <- parameter(0.1)
eta <- parameter()
rel_susceptibility <- parameter()

# Dimensions of arrays
n_age <- parameter()
n_vax <- parameter()
dim(S, S0, n_SI, p_SI) <- c(n_age, n_vax)
dim(I, I0, n_IR) <- c(n_age, n_vax)
dim(R) <- c(n_age, n_vax)
dim(m, s_ij) <- c(n_age, n_age)
dim(lambda) <- n_age
dim(eta) <- c(n_age, n_vax)
dim(rel_susceptibility) <- c(n_vax)
dim(p_vax, n_S_vax, new_S) <- c(n_age, n_vax)
