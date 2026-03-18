# Equations for transitions between compartments by age group
update(S[]) <- S[i] - n_SI[i]
update(I[]) <- I[i] + n_SI[i] - n_IR[i]
update(R[]) <- R[i] + n_IR[i]
update(incidence) <- incidence + sum(n_SI)

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
initial(incidence, zero_every = 1) <- 0

S0 <- parameter()
I0 <- parameter()
beta <- parameter(0.2)
gamma <- parameter(0.1)

n_age <- parameter()
dim(S, S0, n_SI, p_SI) <- n_age
dim(I, I0, n_IR) <- n_age
dim(R) <- n_age
dim(m, s_ij) <- c(n_age, n_age)
dim(lambda) <- n_age
