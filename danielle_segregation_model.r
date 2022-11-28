# Author(s): Danielle Stone
# Topic: Initial Segregation Code
# File Created On  Zelner and Stephanie Choi Work

###* Call in Libraries ###

library(deSolve)
library(ggplot2)
library(plotly)
library(ppcor)
library(lhs)

###* Call in Libraries ###

#test

########################################################################
############################# TODO: ODE ################################
########################################################################

################### Model Function #################
two_group_sir_model <- function(t, x, params) {
  p1 <- 0.49
  p2  <- 0.51
  #### Parameters ####
  beta <- params[1] # Infections, Case/Day
  gamma <- 14 # Duration of Infectiousness, Days
  tau <- params[2] # Minority Group Degree of Isolation
  rho <- params[3] # Minority Group Degree of Concentration


  #### x[i] #### nolint
  ## Minority Group
  S1 <- x[1]
  I1 <- x[2]
  R1 <- x[3]
  # p1 <- sum(x[1:3]) # Minority Proportion of Total Population
  ## Majority Group
  S2 <- x[4]
  I2 <- x[5]
  R2 <- x[6]
  # p2 <- 1 - p1 # Majority Proportion of Total Population


  #### Contact (Between/Within) #### nolint
  g1_within <- tau / p1
  g1_between <- (1 - tau) / p2
  g2_within <- (1 - (p1 * g1_between)) / p2
  g2_between <- (p1 * g1_between) / p1


  #### Infections ####
  lambda_1 <- rho * beta * S1 * (g1_within * I1 + g1_between * I2)
  lambda_2 <- beta * S2 * (g2_within * I2 + g2_between * I1)


  #### ODE Model ####
  ## Low SES
  dS1dt <- -(lambda_1)
  dI1dt <- (lambda_1) - (gamma * I1)
  dR1dt <- (gamma * I1)
  ## High SES
  dS2dt <- -(lambda_2)
  dI2dt <- (lambda_2) - (gamma) * (I2)
  dR2dt <- (gamma) * (I2)


  #### Results Vector ####
  dxdt <- c(dS1dt, dI1dt, dR1dt, dS2dt, dI2dt, dR2dt)


  #### Return Results Vector ####
  return(list(dxdt))
}



########################################################################
############################ TODO: Setup ###############################
########################################################################

###? Define Sweeping Parameters ########################################
params <- data.frame(
  "param" = c("beta", "tau", "rho"),
  "lower" = c(0.01, 0.0, 1.0),
  "upper" = c(1.0, 1.0, 1.0),
  "default" = c(0.1, 0.5, 1.0)
)

###? Sampling Setup ####################################################
num_samples <- 1000
num_params <- 3
#rand_param_sample <- matrix(
#                        runif(
#                          (num_samples * num_params), ncol = num_params))
lhs_param_sample <- randomLHS(num_samples, num_params)

###? Time Variables ####################################################
t_0 <- 1
t_end <- 365
t_step <- 1

times <- seq(t_0, t_end, by = t_step)

###? Initial Conditions ################################################
S1_0 <- 0.999 * p1_0
I1_0 <- 0.001 * p1_0
R1_0 <- 0.000
S2_0 <- 0.999 * p2_0
I2_0 <- 0.001 * p2_0
R2_0 <- 0.000

init_conditions <- c(S1_0, I1_0, R1_0, S2_0, I2_0, R2_0)

###? Create List of Data Frames to Fill ################################

# for (i in 1:num_samples) {
#   sweep_df <- data.frame(
#                 row.names = c("S1", "S2", "I1", "I2", "R1", "R2"))
#   sweep_df_list[i] <- sweep_df
# }

sweep_df_list <- list()

paramsample_scale <- matrix(NA, num_samples, num_params, TRUE)

# Scale our params to the min/max values
paramsample_scale <- (params$lower +
                        (params$upper - params$lower) *
                            lhs_param_sample)


########################################################################
############################ TODO: Sweeps ##############################
########################################################################

for(i in 1:num_samples) {

  # Temp Data Frame
  temp <- matrix(data = NA, nrow = length(times), ncol = 6)
  colnames(temp) <- c("S1", "S2", "I1", "I2", "R1", "R2")

  # Run the ODE with current params
  temp <- ode(y = init_conditions,
                    times = times,
                      func = two_group_sir_model,
                        parms = paramsample_scale[i, ],
                          method = lsode
            )

  append(sweep_df_list, c(temp))
  print(i)
}



########################################################################
########################### TODO: Visualize ############################
########################################################################
