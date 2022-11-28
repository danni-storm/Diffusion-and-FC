# Danielle Stone
# Initial Segregation Code
# 10/31/2022


###* Call in Libraries ###

library(deSolve)
library(ggplot2)

###* Call in Libraries ###


###* Define Non-Sweeping Variables ###

#### Time ####
t_0 <- 1
t_end <- 365
t_step <- 1

#### Initial Conditions ####
S1_0 <- 1
S2_0 <- 1
prot_0 <- 1
I_0 <- 1
R_0 <- 1

###* Define Non-Sweeping Variables ###




########################################################################
############################# TODO: ODE ################################
########################################################################

################### Model #################
basic_SES_model <- function(t, x, params) {

    #### Parameters ####
    beta <- params[1] # Rate of transmission
    rho <- params[2] # Reduction of risk by protection
    tau <- params[3] # Increased rate of high SES access
    eta <- params[4] # Rate of access (>= 1) modifier
    gamma <- params[5] # Recovery Rate


    #### x[i] ####
    S1 <- x[1] # High SES
    S2 <- x[2] # Low SES
    prot <- x[3] # Protected
    I <- x[4] # Infected
    R <- x[5] # Recovered

    #### ODE ####
    dxdt <- numeric(length(x))
    
    # High SES
    dxdt[1] <- -(tau * eta * S1) - (beta * I * S1)
    # Low SES
    dxdt[2] <- -(eta * S2) - (beta * I * S2)
    # Protected
    dxdt[3] <- (tau * eta * S1) + (eta * S2) - (rho * beta * I * S2)
    # Infected
    dxdt[4] <- (rho * beta * I * S2) + (beta * I * S1) + (beta * I * S2) - (gamma * I)
    # Recovered
    dxdt[5] <- (gamma * I)


    return(list(dxdt))
}

########## Non-Sweeping Variables ###########
x0 <- c(S1_0, S2_0, prot_0, I_0, R_0)
times <- c(t_0, t_end, t_step)

########################################################################
############################ TODO: Sweeps ##############################
########################################################################