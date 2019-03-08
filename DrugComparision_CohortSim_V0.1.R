#### Markov Cycle Tree Simulation
#### Tue, 5 March 2019

rm(list = ls())               # removing any variables in R's memory

#Model Inputs
ini_age <- 40                 #Start Age of Cohort
end_age <- 100                #Simulation End Age
cl <- 1                       #Cycle Length in years
v.n <- c("W","E","U","D")     #Model States: W-Well, E-Event, U-Unwell, D-Dead
n.s <- length(v.n)            #Number of model states
v.M_1 <- c(2000,0,0,0)        #Everyone is Well at the start
v.Str <- c("Drug A","Drug B") #Storing the Strategy Indicator
n.t <- end_age - ini_age      #Simulation cycles

#Initialization of Transition Probabilities
p.WE <- 0.05                  #Probability of event when well
p.WU <- 0                     #Probability of unwell when well
p.WD <- 0.005                   #Probability of death when well
p.ED <- 0.4                   #Probability of unwell when event
##p.UD <- 0.3                   #Probability of death when unwell attributed to event
#Probability of death due to natural causes is same of Well state as well as Unwell State
#Event state is a transitory state, and p.EE = 0; Also, the patient never recovers to well state if once unwell or event

#Intiatlization of Cost for each state
c.W <- 0                      #cost of being well
c.E <- 3000                   #cost of event occuring
##c.U <- 2000                   #cost of staying unwell (treatment)
c.D <- 0                      #cost of being dead

#Intialization of Utility for each state
u.W <- 1                      #utility of being well
u.E <- 0.3                    #utility during event state
u.U <- 0.9                    #utility of staying unwell
u.D <- 0                      #utility of being dead

#Function to simulate the markov cycle tree
Markov <- function(v.M_1,n.t,p.UD,c.U) #REVISIT THIS
{ 
#Arguments:
          #v.m_1:   Initial allocation of cohort across states
          #n.t:     Total number of cycles to run the model
          #p.UD:    Probability of death when unwell attributed to event
          #c.U:     Cost of staying unwell due to cost of drugs

          v.dwc <- 1 ^ (0:n.t)            #Vector Multiplication for calculating costs
          v.dwu <- 1 ^ (0:n.t)            #Vector Multiplication for calculating QALY
                    
#Transition Probabibility Matrix
m.P <- matrix(c(1-p.WE-p.WU-p.WD,p.WE,p.WU,p.WD,
                0,0,1-(p.ED+p.WD),p.ED+p.WD,
                0,0,1-(p.UD+p.WD),p.UD+p.WD,
                0,0,0,1),
              nrow = n.s, ncol = n.s, byrow = T,
              dimnames = list(v.n,v.n))

#Transition Trace Matrix to capture the proportion of cohort in each state at each time point
m.TR <- matrix(0,nrow=n.t+1,ncol=n.s,
               dimnames = list(paste("cycle", 0:n.t, sep = ""), v.n))

m.TR[1,] <- v.M_1             #The initial health state

v.c <- c(c.W, c.E, c.U, c.D)
v.u <- c(u.W, u.E, u.U, u.D)

for (i in 2:(n.t+1))
{
          #calculating the proportion of cohort in each state at time t
          m.TR[i,] <- t(m.TR[i-1,]) %*% m.P
} #closing loop of single run of cohort

tc <- m.TR %*% v.c            #Calculating the Costs
tu <- m.TR %*% v.u            #Calculating the QALYs

t_tc <- t(tc) %*% v.dwc       #Total Cost
t_tu <- t(tu) %*% v.dwu       #Total QALY

results <- list(m.TR = m.TR,t_tc = t_tc, t_tu = t_tu)
return(results)
}

#### Running the simulation
p.UDA <- 0.02
c.UA <- 2000
sim_markov_drugA <- Markov(v.M_1,n.t,p.UDA,c.UA)

p.UDB <-0.4
c.UB <- 1800
sim_markov_drugB <- Markov(v.M_1,n.t,p.UDB,c.UB)

v.C <- c(sim_markov_drugA$t_tc,sim_markov_drugB$t_tc)       #Vector for capturing cost for both
v.U <- c(sim_markov_drugA$t_tu,sim_markov_drugB$t_tu)       #Vector for capturing QALY for both

ICER <- (v.C[2]-v.C[1])/(v.U[2]-v.U[1])                     #Incremental Cost Effectiveness Ratio

#Data Presentation
table_markov <- data.frame(
          round(v.C, 0),              # costs per arm
          round(v.U, 0),              # health outcomes per arm
          c("", round(ICER, 3))       # ICER
)
rownames(table_markov) = v.Str  # name the rows
colnames(table_markov) = c("Costs", "QALYs", "ICER") # name the columns
table_markov                    # print the table 