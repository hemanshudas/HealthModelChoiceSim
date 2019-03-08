#### Markov Cycle Tree Simulation
#### Tue, 8 March 2019
#### The code is an implementation of the microsimulation model to emulate individuals through state-level transitions


rm(list = ls())               #removing any variables in R's memory

#Model Inputs
n.i <- 2000                 #Number of individuals
ini_age <- 40                 #Start Age of Cohort
end_age <- 100                #Simulation End Age
cl <- 1                       #Cycle Length in years
v.n <- c("W","E","U","D")     #Model States: W-Well, E-Event, U-Unwell, D-Dead
n.s <- length(v.n)            #Number of model states
v.M_1 <- rep("W",n.i)       #Everyone is Well at the start
v.Str <- c("Drug A","Drug B") #Storing the Strategy Indicator
n.t <- end_age - ini_age      #Simulation cycles

#Initialization of Transition Probabilities
p.WE <- 0.05                  #Probability of event when well
p.WU <- 0                 #Probability of unwell when well
p.WD <- 0.005                  #Probability of death when well
p.ED <- 0.4                   #Probability of death when event
##p.UD <- 0.3                   #Probability of death when unwell attributed to event
#Probability of death due to natural causes is same of Well state as well as Unwell State
#Event state is a transitory state, and p.EE = 0; Also, the patient never recovers to well state if once unwell or event

#Intiatlization of Cost for each state
c.W <- 0                      #cost of being well
c.E <- 5000                   #cost of event occuring
##c.U <- 2000                 #cost of staying unwell (treatment)
c.D <- 0                      #cost of being dead

#Intialization of Utility for each state
u.W <- 1                      #utility of being well
u.E <- 0.3                    #utility during event state
u.U <- 0.9                    #utility of staying unwell
u.D <- 0                      #utility of being dead

#Characteristics of Drug A
p.UD_A <- 0.02
c.U_A <- 2000

#Characteristics of Drug B
p.UD_B <- 0.4
c.U_B <- 1800

####Function to simulate the markov cycle tree
MicroSim <- function(p.UD,c.U,TS.out = TRUE,TR.out=TRUE,seed=1) { 
#Arguments:
          #TS.out:  Flag for matrix of transitions between states
          #TR.out:  Flag for microsimulation trace
          #p.UD:    Probability of death when unwell attributed to event
          #c.U:     Cost of staying unwell due to cost of drugs

          v.dwc <- 1 ^ (0:n.t)            #Vector Multiplication for calculating costs
          v.dwu <- 1 ^ (0:n.t)            #Vector Multiplication for calculating QALY

#Matrixes to capture state, cost and health outcomes for all individuals at any time point
m.M <- m.C <- m.U <- matrix(nrow = n.i, ncol=n.t+1,
                            dimnames = list(paste("ind",1:n.i,sep = " "),
                                            paste("cycle",0:n.t, sep = " ")))
m.M[,1] <- v.M_1              #Initial health state

for (i in 1:n.i) {
          set.seed(seed+i)
          m.C[i,1] <- Cost(m.M[i,1],c.U)              #Costs per individual for initial health state
          m.U[i,1] <- Util(m.M[i,1])              #QALY per individual for initial health state
          
          for (t in 1:n.t) {
                    v.p <- Probs(m.M[i,t],p.UD)        #Transition probabilities at cycle t
                    
                    m.M[i, t+1] <- sample(v.n, prob = v.p, size =1)   #Sample the next health state
                    m.C[i, t+1] <- Cost(m.M[i,t+1],c.U)               #Costs per individual at cycle t+1
                    m.U[i, t+1] <- Util(m.M[i,t+1])                   #QALY per individual at cycle t+1
          }                                       #Closing loop for cycles
          if (i/100 == round (i/100,0)) {
                    cat('\r', paste(i/n.i * 100, "% done", sep = " "))
          }                                       #Closing loop for progress display
}


tc <- m.C %*% v.dwc           #Calculating the Costs
tu <- m.U %*% v.dwu           #Calculating the QALYs

tc_hat <- mean(tc)            #Average Cost
tu_hat <- mean(tu)            #Average QALY

#Optional Matrix of Transitions between states
if (TS.out == TRUE) {
          TS <- paste(m.M, cbind(m.M[,-1],NA), sep = "->")            #Transitions from one state to another
          TS <- matrix(TS,nrow=n.i)
          rownames(TS) <- paste("Ind", 1:n.i, sep = " ")              #Naming the rows
          colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")            #Naming the columns
} else {
          TS <- NULL
}

#Optional Output Trace
if (TR.out == TRUE) {
          TR <- t(apply(m.M,2,function(x) table(factor(x,levels=v.n, ordered = TRUE))))
          TR <- TR/n.i                                                #Creating a Distribution Trace
          rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")            #Naming the rows
          colnames(TR) <- v.n                                         #Naming the columns
} else {
          TR <- NULL
}

results <- list(m.M = m.M, m.C = m.C, m.U = m.U, tc_hat = tc_hat, tu_hat = tu_hat, tu=tu, TS = TS, TR = TR)
return(results)
}                             #End of the MicroSim Function

####Probability Function to update the transition probabilities of every cycle
Probs <- function(M.it,p.UD) {
          #M.it:    Health state occupied by individual i at cycle t
          #p.UD:    Probability of death when unwell due to event
          
          v.p.it <- rep(NA,n.s)                                       #Vector of transition probabilities
          names(v.p.it) <- v.n                                        #Naming the vector
          
          v.p.it[M.it == "W"] <- c(1-(p.WE+p.WU+p.WD),p.WE,p.WU,p.WD) #tranistion probability when well
          v.p.it[M.it == "E"] <- c(0,0,1-(p.ED+p.WD),p.ED+p.WD)       #transition probability when event occurs
          v.p.it[M.it == "U"] <- c(0,0,1-(p.UD+p.WD),p.UD+p.WD)       #transition probability when unwell
          v.p.it[M.it == "D"] <- c(0,0,0,1)                           #transition probability when death
          return(v.p.it)                                              #returning probabilities
}

####Cost function estimates the cost at every cycle
Cost <- function (M.it, c.U) {
          #M.it:    Heath state occupied by individual i at cycle t
          #c.U:     Cost of unwell state due to continued medication
          
          c.it <- 0
          c.it[M.it == "W"] <- c.W                #Cost if well
          c.it[M.it == "E"] <- c.E                #Cost if event
          c.it[M.it == "U"] <- c.U                #Cost if unwell
          c.it[M.it == "D"] <- c.D                #Cost if dead
          return(c.it)                            #Returning costs
}

#Util function estimates the QALY at every cycle
Util <- function (M.it) {
          #M.it:    Health state occupied by individual i at cycle t
          
          u.it <- 0
          u.it[M.it == "W"] <- u.W                #QALY if well
          u.it[M.it == "E"] <- u.E                #QALY if event
          u.it[M.it == "U"] <- u.U                #QALY if unwell
          u.it[M.it == "D"] <- u.D                #QALY if dead
          QALY <- u.it * cl
          return(QALY)                            #Returning costs
}

#### Running the simulation
sim_drugA <- MicroSim(p.UD_A,c.U_A)
sim_drugB <- MicroSim(p.UD_B,c.U_B)

#Calculating the comparision between the drugs
v.C <- c(sim_drugA$tc_hat,sim_drugB$tc_hat)
v.U <- c(sim_drugA$tu_hat,sim_drugB$tu_hat)

ICER <- (v.C[2]-v.C[1])/(v.U[2]-v.U[1])                     #Incremental Cost Effectiveness Ratio

#Data Presentation
table_microsim <- data.frame(
          round(v.C, 0),              # costs per arm
          v.U,              # health outcomes per arm
          c("", round(ICER, 3))       # ICER
)
rownames(table_microsim) = v.Str  # name the rows
colnames(table_microsim) = c("Costs", "QALYs", "ICER") # name the columns
table_microsim                    # print the table 