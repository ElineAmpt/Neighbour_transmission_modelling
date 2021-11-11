### Disease Transmission Models Ampt et al 2022 New Phytologist####

####################
# This script performs the model fitting for all models as described in Ampt et al 2022; New Phytologist. 
# Model 1-5 are fitted with an iterative approach and models 1-4 and 6-8 are fitted with an numerical approach.  
# Experimental disease transmission data contains metadata (column 1:7) and measurements of the furthest infected plant in the row on each timepoint (0-17days) (column 8:16)
# Each row in the data represents a single replicate in the experiment. 
# Replicates where the inoculated (i.e. 1st) seedling did not show disease symptoms were removed from this data. 
####################

#### General prep ####

# Load experimental disease transmission data

setwd('') 
data.raw <- read.csv("Experimental_transmission_Data_Ampt_et_al.csv", header = T,sep=";",check.names=F)

#create vector of row types (i.e. experimental treatments)
row.vector<-unique(data.raw$Row) 
row.vector<-sort(row.vector,decreasing=FALSE)

#set output working directory to save results of searches
setwd('') 

set.seed(100)

#Model evaluation per row type ####

for(t in 1:length(row.vector)){ #loop through row types
  
  #Prepare experimental data of single row type
  row.index <- which(data.raw$Row ==  
                       row.vector[t])
  data.subset <- data.raw[row.index,]# take one row type per loop iteration
  data.exp<-data.subset[,-c(1:7)]  # remove all metadata
  cn<-as.integer(colnames(data.exp)) # extract day numbers (column names) 
  data.exp<-rbind(cn,data.exp) #add day numbers as top row
  
  # Reformat experimental data for assessing fit of model predictions
  day.exp = data.exp[1,]
  num.day.exp = length(day.exp) 	# Get the dimensions for an array for experimental data
  num.day.tot = max(day.exp) 
  data.exp = data.exp[-1,] 	# Take everything except the top row with the days
  plant.num = 7 			# The total number of focal plants in a row
  exp.reps = nrow(data.exp) 
  
  overview.exp = array(NA, dim=c(plant.num,num.day.exp))
  for(i in 1:plant.num) { for(j in 1:num.day.exp) {
    overview.exp[i,j] = length(which(data.exp[,j] == i))
  } } # end of i and j loops
  
  # Iterative evaluation ####
  
  #Model 1 - Iterative evaluation #### 
  
  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain 
  model.parameters = 1	# The total number of free model parameters to fit.
  par.unit <- c(0.01)  	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 50	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      new.par[choose.par] = model.par[choose.par] + mutate.up.or.down*par.unit[choose.par]
      
      # Limit the range of log transmission probability to <= 0.
      if(new.par[1] > 0) new.par[1] = 0 
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      
      # Run the simulation model.
      data.sim = array(NA, dim=c(sim.num, num.day.tot)) # Array to store simulation data. 
      
      # Loop to run the model
      for(i in 1:sim.num) {
        p.trans.day=p.trans.start
        plant.inf = 1	# Per definition the first plant is infected.
        for(j in 1:num.day.tot) {
          if(j>2) new.infect = rbinom(n = 1, size = 1, prob = p.trans.day) else new.infect=0 # Start of transmission on day 3.
          if(new.infect == 1 & plant.inf < plant.num) plant.inf = plant.inf + 1
          data.sim[i,j] = plant.inf
        } # end of j loop
      } # end of i loop
      
      # Make an overview of the simulation results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.sim = array(NA, dim=c(plant.num,num.day.exp))
      overview.sim[,1] = 0
      overview.sim[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.sim[i,j] = 0 
        overview.sim[i,j] = (length(which(data.sim[,day.is.now] == i))+1)/(sim.num+2)   
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.sim[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up simulation data
      rm(data.sim)
      rm(overview.sim)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=2)
  write.table(final.data2, file=paste0("Iter Model1 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #clean up before starting new modeltype evaluation
  rm(final.data2)
  
  
  #Model 2 - Iterative evaluation #### 
  
  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain 
  model.parameters = 2	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.1) 	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[3])) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 65	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Psi") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if(choose.par==2&&new.par[2]>1) {par.unit1=par.unit[3]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] = model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameters 
      if(new.par[1] > 0) new.par[1] = 0     #to <= 0 (log scale)
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < 0) new.par[2] = 0  #to => 0
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      psi = new.par[2]
      
      # Run the simulation model.
      data.sim = array(NA, dim=c(sim.num, num.day.tot)) # Array to store simulation data. 
      
      # Loop to run the model
      for(i in 1:sim.num) {
        plant.inf = 1        # Per definition the first plant is infected.
        for(j in 1:num.day.tot) {
          if(j>2) p.trans.day=p.trans.start^exp(-psi*(plant.inf)) # Evaluate p.trans every day, start of transmission on day 3.
          if(j>2) new.infect = rbinom(n = 1, size = 1, prob = p.trans.day) else new.infect=0
          if(new.infect == 1 && plant.inf < plant.num) plant.inf = plant.inf + 1
          data.sim[i,j] = plant.inf
        } # end of j loop
      } # end of i loop
      
      # Make an overview of the simulation results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.sim = array(NA, dim=c(plant.num,num.day.exp))
      overview.sim[,1] = 0
      overview.sim[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.sim[i,j] = 0 
        overview.sim[i,j] = (length(which(data.sim[,day.is.now] == i))+1)/(sim.num+2)   
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.sim[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up simulation data
      rm(data.sim)
      rm(overview.sim)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=10)
  final.data2$Psi<-round((final.data2[,5]),digits=5)
  write.table(final.data2, file=paste0("Iter Model2 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #clean up before starting new modeltype evaluation
  rm(final.data2)
  
  #Model 3 - Iterative evaluation #### 
  
  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain 
  model.parameters = 2	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.1)   	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1])  #start transmission probability on a log scale
  par.range.2 <- c(seq(-10,-1,by=par.unit[3]),seq(-0.99,-0.01, by = par.unit[2]),seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[3])) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 65	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Gamma") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # Stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if(choose.par==2&&(new.par[2]>1||new.par[2]< -1)) {par.unit1=par.unit[3]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameter values.
      if(new.par[1] > 0) new.par[1] = 0   #to <= 0   (log scale)
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < -10) new.par[2] = -10  #to => 0
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      gamma=new.par[2]
      
      # Run the simulation model.
      data.sim = array(NA, dim=c(sim.num, num.day.tot)) # Array to store simulation data. 
      
      # Loop to run the model
      for(i in 1:sim.num) {
        plant.inf = 1        # Per definition the first plant is infected.
        for(j in 1:num.day.tot) {
          if(j>2) p.trans.day = p.trans.start^exp(-gamma*(j))  # Evaluate p.trans every day, start of transmission on day 3.
          if(j>2) new.infect = rbinom(n = 1, size = 1, prob = p.trans.day) else new.infect=0
          if(new.infect == 1 && plant.inf < plant.num) plant.inf = plant.inf + 1
          data.sim[i,j] = plant.inf
          
        } # end of j loop
      } # end of i loop
      
      # Make an overview of the simulation results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.sim = array(NA, dim=c(plant.num,num.day.exp))
      overview.sim[,1] = 0
      overview.sim[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.sim[i,j] = 0 
        overview.sim[i,j] = (length(which(data.sim[,day.is.now] == i))+1)/(sim.num+2)   
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.sim[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up simulation data
      rm(data.sim)
      rm(overview.sim)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=10)
  final.data2$Gamma<-round((final.data2[,5]),digits=5)
  write.table(final.data2, file=paste0("Iter Model3 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #clean up before starting new modeltype evaluation
  rm(final.data2)
  
  
  #Model 4 - Iterative evaluation #### 
  
  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain 
  model.parameters = 3	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.01,0.1)    	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[4]))
  par.range.3 <- c(seq(-10,-1,by=par.unit[4]),seq(-0.99,-0.01, by = par.unit[3]),seq(0, 0.99, by = par.unit[3]),seq(1,10,by=par.unit[4])) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1200		# The maximum number of steps to try in the search
  max.iter.stuck = 80	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Psi","Gamma") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1), sample(par.range.3, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if((choose.par==2&&new.par[2]>1)||(choose.par==3&&(new.par[3]>1||new.par[3]< -1))) {par.unit1=par.unit[4]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameter values
      if(new.par[1] > 0) new.par[1] = 0   #to <= 0   (log scale)   #(log scale)                    
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < 0) new.par[2] = 0  #to => 0
      if(new.par[3] > 10) new.par[3] = 10  #to <= 10
      if(new.par[3] < -10) new.par[3] = -10  #to => -10
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      psi = new.par[2]
      gamma=new.par[3]
      
      # Run the simulation model.
      data.sim = array(NA, dim=c(sim.num, num.day.tot)) # Array to store simulation data. 
      
      # Loop to run the model
      for(i in 1:sim.num) {
        plant.inf = 1        # Per definition the first plant is infected.
        for(j in 1:num.day.tot) {
          if(j>2)p.trans.day = p.trans.start^exp(-psi*(plant.inf)-gamma*(j)) # Evaluate p.trans every day, start of transmission on day 3.
          if(j>2) new.infect = rbinom(n = 1, size = 1, prob = p.trans.day) else new.infect=0
          if(new.infect == 1 && plant.inf < plant.num) plant.inf = plant.inf + 1
          data.sim[i,j] = plant.inf
        } # end of j loop
      } # end of i loop
      
      # Make an overview of the simulation results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.sim = array(NA, dim=c(plant.num,num.day.exp))
      overview.sim[,1] = 0
      overview.sim[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.sim[i,j] = 0 
        overview.sim[i,j] = (length(which(data.sim[,day.is.now] == i))+1)/(sim.num+2)   
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.sim[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up simulation data
      rm(data.sim)
      rm(overview.sim)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=10)
  final.data2$Psi<-round((final.data2[,5]),digits=5)
  final.data2$Gamma<-round((final.data2[,6]),digits=5)
  write.table(final.data2, file=paste0("Iter Model4 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #Model 5 - Iterative evaluation #### 
  
  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain 
  model.parameters = 2	# The total number of free model parameters to fit.
  par.unit <- c(0.1, 0.1) 	# The step size used when mutating the parameter. 
  par.range.1 <- seq(0.1, 2, by = par.unit[1]) # Range of values to be used for randomly picking a starting point of the search
  par.range.2 <- seq(0.1, 2, by = par.unit[2])
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 50	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "Alpha", "Beta") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1)) 
    
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      new.par[choose.par] = model.par[choose.par] + mutate.up.or.down*par.unit[choose.par]
      
      # Limit the range of parameter values
      if(new.par[1] < 0) new.par[1] = par.unit[1] 
      if(new.par[2] < 0) new.par[2] = par.unit[2] 
      
      # Antilog transform model parameters and rename
      alpha = 10^new.par[1]
      beta = 10^new.par[2]
      
      # Run the simulation model.
      data.sim = array(NA, dim=c(sim.num, num.day.tot)) # Array to store simulation data. 
      
      # Loop to run the model
      for(i in 1:sim.num) {
        # Draw the probability values to represent the susceptibility of plants in the row
        p.vector = rbeta(n = (plant.num), shape1 = alpha, shape2 = beta) #include 7 plants to avoid rbinom warnings when plant.inf==7; (transmission stops at this point) 
        plant.inf = 1	# Per definition the first plant is infected.
        for(j in 1:num.day.tot) {
          p.trans.now = p.vector[plant.inf] #pick probability next plant - p.vector[1] represents susceptibility of 2nd plant in the row (i.e. when plant.inf=1) and so on
          if(j>2) new.infect = rbinom(n = 1, size = 1, prob = p.trans.now) else new.infect=0
          if(new.infect == 1 & plant.inf < plant.num) plant.inf = plant.inf + 1
          data.sim[i,j] = plant.inf
        } # end of j loop
      } # end of i loop
      
      # Make an overview of the simulation results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.sim = array(NA, dim=c(plant.num,num.day.exp))
      overview.sim[,1] = 0
      overview.sim[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.sim[i,j] = 0 
        overview.sim[i,j] = (length(which(data.sim[,day.is.now] == i))+1)/(sim.num+2)   
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.sim[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up simulation data
      rm(data.sim)
      rm(overview.sim)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL=round(final.data2[,3],digits=2)
  final.data2$Alpha=round((10^final.data2[,4]),digits=2)
  final.data2$Beta=round((10^final.data2[,5]),digits=2)
  write.table(final.data2, file=paste0("Iter Model5 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  
  #Numerical evaluation####
  #Model 1 - Numerical evaluation #### 
  
  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 1	# The total number of free model parameters to fit.
  par.unit <- c(0.01)  	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 50	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      new.par[choose.par] = model.par[choose.par] + mutate.up.or.down*par.unit[choose.par]
      
      # Limit the range of log transmission probability to <= 0.
      if(new.par[1] > 0) new.par[1] = 0 
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)] - p.trans.day*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.day*overview.mod[(r-1),(q-1)] - p.trans.day*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.day*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=2)
  write.table(final.data2, file=paste0("Num Model1 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  
  #Model 2 - Numerical evaluation #### 

  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 2	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.1) 	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[3])) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 65	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Psi") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if(choose.par==2&&new.par[2]>1) {par.unit1=par.unit[3]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] = model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameters 
      if(new.par[1] > 0) new.par[1] = 0     #to <= 0 (log scale)
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < 0) new.par[2] = 0  #to => 0
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      psi = new.par[2]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          p.trans.before = p.trans.day^exp(-psi*(r-1))
          p.trans.now = p.trans.day^exp(-psi*(r))
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)]                                            - p.trans.now*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=10)
  final.data2$Psi<-round((final.data2[,5]),digits=5)
  write.table(final.data2, file=paste0("Num Model2 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #Model 3 - Numerical evaluation #### 

  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain  #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 2	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.1)   	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1])  #start transmission probability on a log scale
  par.range.2 <- c(seq(-10,-1,by=par.unit[3]),seq(-0.99,-0.01, by = par.unit[2]),seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[3])) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 65	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Gamma") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # Stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if(choose.par==2&&(new.par[2]>1||new.par[2]< -1)) {par.unit1=par.unit[3]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameter values.
      if(new.par[1] > 0) new.par[1] = 0   #to <= 0   (log scale)
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < -10) new.par[2] = -10  #to => 0
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      gamma=new.par[2]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          p.trans.now = p.trans.day^exp(-gamma*(q)) # Could be evalutad in q loop, but for consistency between models done here
          p.trans.before = p.trans.now # Not necessary here, but for consistency has been kept in.
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)]                                            - p.trans.now*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=10)
  final.data2$Gamma<-round((final.data2[,5]),digits=5)
  write.table(final.data2, file=paste0("Num Model3 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #Model 4 - Numerical evaluation #### 

  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain  #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 3	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.01,0.1)    	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[4]))
  par.range.3 <- c(seq(-10,-1,by=par.unit[4]),seq(-0.99,-0.01, by = par.unit[3]),seq(0, 0.99, by = par.unit[3]),seq(1,10,by=par.unit[4])) # Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1200		# The maximum number of steps to try in the search
  max.iter.stuck = 80	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Psi","Gamma") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1), sample(par.range.3, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if((choose.par==2&&new.par[2]>1)||(choose.par==3&&(new.par[3]>1||new.par[3]< -1))) {par.unit1=par.unit[4]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameter values
      if(new.par[1] > 0) new.par[1] = 0   #to <= 0   (log scale)   #(log scale)           
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < 0) new.par[2] = 0  #to => 0
      if(new.par[3] > 10) new.par[3] = 10  #to <= 10
      if(new.par[3] < -10) new.par[3] = -10  #to => -10
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      psi = new.par[2]
      gamma=new.par[3]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          gamma.now = gamma
          psi.now = psi
          p.trans.before = p.trans.day^exp(-psi.now*(r-1)-gamma.now*(q))
          p.trans.now = p.trans.day^exp(-psi.now*(r)-gamma.now*(q))
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=10)
  final.data2$Psi<-round((final.data2[,5]),digits=5)
  final.data2$Gamma<-round((final.data2[,6]),digits=5)
  write.table(final.data2, file=paste0("Num Model4 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #Model 6 - Numerical evaluation #### 

  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 3	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,1,0.1)  	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[4]))
  par.range.3 <- seq(0, 5, by = par.unit[3])# Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 65	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Psi","omega") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1), sample(par.range.3, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if(choose.par==2&&new.par[2]>1) {par.unit1=par.unit[4]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of log transmission probability to <= 0.
      if(new.par[1] > 0) new.par[1] = 0      #to <= 0  (log scale)                  
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < 0) new.par[2] = 0  #to => 0
      if(new.par[3] > 5) new.par[3] = 5 #to <= 5
      if(new.par[3] < 0) new.par[3] = 0 #to => 0
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      psi = new.par[2]
      omega=new.par[3]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          if((r-1) <= omega) p.trans.before = p.trans.day else p.trans.before = p.trans.day^exp(-psi*(r-1-omega))
          if(r <= omega)     p.trans.now = p.trans.day else p.trans.now = p.trans.day^exp(-psi*(r-omega))
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=2)
  final.data2$Psi<-round((final.data2[,5]),digits=2)
  final.data2$omega<-round(final.data2[,6],digits=2)
  write.table(final.data2, file=paste0("Num Model6 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #Model 7 - Numerical evaluation #### 

  # Parameters for model fitting
  sim.num = 10000  # The number of times to simulate a transmission chain  #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 3	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,1,0.1)  	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(-10,-1,by=par.unit[4]),seq(-0.99,-0.01, by = par.unit[2]),seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[4]))
  par.range.3 <- seq(2, 16, by = par.unit[3])# Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1000		# The maximum number of steps to try in the search
  max.iter.stuck = 65	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Gamma","Delta") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1), sample(par.range.3, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if(choose.par==2&&(new.par[2]>1||new.par[2]< -1)) {par.unit1=par.unit[4]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameter values.
      if(new.par[1] > 0) new.par[1] = 0 #to <= 0 (log scale)
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < -10) new.par[2] = -10  #to => -10
      if(new.par[3] > 17) new.par[3] = 17  #to <= 17 
      if(new.par[3] < 2) new.par[3] = 2 #to => 2
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      gamma=new.par[2]
      Delta=new.par[3]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          if(q <= Delta ) p.trans.now = p.trans.day else p.trans.now = p.trans.day^exp(-gamma*(q-Delta)) # Could be evaluated in q loop, but for consistency between models done here
          p.trans.before = p.trans.now # Not necessary here, but for consistency has been kept in.
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=2)
  final.data2$Gamma<-round((final.data2[,5]),digits=2)
  final.data2$Delta<-round(final.data2[,6],digits=2)
  
  write.table(final.data2, file=paste0("Num Model7 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  #Model 8 - Numerical evaluation #### 

  # Parameters for model fitting
  sim.num = 10000   	# The number of times to simulate a transmission chain  #!! In numerical evaluation only used for Laplace Binomial Point Estimator
  model.parameters = 5	# The total number of free model parameters to fit.
  par.unit <- c(0.01, 0.01,0.01,1,1,0.1)   	# The step size used when mutating the parameter. 
  par.range.1 <- seq(-2, -0.1, by = par.unit[1]) 
  par.range.2 <- c(seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[6]))
  par.range.3 <- c(seq(-10,-1,by=par.unit[6]),seq(-0.99,-0.01, by = par.unit[2]),seq(0, 0.99, by = par.unit[2]),seq(1,10,by=par.unit[6]))
  par.range.4 <- seq(0, 5, by = par.unit[4])
  par.range.5 <- seq(2, 16, by = par.unit[5])# Range of values to be used for randomly picking a starting point in the search.
  
  # Search parameters
  max.iter = 1200		# The maximum number of steps to try in the search
  max.iter.stuck = 80	# The number of times to mutate parameters and try to improve fit before quitting the search.
  repetitions = 1750	# Number of times to repeat the stochastic hill climbing process
  
  # Generate an array to store final data; column 1 is the total iterations, 2 is stuck counter,
  # 3 is the NLL, and 4 is the model parameter estimate.
  final.data <- array(NA, dim=c(repetitions,(3+model.parameters)))
  colnames(final.data) <- c("Iterations", "Stuck_Counter", "NLL", "P_start", "Psi","Gamma","omega","Delta") 
  bs.counter = 1 		# A counter that is printed to show the search is still running.
  
  # stochastic hill climbing loop 
  for(a in 1:repetitions) {
    
    # Randomly pick start values for the search. 
    model.par = c(sample(par.range.1, size = 1), sample(par.range.2, size = 1), sample(par.range.3, size = 1), sample(par.range.4, size = 1), sample(par.range.5, size = 1)) 
    
    # Housekeeping: set iterations to 1, exit.flag indicates when to end the search, and set
    # NLL to a very high value so that the first value obtained in the search is always accepted.
    iter = 1
    exit.flag = 0
    nll = 1e15
    stuck.counter = 0
    
    while(exit.flag == 0) {
      
      # Choose model parameter to mutate and then mutate it
      new.par = model.par		
      choose.par = sample(1:model.parameters, size = 1)
      mutate.up.or.down = sample(c(1, -1), size = 1)
      if((choose.par==2&&new.par[2]>1)||(choose.par==3&&(new.par[3]>1||new.par[3]< -1))) {par.unit1=par.unit[6]} else {par.unit1=par.unit[choose.par]}
      new.par[choose.par] =  model.par[choose.par] + mutate.up.or.down*par.unit1
      
      # Limit the range of parameter values
      if(new.par[1] > 0) new.par[1] = 0 #to => 0 (log scale)
      if(new.par[2] > 10) new.par[2] = 10  #to <= 10
      if(new.par[2] < 0) new.par[2] = 0  #to => 0
      if(new.par[3] > 10) new.par[3] = 10  #to <= 10
      if(new.par[3] < -10) new.par[3] = -10  #to => -10
      if(new.par[4] > 5) new.par[4] = 5 #to <= 5
      if(new.par[4] < 0) new.par[4] = 0 #to => 0
      if(new.par[5] > 17) new.par[5] = 17#to <= 17
      if(new.par[5] < 2) new.par[5] = 2#to => 2
      
      # Antilog transform model parameter and rename
      p.trans.start = 10^new.par[1]
      psi = new.par[2]
      gamma=new.par[3]
      omega=new.par[4]
      Delta=new.par[5]
      
      # Generate an array for storing the results
      overview.mod = array(NA, dim=c(plant.num,num.day.tot))
      # Add results for the first and second day, when only the inoculated plant is infected by the pathogen.
      overview.mod[,1:2] = c(1, rep(0, (plant.num-1)))
      
      # Now calculate each frequency of the furthest infected plant of every day after the inoculation day.
      p.trans.day=p.trans.start
      for(q in 3:(num.day.tot)) {
        for(r in 1:plant.num) {
          if(q <= Delta ) gamma.now=  0 else gamma.now = gamma
          if((r-1) <= omega) psi.now = 0 else psi.now = psi
          p.trans.before = p.trans.day^exp(-psi.now*(r-1-omega)-gamma.now*(q-Delta))
          if(r <= omega)  psi.now = 0 else psi.now=psi   
          p.trans.now = p.trans.day^exp(-psi.now*(r-omega)-gamma.now*(q-Delta))
          if(r==1) 		overview.mod[r,q] = overview.mod[r,(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r>1 & r<plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)] - p.trans.now*overview.mod[r,(q-1)]
          if(r==plant.num)	overview.mod[r,q] = overview.mod[r,(q-1)] + p.trans.before*overview.mod[(r-1),(q-1)]
        } # end of r loop
      } # End of q loop
      
      # Make an overview of the model results, as was done for the experimental data.
      # Note that only those days for which there are experimental data should be sent to
      # this array. Use Laplace Binomial Point Estimator for the predicted
      # model frequency. The first plant is assumed to be infected on day 0.
      
      overview.mod2 = array(NA, dim=c(plant.num,num.day.exp))
      overview.mod2[,1] = 0
      overview.mod2[1,1] = 1
      for(i in 1:plant.num) { for(j in 2:num.day.exp) {  # "j in 2..." because in col 1 only the first plant is assumed to be infected on day 0.
        day.is.now = as.numeric(day.exp[j])
        overview.mod2[i,j] = 0 
        overview.mod2[i,j] = ((overview.mod[i,day.is.now]*sim.num)+1)/(sim.num+2)
      } } # end of i and j loops
      
      # Compare the model to the experimental data using the NLL derived from the multinomial 
      # distribution for each day, then sum over days.
      nll.day = numeric(num.day.exp) # Vector to store data
      for(i in 1:num.day.exp) nll.day[i] = dmultinom(x = overview.exp[,i], prob = overview.mod2[,i], log = TRUE) # Likelihood calculation
      new.nll = -sum(nll.day)
      
      # Here check if the new parameter values lead to improved fit. If so, accept the new
      # parameter values and NLL, and reset the StuckCounter. Otherwise, repeat the process
      # with the same values used previously.
      stuck.counter = stuck.counter + 1
      
      if (new.nll < nll) {
        model.par = new.par
        nll = new.nll
        stuck.counter = 0
      }   # end of if statement
      
      # Check whether to continue the search or exit
      iter = iter + 1
      if (stuck.counter == max.iter.stuck) exit.flag = 1
      if (iter == max.iter) exit.flag = 1 
      
      #clean-up model data
      rm(overview.mod)
      rm(overview.mod2)
      
    }  # End of while loop for search  
    
    #save outcome of search
    final.data[a,] = c(iter, stuck.counter, nll, model.par)
    bs.counter = bs.counter + 1
    
  } # end of a loop 
  
  #save results of searches and write to file
  final.data2<-as.data.frame(final.data)
  final.data2=final.data2[order(final.data2[,3]),]
  final.data2$NLL<-round(final.data2[,3],digits=2)
  final.data2$P_start<-round((10^final.data2[,4]),digits=2)
  final.data2$Psi<-round((final.data2[,5]),digits=2)
  final.data2$Gamma<-round((final.data2[,6]),digits=2)
  final.data2$omega<-round(final.data2[,7],digits=2)
  final.data2$Delta<-round(final.data2[,8],digits=2)
  write.table(final.data2, file=paste0("Num Model8 ",row.vector[t], ".txt"),sep=";",dec=",")
  
  
  } # end of r loop
