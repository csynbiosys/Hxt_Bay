files<-c("20170302Hxt1")


data_extraction_multiexperiment <- function (fileNamesVector){
  
  ######################## OBSERVABLES ########################
  
  mdp <- c() # Maximum data point for each file
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileObs <- paste(i, "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
    mdp = c(mdp, length(obs[,1]))
  }
  mrow <- max(mdp) # Maximum number of rows
  mcol <- length(fileNamesVector) # Maximum number of columns
  
  # Definition of the matrices with 0s on them
  samplingT <- matrix(data=0, nrow=mrow, ncol=mcol)
  GFPmean <- matrix(data=0, nrow=mrow, ncol=mcol)
  GFPstd <- matrix(data=0, nrow=mrow, ncol=mcol)

  
  for(i in 1:length(fileNamesVector)){ # Loop that opens a file for each iteration to extract the data
    fileObs <- paste(fileNamesVector[i], "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
      for(j in 1:length(obs[,1])){
        samplingT[j,i] <- round(obs[j,1])
        GFPmean[j,i] <- obs[j,2]
        GFPstd[j,i] <- obs[j,3]
      }
  }
  mdp2 <- matrix(data=mdp, nrow=1, ncol = mcol)
  
  ######################## INPUTS ########################
  
  mdpI <- c() # Maximum data point for each file
  mtimes <- c() # Maximum time point for each file
  
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileInp <- paste(i, "_Events_Inputs.csv", sep="")
    inp <- read.csv(file=fileInp, header=TRUE, sep=",")
    mdpI = c(mdpI, length(inp[,1]))
    mtimes = c(mtimes, round(inp[1,2]))
    
  }
  mrowI <- max(mdpI) # Maximum number of rows
  mcolI <- length(fileNamesVector) # Maximum number of columns
  mt <- max(mtimes) # Maximum number of time points
  
  mdpI2 <- matrix(data=mdpI, nrow=1, ncol = mcol)
  
  evnT <- matrix(data=0, nrow=mrowI+1, ncol=mcolI)
  u_Glu <- matrix(data=0, nrow=mrowI, ncol=mcolI)
  preGlu <- matrix(data=0, nrow=1, ncol=mcolI)#<<-inputs[1,3]
  time <- matrix(data=0, nrow=mt+1, ncol=mcolI)#<<- seq(1e-9, round(inputs[1,2]), length=round(inputs[1,2]))
  inps <- matrix(data=0, nrow=mrowI, ncol=mcolI)#<<-c()
  
  ltimes <- c() # Length of time series
  for(i in 1:length(fileNamesVector)){ # Loop that opens a file for each iteration to extract the data
    fileInp <- paste(fileNamesVector[i], "_Events_Inputs.csv", sep="")
    inp <- read.csv(file=fileInp, header=TRUE, sep=",")
    tempT <- seq(1e-9, round(inp[1,2]), length=round(inp[1,2])+1)
    ltimes = c(ltimes, length(tempT))
    ind = 1
    for(j in 1:length(inp[,1])){
      evnT[j,i] <- round(inp[j,1])
      u_Glu[j,i] <- inp[j,4]
      inps[j,i] <- inp[j,4]
    }
    
    for(l in 1:length(tempT)){
      
      time[l,i] = tempT[l]
    }
    evnT[(length(inp[,1])+1),i] = round(inp[1,2])
    preGlu[,i] <- inp[1,3]
     
  }
  ltimes2 <- matrix(data=ltimes, nrow=1, ncol = length(ltimes))
 
    toni <- seq(1e-9, 24*60) # Over night incuvation time
    
   data_multi <<- list (
          
          elm = mrowI, # Maximum length of the rows of the matrices except for time and evnT and pres
          tml = trunc(mt+1), # Maximum length of the rows for the time matrix
          ts = time+(1e-9),
          tsl = ltimes2, # length of time series per event
          tsmax = round(mtimes), # maximum time per event
          preGlu = preGlu,
          Glu = u_Glu,
          Nsp = (mdpI2+1), # length(evnT),
          inputs = inps+1e-7,
          evnT = round(evnT),
          m = mcol, # Number of time series
          stsl = mdp2, # Number of elements at each time series
          stslm = mrow,
          sts = round(samplingT),
          GFPmean = GFPmean,
          GFPstd = GFPstd,
          toni = toni,
          tonil = length(toni)
  )
}






