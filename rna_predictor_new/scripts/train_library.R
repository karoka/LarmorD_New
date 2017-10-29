SLR <- function(X){
  # Helper function to compute the  Sum of Logarithms of Ranks (SLR) metric [Early recognition Metric]
  # See: Comprehensive Comparison of Ligand-Based Virtual Screening Tools Against the DUD Data set Reveals Limitations of Current 3D Methods
  # Args:
  #   X: input sorted data-frame with column named 'status' that takes values of 0 (inactive) and 1 (active)
  # Returns:
  #   returns: SLR score
  
  ri <- which(X$status==1)
  N <- nrow(X)
  i <-  1:length(ri)
  SLRmax <- -sum(log(i/N))
  return(-sum(log(ri/N))/SLRmax)
} 

mae <- function(x,y){
  # Helper function to compute mean absolute error by two vectors
  # Args:
  #   x: numeric vector 1
  #   y: numeric vector 2
  # Returns:
  #   returns: SLR score
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y)) 
    stop("'x' and 'y' must have the same number of elements")
  mean(abs(x-y))
}

popd <- function(X){
  # Helper function to remove the last row from a data frame or matrix
  # Args:
  #   X: data frame
  # Returns:
  #   returns: x without the last row
  X[-nrow(X),] 
}

popv <- function(x){
  # Helper function to remove the last element in a vector
  # Args:
  #   x: vector
  # Returns:
  #   returns: x without the element
  x[-length(x)] 
}

compute_chemical_shifts <- function(w,X){
  # Helper function to compute the chemical shifts
  # Args:
  #   w: weights -- vector 
  #   X: data matrix of features -- data matrix
  # Returns:
  #   returns: the computed chemical shifts shifts as weighted sum
  MX = nrow(X)
  NX = ncol(X)
  .rowSums(matrix(c(w),byrow=T,nrow=MX,ncol=NX)*X,MX,NX)
}

compute_mae_score <- function(w,y,X){
  # Helper function to compute error score
  # Args:
  #   w: weights -- vector
  #   y: reference shifts -- vector
  #   X: data matrix of features -- data matrix
  # Returns:
  #   returns: the computed chemical shifts shifts as weighted sum -- double
  yy <- compute_chemical_shifts(w,X)
  1/mae(y,yy)
}

compute_slr_score <- function(w, y, X, INFO, rmsd_threshold = 2.5, testing = FALSE, key = "Training"){
  # Helper function to compute SLR score
  # Args:
  #   w: weights -- vector
  #   trainy: actual chemical shifts for
  #   y: actual chemical shifts for models in the conformational decoy pool -- vector
  #   X: structure features used to compute chemical shifts for models in conformational decoy pool -- data matrix
  #   INFO: information about model in the decoy pools -- data frame
  #   rmsd_thresold: threshold value used to designate native and non-native decoys -- double
  #   testing: whether to run in validation mode (if yes, return the SLR data over individual systems) -- bool
  # Returns:
  #   returns: the mean SLR obtained when using the MAE between actual and computed chemical shifts to separate models in the decoy pools
  INFO$expCS <- y
  INFO$predCS <- compute_chemical_shifts(w,data.matrix(X))
  results <- ddply(.dat=INFO,.var=c("ID","frame","rmsd"),.fun=function(x){mae(x$predCS,x$expCS)})
  results <- results[order(results$V1),]
  results$status <- 0
  results$status[results$rmsd<rmsd_threshold] <- 1
  if(testing){
    results <- ddply(.dat=results,.var=c("ID"),.fun=function(x){summarize_rmsd_slr(x)})
    mean <- data.frame(ID = "MEAN", SLR = mean(results$SLR), min_rmsd = mean(results$min_rmsd), mean_rmsd = mean(results$mean_rmsd) )
    results <- rbind(results, mean)
    results$KEY <- key
    results
  } else {
    results <- ddply(.dat=results,.var=c("ID"),.fun=function(x){SLR(x)})
    mean(results$V1)  
  }
}

summarize_rmsd_slr <- function(x){
  # Helper function that report statisitics summarizing the ability to resolution native and non-native decoys
  # Args:
  #   x: comformational decoy data -- data frame
  # Returns:
  #   returns:  SLR, minimum rmsd, and mean rmsd -- data frame
  mean_rmsd <- mean(head(x,5)$rmsd)
  min_rmsd <- mean(head(x,1)$rmsd)
  SLR <- SLR(x)
  return(data.frame(SLR=round(SLR,3), min_rmsd = round(min_rmsd,3), mean_rmsd = round(mean_rmsd,3)))
}

fitness <- function(w,trainy,trainx,y,X,INFO,rmsd_threshold=2.5){
  # Helper function to the overall fitness for GA
  # Args:
  #   w: weights -- vector
  #   trainy: actual chemical shifts for the training set
  #   trainx: structure features used to compute chemical shifts for the training set -- data matrix
  #   y: actual chemical shifts for models in the conformational decoy pool -- vector
  #   X: structure features used to compute chemical shifts for models in conformational decoy pool -- data matrix
  #   INFO: information about model in the decoy pools -- data frame
  #   rmsd_threshold: threshold value used to designate native and non-native decoys (passed to compute_slr_score) -- double
  # Returns:
  #   returns: the sum of the  of the MAE and SLR scores 
  #compute_slr_score(w,y,X,INFO,rmsd_threshold) + compute_mae_score(w,trainy,trainx)
  if(is.null(INFO)){
  	compute_mae_score(w,trainy,trainx)
  } else {
  	compute_slr_score(w,y,X,INFO,rmsd_threshold)
  }
}

load_training_data <- function(datafile,skip=0){
  # Helper function to load training data used to derive parameters for model
  # Args:
  #   datafile: path to file name -- character string
  #   filter_nucleus: nucleus type to retain -- character string
  #   skip: number of lines to skip when reading data -- integer
  # Returns:
  #   returns: data for the file -- data frame
  names <- c("expCS","randCS","resnameG","resnameA","resnameC","resnameU","ringCurrent","GUA.C1p","GUA.C2p","GUA.C3p","GUA.C4p","GUA.C5p","GUA.P","GUA.O5p","GUA.O3p","GUA.C2","GUA.C4","GUA.C5","GUA.C6","GUA.C8","GUA.N1","GUA.N2","GUA.N3","GUA.N7","GUA.N9","GUA.O6","ADE.C1p","ADE.C2p","ADE.C3p","ADE.C4p","ADE.C5p","ADE.P","ADE.O5p","ADE.O3p","ADE.C2","ADE.C4","ADE.C5","ADE.C6","ADE.C8","ADE.N1","ADE.N3","ADE.N6","ADE.N7","ADE.N9","URA.C1p","URA.C2p","URA.C3p","URA.C4p","URA.C5p","URA.P","URA.O5p","URA.O3p","URA.C2","URA.C4","URA.C5","URA.C6","URA.N1","URA.N3","URA.O4","CYT.C1p","CYT.C2p","CYT.C3p","CYT.C4p","CYT.C5p","CYT.P","CYT.O5p","CYT.O3p","CYT.C2","CYT.C4","CYT.C5","CYT.C6","CYT.N1","CYT.N3","CYT.N4","CYT.O2")
  cnames <- c("frame","ID","resname","resid","nucleus",names)
  tmp <- read.table(datafile,skip=skip,col.names=cnames)[,c(-1:-5)]
  tmp2 <- matrix(as.numeric(as.matrix(tmp)),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
  colnames(tmp2) <- names
  tmp2 <- as.data.frame(tmp2)
  tmp2$expCS <- tmp2$expCS - tmp2$randCS
  return(tmp2)
}

get_ga_train_data <- function(train){
  # Helper function to generate the train data need for the GA
  # Args:
  #   train: training data initially used for Bayesian modeling  -- data frame
  # Returns:
  #   returns: the actual chemical shifts for the training data (--vector) and structure features for the training data (-- data frame)  -- list
  return(list(trainy=train[,1],trainx=data.matrix(train[,c(-1,-2)])))
}

xnames <- c("resnameG","resnameA","resnameC","resnameU","ringCurrent","GUA.C1p","GUA.C2p","GUA.C3p","GUA.C4p","GUA.C5p","GUA.P","GUA.O5p","GUA.O3p","GUA.C2","GUA.C4","GUA.C5","GUA.C6","GUA.C8","GUA.N1","GUA.N2","GUA.N3","GUA.N7","GUA.N9","GUA.O6","ADE.C1p","ADE.C2p","ADE.C3p","ADE.C4p","ADE.C5p","ADE.P","ADE.O5p","ADE.O3p","ADE.C2","ADE.C4","ADE.C5","ADE.C6","ADE.C8","ADE.N1","ADE.N3","ADE.N6","ADE.N7","ADE.N9","URA.C1p","URA.C2p","URA.C3p","URA.C4p","URA.C5p","URA.P","URA.O5p","URA.O3p","URA.C2","URA.C4","URA.C5","URA.C6","URA.N1","URA.N3","URA.O4","CYT.C1p","CYT.C2p","CYT.C3p","CYT.C4p","CYT.C5p","CYT.P","CYT.O5p","CYT.O3p","CYT.C2","CYT.C4","CYT.C5","CYT.C6","CYT.N1","CYT.N3","CYT.N4","CYT.O2")

load_decoy_data <- function(datafile,skip=0){
  # Helper function to load training data used to derive parameters for model
  # Args:
  #   datafile: path to file name -- character string
  #   skip: number of lines to skip when reading data -- integer  
  # Returns:
  #   returns: data for the file -- data frame
  names <- c("expCS","randCS","resnameG","resnameA","resnameC","resnameU","ringCurrent","GUA.C1p","GUA.C2p","GUA.C3p","GUA.C4p","GUA.C5p","GUA.P","GUA.O5p","GUA.O3p","GUA.C2","GUA.C4","GUA.C5","GUA.C6","GUA.C8","GUA.N1","GUA.N2","GUA.N3","GUA.N7","GUA.N9","GUA.O6","ADE.C1p","ADE.C2p","ADE.C3p","ADE.C4p","ADE.C5p","ADE.P","ADE.O5p","ADE.O3p","ADE.C2","ADE.C4","ADE.C5","ADE.C6","ADE.C8","ADE.N1","ADE.N3","ADE.N6","ADE.N7","ADE.N9","URA.C1p","URA.C2p","URA.C3p","URA.C4p","URA.C5p","URA.P","URA.O5p","URA.O3p","URA.C2","URA.C4","URA.C5","URA.C6","URA.N1","URA.N3","URA.O4","CYT.C1p","CYT.C2p","CYT.C3p","CYT.C4p","CYT.C5p","CYT.P","CYT.O5p","CYT.O3p","CYT.C2","CYT.C4","CYT.C5","CYT.C6","CYT.N1","CYT.N3","CYT.N4","CYT.O2")
  cnames <- c("index","ID","frame","rmsd","resname","resid","nucleus",names)
  xnames <- c("resnameG","resnameA","resnameC","resnameU","ringCurrent","GUA.C1p","GUA.C2p","GUA.C3p","GUA.C4p","GUA.C5p","GUA.P","GUA.O5p","GUA.O3p","GUA.C2","GUA.C4","GUA.C5","GUA.C6","GUA.C8","GUA.N1","GUA.N2","GUA.N3","GUA.N7","GUA.N9","GUA.O6","ADE.C1p","ADE.C2p","ADE.C3p","ADE.C4p","ADE.C5p","ADE.P","ADE.O5p","ADE.O3p","ADE.C2","ADE.C4","ADE.C5","ADE.C6","ADE.C8","ADE.N1","ADE.N3","ADE.N6","ADE.N7","ADE.N9","URA.C1p","URA.C2p","URA.C3p","URA.C4p","URA.C5p","URA.P","URA.O5p","URA.O3p","URA.C2","URA.C4","URA.C5","URA.C6","URA.N1","URA.N3","URA.O4","CYT.C1p","CYT.C2p","CYT.C3p","CYT.C4p","CYT.C5p","CYT.P","CYT.O5p","CYT.O3p","CYT.C2","CYT.C4","CYT.C5","CYT.C6","CYT.N1","CYT.N3","CYT.N4","CYT.O2")
  shifts <- read.table(datafile,skip=skip,col.names=cnames)
  y <- shifts$expCS - shifts$randCS
  X <- shifts[,xnames]
  INFO <- shifts[,c("ID","frame","rmsd","nucleus")]
  return(list(INFO=INFO,X=X,y=y))
}

get_bayesian_parameters <- function(mcmc){
  # Helper function to retrive the model parameters from MCMC object
  # Args:
  #   mcmc: MCMC fitting object  -- MCMC object
  # Returns:
  #   returns: parameters -- data frame
  paras <- as.data.frame(summary(mcmc)$statistics[xnames,c("Mean","SD")])
  sigma <- sqrt(summary(mcmc)$statistics["sigma2","Mean"])[1]
  #print(paras)
  paras[,3] <- rownames(paras)
  colnames(paras) <- c("para","spread","name")
  rownames(paras) <- 1:nrow(paras)
  paras <- paras[,c("name","para","spread")]
  paras$error <- sigma
  return(paras)
}

runBayesian <- function(train, cycles, seed, formula=expCS~ringCurrent+resnameG+resnameA+resnameC+resnameU+GUA.C1p+GUA.C2p+GUA.C3p+GUA.C4p+GUA.C5p+GUA.P+GUA.O5p+GUA.O3p+GUA.C2+GUA.C4+GUA.C5+GUA.C6+GUA.C8+GUA.N1+GUA.N2+GUA.N3+GUA.N7+GUA.N9+GUA.O6+ADE.C1p+ADE.C2p+ADE.C3p+ADE.C4p+ADE.C5p+ADE.P+ADE.O5p+ADE.O3p+ADE.C2+ADE.C4+ADE.C5+ADE.C6+ADE.C8+ADE.N1+ADE.N3+ADE.N6+ADE.N7+ADE.N9+URA.C1p+URA.C2p+URA.C3p+URA.C4p+URA.C5p+URA.P+URA.O5p+URA.O3p+URA.C2+URA.C4+URA.C5+URA.C6+URA.N1+URA.N3+URA.O4+CYT.C1p+CYT.C2p+CYT.C3p+CYT.C4p+CYT.C5p+CYT.P+CYT.O5p+CYT.O3p+CYT.C2+CYT.C4+CYT.C5+CYT.C6+CYT.N1+CYT.N3+CYT.N4+CYT.O2+0 ){
  # Helper function that runs a Bayesian regression to general initial parameters for model
  # Args:
  #   train: training database -- data frame
  #   cycles: number of MCMC optimization cycles -- integer
  #   seed: random number used by MCMC optimization -- integer
  # Returns:
  #   returns: resulting model parameters -- data frame
  start <- rep(0,ncol(train)-2)
  start[1:5] <- 1 # make more general
  mcmc <- MCMCregress(formula, data = as.data.frame(train), verbose = TRUE, marginal.likelihood = "Chib95", b0 = 0, B0 = 0.1, c0 = 2, d0 = 0.11, beta.start = start, mcmc = cycles, burnin = 5000, seed = seed)
  get_bayesian_parameters(mcmc)
}

runGA <- function(trainy, trainx, initial_model, y=NULL, X=NULL, INFO=NULL, cycles = 100, population_size = 25, std = 2, rmsd_threshold = 2.5, seed = 12345, verbose = FALSE){
  # Helper function that runs GA regression to general refined parameters for model
  # Args:
  #   trainy: actual chemical shifts for training data -- vector
  #   y: actual chemical shifts for data in the decoy pools -- vector
  #   trainx: structure features for training data -- data frame
  #   X: structure features for decoy pools -- data frame
  #   bayesian_model: the initial Bayesian parameters -- data frame
  #   cycles (optional): optimization cycles -- integer
  #   populationSize (optional): population or number of solutions to evolve -- integer
  #   rmsd_threshold: threshold value used to designate native and non-native decoys (passed to compute_slr_score) -- double
  #   seed: random number seed -- double  
  # Returns:
  #   returns: resulting model parameters -- data frame
  ncomponents <- ncol(trainx)
  M <- nrow(trainx)
  N <- ncomponents
  lower <- rep(-1000, ncomponents) 
  upper <- rep(1000, ncomponents) 
  initialSolution <- matrix(initial_model$para, ncol = ncomponents, nrow = population_size, byrow = T)
  GA <- ga(type = "real", fitness = fitness, suggestions = initialSolution, min = lower, max = upper, seed = seed, popSize = population_size, maxiter = cycles, keepBest = TRUE, monitor = verbose, trainy = trainy, trainx = trainx, y = y, X = X, INFO = INFO, rmsd_threshold = rmsd_threshold)
  initial_model$ga <- GA@solution[1,]
  return(initial_model)
}