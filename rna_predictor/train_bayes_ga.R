#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
  make_option(c("--output_parameters"), type="character",default="bayesian_ga_model.txt",
              help="identification code [default %default]"),
  make_option(c("--mcmc_cycles"), type="integer", default=20000,
              help="number of mcmc optimization cycles [default %default]"),
  make_option(c("--ga_cycles"), type="integer", default=100,
              help="number of GA optimization cycles [default %default]"),
  make_option(c("--ga_populations"), type="integer", default=25,
              help="number of GA populations to evolve [default %default]"),
  make_option(c("--ga_std"), type="double", default=2.0,
              help="the number of standard deviation (in the mcmc) parameters that is to be used to limit the GA search space around the bayesian parameters [default %default]"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="print header and progress information [default %default]")
)

parser <- OptionParser(usage = "%prog [options] cs_training_database cs_decoy_database", option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

if(length(arguments$args) != 2) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  if (options$verbose){
    cat("Compute the Sum of Logirthims of Ranks (SLR)\n")
    cat("Author: Aaron T. Frank\n")
    cat(sprintf("%s\n",date()))
  }    
  # initialize variables
  traindatafile <- arguments$args[1]
  decoydatafile <- arguments$args[2]
  output_parameters <- options$output_parameters
  mcmc_cycles <- options$mcmc_cycles
  ga_cycles <- options$ga_cycles
  ga_std <- options$ga_std
  ga_populations <- options$ga_populations

  # load needed libraries
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("MCMCpack"))
  suppressPackageStartupMessages(library("GA"))
  source("train_library.R")
  
  # load training data
  train <- load_training_data(traindatafile,0)
  train_ga <- get_ga_train_data(train)
  trainy <- train_ga$trainy
  trainx <- train_ga$trainx
  
  # load decoy data
  decoy_data <- load_decoy_data(decoydatafile,0)
  y <- decoy_data$y
  X <- decoy_data$X
  INFO <- decoy_data$INFO
  
  # train preliminary bayesian model by minimizing the error between actual and measured chemical shifts
  bayesian_model <- runBayesian(train,mcmc_cycles)
  
  # refine preliminary bayesian model maximizing the ability to resolve native from non-native conformations in a decoy pool
  refined_bayesian_model <- runGA(trainy, y, trainx, X, bayesian_model, INFO, ga_cycles = ga_cycles, population_size = ga_populations, std = ga_std)
  
  # output model to file
  write.table(refined_bayesian_model, file = output_parameters, col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(options$verbose)
    print(refined_bayesian_model)
}
