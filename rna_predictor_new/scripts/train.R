#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
  make_option(c("--output_parameters"), type="character",default="bayesian_ga_model.txt",
              help="identification code [default %default]"),
  make_option(c("--mcmc_cycles"), type="integer", default=20000,
              help="number of mcmc optimization cycles [default %default]"),
  make_option(c("--mcmc_seed"), type="integer", default=12345,
              help="random number seed used for MCMC optimization [default %default]"),
  make_option(c("--ga_cycles"), type="integer", default=100,
              help="number of GA optimization cycles [default %default]"),
  make_option(c("--ga_populations"), type="integer", default=10,
              help="number of GA populations to evolve [default %default]"),
  make_option(c("--ga_std"), type="double", default=3.0,
              help="the number of standard deviation (in the mcmc) parameters that is to be used to limit the GA search space around the bayesian parameters [default %default]"),
  make_option(c("--ga_rmsd_threshold"), type="double", default=2.5,
              help="threshold value used to designate native and non-native decoys [default %default]"),
  make_option(c("--ga_seed"), type="integer", default=12345,
              help="random number seed used for GA optimization [default %default]"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="print header and progress information [default %default]")
)

parser <- OptionParser(usage = "%prog [options] cs_training_database", option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

if(length(arguments$args) != 1) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  if (options$verbose){
    cat("Training LARMORD\n")
    cat("Author: Aaron T. Frank\n")
    cat("Author: Blair Whithington\n")
    cat(sprintf("%s\n",date()))
  }  
  # initialize variables
  traindatafile <- arguments$args[1]
  output_parameters <- options$output_parameters
  mcmc_cycles <- options$mcmc_cycles
  mcmc_seed <- options$mcmc_seed
  ga_cycles <- options$ga_cycles
  ga_std <- options$ga_std
  ga_populations <- options$ga_populations
  ga_rmsd_threshold <- options$ga_rmsd_threshold
  ga_seed <- options$ga_seed

  # load needed libraries
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("MCMCpack"))
  suppressPackageStartupMessages(library("GA"))
  source("scripts/train_library.R")
  
  # load training data
  train <- load_training_data(traindatafile,0)
  train_ga <- get_ga_train_data(train)
  trainy <- train_ga$trainy
  trainx <- train_ga$trainx[,-1:-5]
  
  # load decoy data
  y <- NULL
  X <- NULL
  INFO <- NULL
  
  # load initial  
	initial_model <- "data/initial_parameters_tmp.txt"
	initial_model <- read.table(initial_model, col.names = c("nucleus","atompair_id","para"))
	rownames(initial_model) <- initial_model$atompair_id
      
  # refine preliminary bayesian model maximizing the ability to resolve native from non-native conformations in a decoy pool
  initial_model <- runGA(trainy, trainx, initial_model, y=y, X=X, INFO=INFO, cycles = ga_cycles, population_size = ga_populations, std = ga_std, rmsd_threshold = ga_rmsd_threshold, seed = ga_seed, verbose = options$verbose)
  
  # output model to file
  write.table(initial_model, file = output_parameters, col.names = TRUE, row.names = FALSE, quote = FALSE)
  if(options$verbose)
    print(initial_model)
}
