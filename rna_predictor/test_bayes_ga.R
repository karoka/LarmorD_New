#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
  make_option(c("--comparison_fparameters"), type="character",default=NULL,
              help="additional parameter set [default %default]"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="print header and progress information [default %default]")
)

parser <- OptionParser(usage = "%prog [options] refined_bayesian_model cs_decoy_database", option_list=option_list)

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
  refined_bayesian_model <- arguments$args[1]
  decoydatafile <- arguments$args[2]
  comparison_fparameters <- options$comparison_fparameters

    # load needed libraries
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("MCMCpack"))
  suppressPackageStartupMessages(library("GA"))
  source("train_library.R")
  
  # load parameters
  refined_bayesian_model <- read.table(refined_bayesian_model, header = T)
  
  # load additional parameter set
  if (!is.null(comparison_fparameters)){
    comparison_parameters <- "data/initial_parameters_tmp.txt"
    comparison_parameters <- read.table(comparison_fparameters, col.names = c("nucleus","atompair_id","para"))
    rownames(comparison_parameters) <- comparison_parameters$atompair_id
    comparison_parameters <- comparison_parameters[xnames[-1:-5], "para"]
    refined_bayesian_model$larmord <- c(rep(0,5), comparison_parameters)
  }
    
  # load decoy data
  decoy_data <- load_decoy_data(decoydatafile, 0)
  y <- decoy_data$y
  X <- decoy_data$X
  INFO <- decoy_data$INFO
  
  # test refined_bayesian_model
  print(compute_slr_score(refined_bayesian_model$ga, y, X, INFO, rmsd_threshold = 2.5, testing = TRUE, key = "Bayesian-GA"))
  print(compute_slr_score(refined_bayesian_model$para, y, X, INFO, rmsd_threshold = 2.5, testing = TRUE, key = "Bayesian"))
  if(!is.null(comparison_fparameters))
    print(compute_slr_score(refined_bayesian_model$larmord, y, X, INFO, rmsd_threshold = 2.5, testing = TRUE, key = "Larmord"))
}
