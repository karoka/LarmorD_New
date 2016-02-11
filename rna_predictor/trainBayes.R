suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("MCMCpack"))
suppressPackageStartupMessages(library("GA"))
source("train_library.R")

train <- load_training_data("train_features_H1p.txt","H1'",0)
train_ga <- get_ga_train_data(train)
trainy <- train_ga$trainy
trainx <- train_ga$trainx

decoy_data <- load_decoy_data("larmord_shifts.txt","H1'",1)
y <- decoy_data$y
X <- decoy_data$X
INFO <- decoy_data$INFO

# train preliminary model
paras <- runBayesian(train)


compute_mae_score(paras$para,trainy,trainx)
# initialize genetic algorithm optimization
refined_paras <- runGA(trainy,y,trainx,X,paras,INFO)