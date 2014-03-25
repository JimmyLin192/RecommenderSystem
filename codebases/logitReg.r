# AUTHOR: JIMMY LIN

library("glmnet")

training_data <- read.csv("./sampleData/joint.csv", head=FALSE, sep=",")
NUM_OF_RECORDS = nrow(training_data)
NUM_OF_FEATURES = ncol(training_data)

X = as.matrix(training_data[,-1])
Y = training_data[,nCols]
cat ("nRecords: ", NUM_OF_RECORDS, "\n")
cat ("nFeatures: ", NUM_OF_FEATURES, "\n")

## PRE-DEFINE THE PARAMETERS FOR THE SETTING
NFOLDS = 10 
ALPHA_PREC = 1e-3;
MODEL = "binomial" # for logistic regression


cat("Basis function: ", MODEL, "\n")
cat("Alpha_precision: ", ALPHA_PREC, "\n")
cat("nFolds: ", NFOLDS, "\n")

# CONFIGURATION BASED ON GIVEN SETTINGS
WEIGHTS = rep(1, nrow(X))  
# REGULARIZATION MIXTURE COEFFICIENT
ALPHAS = seq(0, 1, ALPHA_PREC) 
NALPHAS = length(ALPHAS)
# PRE-COMPUTE THE FOLD-ID (INDICATE WHICH VECTOR IT IS IN)
FOLD_ID = sample (rep(seq(NFOLDS), length = NUM_OF_RECORDS))

cat(nrow(X))
cat("\n")
cat(length(Y))
cat("\n")
Y = matrix(Y, nrow= length(Y), ncol = 1)

fitting = cv.glmnet(X, Y, family=MODEL, foldid=FOLD_ID,
                    WEIGHTS, alpha=0.5, type.measure="deviance")
