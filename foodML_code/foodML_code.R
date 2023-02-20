#========================================================================#
### Script for Predicting Food Intake w/ Random Forest ###

# This analyses is conducted under the comments of 1st submission to
# Journal of Nuitrition.

# Add cross-validated results for classification

# Rerun analysis with LASSO 
# Compare it with pure LASSO, no need for variable selection

# Create dummy data for code reproducibility

# Multi-class ROC

#========================================================================#
# Setup the required R packages (aka libraries) to run the entire script #
library(easypackages)
packages("stringr", "qdapRegex", "glmnet", "vegan", "car", "randomForest", "readxl", "caret",
         "pROC", "ROCR", "changepoint")

############################################################################
######################### Section 0: Data Cleaning #########################

#' The goal of this section is to read and perform an initial cleaning of
#' the raw Silva dataset before computing the difference between end and baseline
#' which were the features we used in the study.

# Set the working directory as the directory with the dataset (**Modify this your computer**)
setwd("./")

# Read raw data (This is the data that is uploaded in Box)
dat = read_xlsx("../Dataset/foodML_data.xlsx")   # Read the excel dataset

dat = dat[,c(2:3, 5:ncol(dat))]    # Retain the relevant information

# Convert all demographic variables to factors
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)
dat$SubjectID = as.factor(dat$SubjectID)
dat$Period = as.factor(dat$Period)

# Replace period information that are missing (NA for not available) as period 1
dat[which(dat$Period == '_'),4] = 1

# Lowercase baseline and end
dat$BaselineEnd = tolower(dat$BaselineEnd)

# Remove missing avocado information
dat = dat[duplicated(dat$SubjectID) | duplicated(dat$SubjectID, fromLast = T),]

# Reorder the dataset according to 1) subject ID, 2) treatment, and 3) baseline/end
dat = dat[order(dat$SubjectID, dat$Period, dat$BaselineEnd),]

# View the current dataset in an excel format to check that it makes sense
# Only view the first 20 columns to save memory
View(dat[,1:20])

# Broccoli is not sorted properly from the above step due to the mislabeling of the period
# Manually sort broccoli
datBroc = dat[str_detect(dat$Food, "Broccoli"),]
datBroc = datBroc[order(datBroc$SubjectID, datBroc$Food, datBroc$BaselineEnd),]
dat = rbind(dat[!str_detect(dat$Food, "Broccoli"),], datBroc)

# Remove almonds with missing period information 
alm = dat[which(dat$Study == "Almond"),]
rownames(alm) = 1:nrow(alm)
alm = alm[-c(151, 154),]   # Subject 4517 has no baseline for period 2 and no end for period 4 

# Remove butter 
alm = alm[which(alm$Treatment != 'Butter'),]
dat = rbind(dat[dat$Study != "Almond",], alm)

# Difference of baseline and end Microbiomes across all subjects
datEnd = dat[which(1:nrow(dat) %% 2 == 0),]    # Even rows are end
datBase = dat[which(1:nrow(dat) %% 2 != 0),]   # Odd rows are baseline
diffDat = cbind(datBase[,1:7], datEnd[,8:ncol(datEnd)] - datBase[,8:ncol(datBase)])

### Create the dataset from 4/5/6 diets and Silva features ###
### (Only working with the difference in abundance starting from here) ###

# Extract the Silva features
SV = diffDat[,c(1:7, which(grepl("SV", colnames(diffDat))))]

# Retain only species level features 
levels = str_count(colnames(SV), "__")    # Count the number of double underlines for all features 
spe.ind = which(levels == 7)   # Species level features have 7 double underlines
SV = SV[,c(1:7, spe.ind)]

# Slightly modify the microbes name to a format that supports random forest
colnames(SV) = gsub("SV_F_16S_","",colnames(SV))
colnames(SV) = gsub(";",".", colnames(SV))
colnames(SV) = make.names(colnames(SV))

### Create the datasets for different number of diets ###

# Function to create the 6 diets dataset
# Grain is split into barley and oats #
split6diets = function(dset) {
  
  grain = factor(dset[which(dset$Food == 'Grains'), 7])
  dset$Food = as.character(dset$Food)
  dset[which(dset$Food == 'Grains'), 2] = as.character(grain)
  dset$Food = factor(dset$Food)
  
  return(dset)
}

#<<No function is needed for the 5 diets>>#

# Function to create the 4 diets dataset 
# Broccoli is removed
split4diets = function(dset) {
  dset.4 = dset[which(dset$Study != 'Broccoli'),]
  dset.4$Food = factor(dset.4$Food)
  return(dset.4)
}

#**Create the 4/5/6 diets dataset**#
SV.6 = split6diets(SV)  
SV.5 = SV
SV.4 = split4diets(SV)

# View table to verify
View(SV.6)
levels(SV.6$Food)    # Check that Barley and Oat are present as a factor level

###################################################################################
###################### Section 1: Classification (part 1) #########################

#' This is the first part of the Classification step in the manuscript.
#' 1. Marginal Screening w/ Krusal-Wallis test for each study separately 
#' 2. Pool the top 10 microbes from each study into a larger set of features, this is our feature matrix X.
#' 3. Fit a Random Forest model on the response Y (the foods) and the pooled feature matrix X 

# The individual functions for each task is presented first.
# A final wrapper function at the end of this section (named 'pipe.line') that combines everything together.

#------------------------------------------------------------------------#
#-------------- Marginal Screening using Kruskal Wallis -----------------#

# This is a non-parameteric method for hypothesis testing, testing whether the deviance of the distribution between two samples
# is significant (significantly different if pvalue < 0.05). Non-parametric here means that we do not assume an underlying
# distribution for the samples. Parametric methods, such as the well-known ANOVA, is the parametric analogy of the 
# Kruskal-Wallis test, where we assume both distributions to be approximately normal. The reason we are doing an 
# non-paramteric version is due to the compositional nature of the data, which does not follow a normal distribution. 

#> dset: input dataset of a SINGLE FOOD STUDY
#> top: how many features we want to retain for each study
ks.screen = function(dset, top) {
  
  # Reorganize the data, so that we have the binary treatment indicator (whether the subject belongs to the treatment group)
  # as the response and the difference in abundance of the microbes as features
  y = ifelse(grepl("No", dset$Food), 0, 1)
  temp.dat = cbind(y, dset[,-c(1:7)])
  
  # Perform the Kruskal-Wallis test for Marginal Screening
  res = c()
  for (i in 2:ncol(temp.dat)) { res = c(res, kruskal.test(temp.dat[,i] ~ y)$p.value) }
  
  # Combind the results 
  feat.res = data.frame(cbind(names(temp.dat[,-1]), res))
  feat.res[,2] = as.numeric(as.character(feat.res[,2]))
  colnames(feat.res) = c("Microbiome", "p_value")
  
  # Sort by p-value 
  feat.res = feat.res[order(feat.res$p_value),]
  
  # Retain top features
  feat.res = feat.res[1:top,]
  
  # Round to 4 significant values
  feat.res[,2] = round(feat.res[,2], 4)
  
  # Output a table with the microbes and their p-value from the KW-test
  return(feat.res)
}

# Screen the microbes across different studies
# Then pool them into one long feature table
# This will be the set of features to serve as input in the random forest model

#> dset: input dataset of ALL SUBJECTs
#> top: how many features we want to retain for each study
screen.all = function(dset, top) {
  
  diets = as.character(unique(dset$Study))
  
  # Extract top significant features from each study
  feat.list = list()
  for (i in 1:length(diets)) {
    temp.set = dset[which(dset$Study == diets[i]),]
    feat.list[[i]] = ks.screen(temp.set, top)
    names(feat.list)[i] = diets[i]
  }
  
  # Pool the top significant features from each diet together
  pool.feats = c()
  for (i in 1:length(feat.list)) {
    pool.feats = c(pool.feats, as.character(feat.list[[i]]$Microbiome))    
  }
  
  # Remove the duplicate features due to pooling
  pool.feats = unique(pool.feats)   
  
  # Output an array of the names of features to be passed into the random forest model
  return (pool.feats)
}

#------------------------------------------------------------------------#
#----------------------- Fit Random Forest Model ------------------------#

# Run a single random forest

#> dset: A cleaned dataset, with the response Y as the different studies and X as the pooled features from marginal screening.
#> nsize: Control 'nodesize', a tuning parameter of random forest. 
rf.fit = function(dset, nsize) {
  
  # Fit the random forest model 
  fit.rf = randomForest(factor(y)~., data = dset, ntree = 500, importance = TRUE, nodesize = nsize)
  
  # Extract the out-of-bag confusion matrix from the trained random forest model 
  oobmat = fit.rf$confusion
  oobmat[,ncol(oobmat)] = 1 - oobmat[,ncol(oobmat)]
  colnames(oobmat)[ncol(oobmat)] = "classAccu"
  
  # Calculate the weighted classification error
  oob.err = fit.rf$confusion[,ncol(fit.rf$confusion)]
  weight.err = sum(prop.table(table(factor(dset$y))) * oob.err)
  
  # Calculate the multi-class AUC
  auc = as.numeric(multiclass.roc(factor(dset$y), fit.rf$votes)$auc)
  
  # Organize the variable importance table
  # Sort the importance table in descending order, with the most important features at the top
  imp.tab = c()
  if (ncol(dset) > 2) {
    for (i in 1:length(unique(dset$y))) {
      curr = sort(fit.rf$importance[,i], decreasing = T)
      curr.tab = cbind(names(curr), round(as.numeric(curr), 5))
      colnames(curr.tab) = c(as.character(unique(dset$y)[i]), "Importance")
      imp.tab = data.frame(cbind(imp.tab, curr.tab))
    }
    all.imp = fit.rf$importance[,length(unique(dset$y))+1]
    all.imp = data.frame(sort(all.imp, decreasing = T))
  } else {
    all.imp = c()
  }
  accu = c(1 - weight.err, auc)
  names(accu) = c("Accuracy", "AUC")
  
  # Output a list of objects, including the confusion matrix, weighted classification accuracy
  # variable importance table for each food, and the overall variable importance table
  return(list(fit = fit.rf, mat = oobmat, accu = accu, 
              imp = imp.tab, all.imp = all.imp))   
}

# Test a bunch of nodesizes for random forest 
# The best result across difference nodesizes will be selected in the pipeline function below.
rf.runode = function(dset) {
  
  # Consider nodesize 1-10
  nsize = c(1:10); res = list(); sum.accu = 0
  for (i in 1:length(nsize)) {
    res[[i]] = rf.fit(dset, nsize[i])
    sum.accu = sum.accu + res[[i]]$accu
  }
  res$mean.accu = sum.accu / length(nsize)
  
  return(res)  
}

#------------------------------------------------------------------------#
#----------------------- KW TEST + RANDOM FOREST ------------------------#

# This is a wrapper function that combines all above functions into one

#> dset: An input difference dataset (the SV.4 / SV.5 / SV.6 dataset we created in the first section)
pipe.line = function(dset, top, seed) {
  
  # Set seed for reproducibility 
  set.seed(seed)
  
  # Obtain pooled feature list from marginal screening 
  pipe.diet = screen.all(dset, top)
  
  # Create dataset to be input into random forest model
  # Retain only the treatment subjects
  all.pipe = cbind(dset$Food, dset[,pipe.diet]); colnames(all.pipe)[1] = "y"
  pipe.treat = all.pipe[which(!grepl("No", all.pipe$y)),]
  
  # Run random forest model for different nodesizes
  all.pipe.rf = rf.runode(pipe.treat)
  
  # Return the mean confusion matrix and error
  # Along with the variable importance 
  accu.vec = all.pipe.rf[[1]]$accu
  for (i in 2:10) { accu.vec = c(accu.vec, all.pipe.rf[[i]]$accu) }
  
  # Select the best result out of all nodesizes
  acc.vec = accu.vec[names(accu.vec) == "Accuracy"]
  best = all.pipe.rf[[which(acc.vec == max(acc.vec))[1]]]
  num = nrow(best$all.imp)
  
  # Returns a list of results, including:
  # 1) the number of features that were inputs of the random forest model
  # 2) confusion matrix from best model
  # 3) classification error from best model
  # 4) variable importance table for each food from best model
  # 5) overall variable importance table from best model 
  return(list(num = num, mat = best$mat, accu = best$accu, imp = best$imp, all.imp = best$all.imp))
}

# Generate the results #
sv6.spe.20 = pipe.line(SV.6, 20, 1031)
sv5.spe.20 = pipe.line(SV.5, 20, 1031)
sv4.spe.20 = pipe.line(SV.4, 20, 1031)

sv5.spe.20$num
sv5.spe.20$mat
sv5.spe.20$accu
sv5.spe.20$imp
sv5.spe.20$all.imp

sv6.spe.20$accu
sv5.spe.20$accu
sv4.spe.20$accu

###################################################################################
###################### Section 2: Classification (part 2) #########################

#' This is the second part of the classification step. Here we take the 
#' top microbes ranked by the variable importance from part 1 and refit
#' a random forest model. The goal here is to identify a smaller set of 
#' microbes that still has an adequate accuracy for classification.

#---------------------------------------------------------------------------------------------------#
#> dset: input dataset
#> rf.res: A random forest result list from part 1
#> num: number of top microbes to retain from the rf.res (default to 10)
#> mode: 1) Take the top important microbes from the overall variable importance table as the features 
#>       2) Pool the top important microbes from each food as the features for the model (This is what we did for our paper)

refit.rf = function(dset, rf.res, num = 10, mode = 2) {
  
  set.seed(1031)
  
  # Extract treatment
  dset = dset[which(!grepl("No", dset$Food)),]
  
  # Construct a new, reduced features dataset according to the specified mode.
  # Mode 2
  dlist = rf.res$imp; feats = c()
  for (i in 1:length(unique(dset$Study))) { feats = c(feats, as.character(dlist[1:num,2*i-1])) }
  feats.uq = unique(feats)
  dat = cbind(factor(dset[,2]), dset[,colnames(dset) %in% feats.uq])

  colnames(dat)[1] = "y"
  
  if (mode == 1) {
    # Refit random forest with updated feature set across different nodesizes
    nsize = c(1:10); res = list()
    for (i in 1:length(nsize)) { res[[i]] = rf.fit(dat, nsize[i]) }
    
    # Return the mean confusion matrix and error
    # Along with the variable importance 
    accu = res[[1]]$accu
    for (i in 2:10) { accu = c(accu, res[[i]]$accu) }
    
    # Select the best result 
    min.ind = which(accu == min(accu))[1]
    min.res = res[[min.ind]]
    min.res$num = length(feats)
    
    return(min.res)    
  } else {
    # Refit a single random forest with nodesize = 2
    resFit = rf.fit(dat, 2)
    return(resFit)
  }
}

#---------------------------------------------------------------------------------------------------#
# This function performs the hypothesis tests to choose the appropriate number of features
# per diet. Both the McNemar and K-fold repeated cross-validation t-test are conducted, but only
# the latter is reported in the manuscript due to the results making more sense.
sigTest = function(mcRes, paramRes) {
  
  pRes = rbind(colMeans(paramRes), apply(paramRes, 2, function(x) sd(x)), -1, -1, -1, -1, -1, -1)
  mcRes[which(mcRes == T)] = 1; mcRes[which(mcRes == F)] = 0
  
  # Perform the hypothesis testing to test the significance 
  for (i in 2:(ncol(paramRes))){
    
    # Sequential Difference
    pRes[3, i] = pRes[1, i] - pRes[1, i-1]
    
    # Kfold Test Sequential
    k1 = t.test(paramRes[,i-1], paramRes[,i], paired = TRUE, alternative = "two.sided")
    pRes[4,i] = k1$p.value
      
    # McNemar's Test Sequential
    m1 = mcnemar.test(mcRes[,i-1], mcRes[,i])
    pRes[5,i] = m1$p.value
    
    # Pair difference w/ 20
    pRes[6,i] = pRes[1, i] - pRes[1, 1]
    
    # Kfold Pair w/ 20
    k2 = t.test(paramRes[,1], paramRes[,i], paired = TRUE, alternative = "two.sided")
    pRes[7,i] = k2$p.value
    
    # McNemar's Test w/ 20
    m2 = mcnemar.test(mcRes[,1], mcRes[,i])
    pRes[8,i] = m2$p.value
  }
  
  rownames(pRes) = c("mean", "se", "seq change", "seq k-fold", "seq mcnemar",
                     "comp20 change", "comp20 k-fold", "comp20 mcnemar")    

  return(round(pRes, 4))
}

#-----------------------------------------------------------------------------------------#
# This generates Figure 3 in our accepted manuscript.

param = c(20, 15, 12, 10, 7, 5, 3, 2) # Define the parameters of interest

refit.allimp = function(dset, param, rf.res, nfolds, string = " ", plot = F) {
  
  set.seed(1031)
  
  # Extract treatment
  dTreat = dset[which(!grepl("No", dset$Food)),]
  dlist = rf.res$imp 
  
  # Use K-fold cross-validated paired t-test
  # Construct 75%-25% split, rerun nfolds times
  cvSplit = createDataPartition(factor(dTreat$Food), times= nfolds, p = 0.75)
  
  # Gather the cross-validation information for classification
  paramRes = c()
  mcRes = c()
  for (p in 1:length(param)) {
    
    print(param[p])
    
    # Construct a new, reduced features dataset according to the specified mode.
    feats = c()
    for (i in 1:length(unique(dset$Study))) { feats = c(feats, as.character(dlist[1:param[p],2*i-1])) }
    feats.uq = unique(feats)
    dat = data.frame(cbind(factor(dTreat[,2]), dTreat[,colnames(dTreat) %in% feats.uq]))
    colnames(dat)[1] = "y"
    
    # Store results for mcnemars test
    mcFit = rf.fit(dat, 2)
    mcRes = cbind(mcRes, mcFit$fit$predicted == dat$y)
    
    # Create array to store the results
    acc = c()
    for (i in 1:length(cvSplit)) {
      trainDat = dat[cvSplit[[i]],]; testDat = dat[-cvSplit[[i]],]
      resFit = rf.fit(trainDat, 2)
      rfPred = predict(resFit$fit, testDat[,-1])
      acc = c(acc, mean(rfPred == testDat$y))
    }
    paramRes = cbind(paramRes, acc)
  }
  colnames(paramRes) = colnames(mcRes) = param
  
  # significance test
  pRes = sigTest(mcRes, paramRes)
  
  if (plot) {
    boxplot(paramRes)
  }

  return(pRes)
}


# Output the accuracy vs. number of important features plot for each dataset
sv.6.spe2 = refit.allimp(SV.6, param, sv6.spe.20, 50) 
sv.5.spe2 = refit.allimp(SV.5, param, sv5.spe.20, 50, T)  # Set the last parameter to T to plot Figure 3
sv.4.spe2 = refit.allimp(SV.4, param, sv4.spe.20, 50) 

# From the plot and our meeting discussions that 10 microbes from each diet is a compact enough
# feature set while also retaining a high classification accuracy.
# Generate the results in the manucript #
sv.6.rf = refit.rf(SV.6, sv6.spe.20, 10, 2)
sv.5.rf = refit.rf(SV.5, sv5.spe.20, 10, 2) 
sv.4.rf = refit.rf(SV.4, sv4.spe.20, 10, 2)

# View the results
sv.5.rf$num
sv.5.rf$mat
sv.5.rf$accu
sv.5.rf$imp
sv.5.rf$all.imp

######################################################################################
### Section 2.1 (Revision) ###

# This part is where we compared the random forest model with the LASSO model (linear) via 
# cross-validation as requested in the review. This part is unimportant and can be skipped to PART 3.

# Classification Models 
fitMod = function(trainTreat, testTreat) {

  # Run LASSO
  laFit = cv.glmnet(data.matrix(trainTreat[,-1]), factor(trainTreat$y), alpha = 1, nfolds = 5,
                    family = "multinomial")
  laRes = calcRes(laFit, testTreat, 1)  
    
  # Run random forest
  rfFit = randomForest(trainTreat[,-1], trainTreat[,1], ntree = 500, nodesize = 3)
  rfRes = calcRes(rfFit, testTreat, 2)
  
  return(list(rf = rfRes, la = laRes))
}

# Calculate confusion matrix
fitObj = rfFit
calcRes = function(fitObj, testTreat, mode) {
  
  # Calculate the AUC
  if (mode == 1) {
    prob = predict(fitObj, data.matrix(testTreat)[,-1], type = "response")
    auc = multiclass.roc(testTreat$y, prob[,,1])
  } else if (mode == 2) {
    prob = predict(fitObj, data.matrix(testTreat)[,-1], type = "prob")
    auc = multiclass.roc(testTreat$y, prob)
  }

  # Calculate the other information
  pred = as.factor(predict(fitObj, data.matrix(testTreat)[,-1], type = "class"))
  levels(pred) = levels(testTreat$y)
  confMat = confusionMatrix(pred, testTreat$y)
  
  # Confusion matrix 
  conMat = t(confMat$table)
  accuVec = rep(0, ncol(conMat))
  for (i in 1:nrow(conMat)) { accuVec[i] = conMat[i,i] / rowSums(conMat)[i] }
  conMat = cbind(conMat, accuVec) 
  colnames(conMat)[ncol(conMat)] = "accu"
  
  # Accuracy
  res = round(c(confMat$overall[1], auc$auc) * 100, 2)
  names(res) = c("Accuracy", "AUC")
  
  return(list(res = res, mat = conMat))
}

# Run classification for LASSO and Random Forest with Cross-Validation 
dset = SV.5
rfObj = sv5.spe.20
nfolds= 10
runClass = function(dset, rfObj, nfolds) {
  
  set.seed(1031)
  
  # Create training set 
  # Extract treatment
  dTreat = dset[which(!grepl("No", dset$Food)),]
  dlist = rfObj$imp 
  
  # Construct a new, reduced features dataset according to the specified mode.
  dat = data.frame(cbind(factor(dTreat[,2]), dTreat[,colnames(dTreat) %in% dlist[,1]]))
  colnames(dat)[1] = "y"
  
  # Use K-fold cross-validated paired t-test
  # Construct 75%-25% split, rerun nfolds times
  cvSplit = createDataPartition(factor(dat$y), times= nfolds, p = 0.75)
  
  # Create array to store the results
  resTab = c()
  resMat = list()
  for (i in 1:length(cvSplit)) {
    trainDat = dat[cvSplit[[i]],]; testDat = dat[-cvSplit[[i]],]
    modelFits = fitMod(trainDat, testDat)
    resTab = rbind(resTab, c(modelFits$rf$res, modelFits$la$res))
    resMat[[i]] = list(rf = modelFits$rf$mat, la = modelFits$la$mat)
  }
  
  resMat[[length(resMat) + 1]] = resTab
  
  return(resMat)
}

sv6.rev.1 = runClass(SV.6, sv6.spe.20, 50)
sv5.rev.1 = runClass(SV.5, sv5.spe.20, 50)
sv4.rev.1 = runClass(SV.4, sv4.spe.20, 50)

resObj = sv6.rev.1
meanSD = function(resObj) {
  
  # Classification results
  tabl = resObj[[length(resObj)]]
  mean = colMeans(tabl)
  sd = apply(tabl, 2, function(x) sd(x))
  res = rbind(mean, sd)
  colnames(res) = paste0(c("RF ", "RF ", "LASSO ", "LASSO "), colnames(res))
  
  # Confusion matrix results
  temp = resObj[-length(resObj)]
  rfmat = temp[[1]]$rf
  lamat = temp[[1]]$la
  for (i in 2:(length(temp))) {
    rfmat = rfmat + temp[[i]]$rf
    lamat = lamat + temp[[i]]$la
  }
  rfmat = rfmat / length(temp)
  lamat = lamat / length(temp)
  
  return(list(accu = res, rfmat = round(rfmat, 2), lamat = round(lamat, 2)))
}

meanSD(sv6.rev.1)
meanSD(sv5.rev.1)
meanSD(sv4.rev.1)

######################################################################################
############################ Section 3: Validation ###################################

#' This part deals with the validation of the results obtained from Section 2 to make sure 
#' that we are predicting the food and not the study. The is done by removing the main control
#' effect from the treatment subjects using. Mathematically, we subtract the top N PCs 
#' (determined by the elbow plot) of the treatment features from the control features. 
#' Then a new random forest model is fitted on the TREATMENT EFFECT FEATURES. We also fit this model
#' on the control values.
#' 
#' The goal is to show that the 1) classification accuracy does not decrease by a lot after removing the control effect
#'                              2) The treatment effect model on the control values does significantly worse than the accuracy
#'                                 of the treatment effect model by itself 

#--------------------------------------------------------------------------------------------------------------------#
# This function shows the elbow plot of the variance explained by each principal components of the control features.
# The goal of this function is to identify how many PCs of the control group ("the main effect") is to be removed. 
# Visually, we want to see where the % of variance explained flattens out on the elbow plot.

#> dset: Input dataset
#> feats: A random forest result object from Section 2
#> top: Number of principal components of the control group to be subtracted from the treatment group
pca.plot = function(dset, rf.res, num) {
  
  # Extract number of studies considered (4/5/6)
  studies = length(unique(study))
  
  # Extract important variables from each diet
  dlist = rf.res$imp
  feats = c()
  for (i in 1:length(unique(dset$Study))) {
    feats = c(feats, as.character(dlist[1:num,2*i-1]))
  }
  feats.uq = unique(feats)
  dset = cbind(factor(dset[,2]), dset[,colnames(dset) %in% feats.uq])
  colnames(dset)[1] = "y"
  
  # Calculate SVD on control group
  cont = dset[which(grepl("No", dset$y)),]
  cont.svd = svd(cont[,-1])  
  eigen = cont.svd$d^2
  
  b = plot(eigen[1:10], main = paste0("Elbow Plot (", studies, " Foods)"), ylab = "% of Variance Explained", 
           xlab = "Principal Component", type = "o", pch = 4)
  return(round(eigen[1:10], 4))
}

# Output the PCA plots #
pca.plot(SV.6, sv.6.rf, 10)   # Title says 5 foods but it should be 6 because barley and oats are regarded as 1 here
# % of variance explained flattens out after the 4th PC, thus we removed the top 4 
pca.plot(SV.5, sv.5.rf, 10)   # Same interpretation as the 6 foods case
pca.plot(SV.4, sv.4.rf, 10)   # % of variance explained flattens out after the 5th PC, thus we removed the top 5

#-----------------------------------------------------------------------------------------------------#
# This function computes the PCs from the control group and removed the top N (found from previous function)
# PCs from the treatment group, to generate the features that contain the treatment effects.

#> dset: Input dataset
#> feats: A random forest result object from Section 2
#> top: Number of principal components of the control group to be subtracted from the treatment group
pca.ex = function(dset, feats, top) {
  
  # Take out the treatment group
  treat = dset[which(!grepl("No", dset$Food)), feats]   
  temp.y = dset[which(!grepl("No", dset$Food)), 2]   # Response variable for treatment group
  
  # SVD of the control group
  cont = dset[which(grepl('No', dset$Food)),feats]
  cont.svd = svd(cont)
  
  # Project the treatment group onto the direction of the control group
  treat.p = as.matrix(treat) %*% cont.svd$v[,1:top] %*% t(cont.svd$v[,1:top])
  treat.p[which(abs(treat.p) < 1e-12)] = 0
  
  # Subtract the control effect from treatment 
  new.treat = cbind(temp.y, treat - treat.p)  
  colnames(new.treat)[1] = "y"
  
  return(new.treat)  
}

#---------------------------------------------------------------------------------------------------#
# Wrapper function for fitting RF on the treatment effect features.
# Input the control values into the treatment effect model.
# Show that the accuracy from the treatment effect model is much better than the accuracy from the treatment 
# effect model fitted on the control values, indicating that our classification result is valid and the treatment
# effect is driving the prediction of the model and not the control effect.

#> dset: Input dataset
#> feats: A random forest result object from Section 2
#> top: Number of principal components of the control group to be subtracted from the treatment group
fitObj = fit
testTreat = cont

resCalc = function(fitObj, testTreat) {
  
  # Transform the name of the control values
  contY = factor(sub("No", "", testTreat[,1]))
  
  # Calculate the auc
  prob = predict(fitObj, data.matrix(testTreat)[,-1], type = "prob")
  pred = as.character(predict(fitObj, data.matrix(testTreat)[,-1], type = "class"))
  if (ncol(prob) == 6) {
    prob[,3] = (prob[,3] + prob[,5]) / 2
    prob = prob[,-5]
    colnames(prob)[3] = "Grains"
    pred[which(pred == "Oat" | pred == "Barley")] = "Grains"
    pred = factor(pred)
  }
  auc = multiclass.roc(contY, prob)
  
  # Calculate the other information
  #levels(pred) = levels(factor(contY))
  confMat = confusionMatrix(pred, contY)
  
  # Confusion matrix 
  conMat = t(confMat$table)
  accuVec = rep(0, ncol(conMat))
  for (i in 1:nrow(conMat)) { accuVec[i] = conMat[i,i] / rowSums(conMat)[i] }
  conMat = cbind(conMat, accuVec) 
  colnames(conMat)[ncol(conMat)] = "accu"
  
  # Accuracy
  res = round(c(confMat$overall[1], auc$auc) * 100, 2)
  names(res) = c("Accuracy", "AUC")
  
  return(list(res = res, mat = conMat))
}

dset = SV.6
rf.res = sv.6.rf
comp.test = function(dset, rf.res, top) {
  
  set.seed(1031)
  
  # Set up the data set, including only the features considered in Section 2
  features = rownames(rf.res$all.imp)
  new.dat = data.frame(cbind(dset[,2], dset[,features]))
  colnames(new.dat)[1] = "y"
  
  # Remove the control effects from treatment 
  treat.pca = pca.ex(dset, features, top)  
  
  # Run random forest on the features with control effect removed
  pca.list = list(); pca.accu = c()
  for (i in 1:10) {
    pca.list[[i]] = rf.fit(treat.pca, i); 
    pca.accu = rbind(pca.accu, pca.list[[i]]$accu)
  }
  
  #---------------------------------------------------------------------------------#
  # Run treatment model on control values and hopefully sees it performs much worse #
  cont = new.dat[which(grepl("No", new.dat$y)),]
  fit = randomForest(factor(y)~., data = treat.pca,
                     ntree = 1000, importance = TRUE, nodesize = which.max(pca.accu[,1]))
  
  # and a list for results for the treatment effect model fitted on the control values
  return(list(pca = pca.list[[which.max(pca.accu[,1])]], cont = resCalc(fit, cont)))
}

# Output the results
sv6.pca = comp.test(SV.6, sv.6.rf, 4)
sv5.pca = comp.test(SV.5, sv.5.rf, 4)
sv4.pca = comp.test(SV.4, sv.4.rf, 4)

# View the results
sv4.pca$pca$mat; sv4.pca$pca$accu        # Confusion matrix of the treatment effect model
sv4.pca$cont;      # Confusion matrix of the treatment effect model fitted on control values




















