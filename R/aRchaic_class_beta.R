

aRchaic_class_beta = function(mat,
                              type = c("mutation", "mutation-flank", "mutation-pos","mutation-flank-pos",
                                  "specific-mutation-pos", "wo-strand", "wo-strand-break"),
                              class_labs = NULL,
                              class_method = c("SVM", "classtpx"),
                              normalize = TRUE,
                              svm.control = list(),
                              classtpx.control = list()){
  
  classtpx.control.default <- list(method="theta.fix", shrink=FALSE,
                                   shrink.method = 1, tol=0.001,
                                   ord=FALSE)
  svm.control.default <- list(scale = TRUE, type = NULL, kernel ="radial", 
                              degree = 3, 
                              coef0 = 0, cost = 1, nu = 0.5,
                              class.weights = NULL, cachesize = 40, tolerance = 0.001, epsilon = 0.1,
                              shrinking = TRUE, cross = 0, fitted = TRUE)
  
  
  classtpx.control <- modifyList(classtpx.control.default, classtpx.control)
  svm.control <- modifyList(svm.control.default, svm.control)
  
  
  if(type == "mutation"){
    mat_reduced <- filter_by_mutation(mat)
  }else if (type == "mutation-pos"){
    mat_reduced <- filter_by_mutation_pos(mat, max_pos, flanking_bases)
  }else if (type == "mutation-flank"){
    mat_reduced <- filter_by_mutation_flank(mat)
  }else if (type == "mutation-flank-pos"){
    mat_reduced <- filter_by_mutation_flank_pos(mat, max_pos, flanking_bases)
  } else if (type == "specific-mutation-pos"){
    if(is.null(pattern)){
      stop("for this option, user needs to provide a mutation type: C->T, C->A C->G, T->A, T->C or T->G")
    }
    mat_reduced <- filter_by_pos_pattern(mat, max_pos, pattern = pattern, flanking_bases)
  }else if (type == "wo-strand"){
    mat_reduced <- filter_out_strand(mat)
  }else if (type == "wo-strand-break"){
    mat_reduced <- filter_out_strand_break(mat)
  }else {
    stop("the type of filtering does not match with the possible options: see documentation")
  }
  
  rownames(mat_reduced) <- rownames(mat)
  
  pooled_data <- mat_reduced
  
  if(length(class_labs) != dim(pooled_data)[1]){
    stop("number of samples not same as number of class labels, aborting")
  } 
  
  test_indices <- which(is.na(class_labs))
  if(length(test_indices)==0){
    stop("no test sample provided, aborting")
  }
  
  train_indices <- which(!is.na(class_labs))
  if(length(train_indices)==0){
    stop("no training sample provided, aborting")
  }
  
  
  if (class_method == "SVM"){
    if(normalize){
      voom_pooled_data <- t(limma::voom(t(pooled_data))$E);
      trainX <- voom_pooled_data[train_indices,]
      testX <- voom_pooled_data[test_indices,]
      y = factor(class_labs[train_indices])
      data <- cbind.data.frame(y, trainX);
      library(e1071)
      model_SVM <- do.call(e1071::svm, append(list(formula = y ~ ., 
                                                   data=data, 
                                                   probability=TRUE), svm.control))
      prob  = predict(model_SVM, testX, probability=TRUE)
    }else{
      trainX <- pooled_data[train_indices,]
      testX <- pooled_data[test_indices,]
      y = factor(class_labs[train_indices])
      data <- cbind.data.frame(y, trainX);
      library(e1071)
      model_SVM <- do.call(e1071::svm, append(list(formula = y ~ ., 
                                                   data=data, 
                                                   probability=TRUE), svm.control))
      
      prob  = predict(model_SVM, testX, probability=TRUE)
    }
    return(list("model" = model_SVM, "test_class_prob" = prob))
  }
  
  if(class_method == "classtpx"){
    known_samples <- train_indices
    class_labs_train <- class_labs[train_indices]
    model_classtpx <- do.call(classtpx::class_topics, append(list(counts = as.matrix(pooled_data),
                                                                  K = length(unique(class_labs_train)),
                                                                  known_samples = known_samples,
                                                                  class_labs = class_labs_train), classtpx.control))
    
    omega_test <- model_classtpx$omega[test_indices,]
    return(list("model" = model_classtpx, "test_class_prob" = omega_test))
  }
  
}