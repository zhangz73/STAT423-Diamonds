library(corrplot)
library(data.table)
library(sets)
library(dequer)
library(MASS)

diamond <- read.csv("diamonds.csv")

get_compressed_df <- function(dataframe, predictor_to_compress, 
                              compressed_vector, bin.width){
  df_new <- c()
  for(v in compressed_vector){
    tmp <- dataframe[dataframe[[predictor_to_compress]] <= v + bin.width & 
                       dataframe[[predictor_to_compress]] >= v - bin.width,]
    tmp[[predictor_to_compress]] <- v
    df_new <- rbind(data.frame(df_new), data.frame(tmp))
  }
  return(df_new)
}

compressed_df <- get_compressed_df(diamond, "carat", c(0.3, 0.5, 0.8, 1, 1.5), 0.03)

non_greedy_model_selection <- function(all_vector, response, init_vector = NULL,
                                       criteria = c("AIC", "BIC", "LOOCV"),
                                       n.best = 3, greedy = T, 
                                       diamond = compressed_df){
  if(!criteria %in% c("AIC", "BIC", "LOOCV")){
    stop("criteria can only be AIC, BIC, LOOCV")
  }
  if(is.null(init_vector)){
    init_vector <- c()
  }
  get_model <- function(vec){
    if(length(vec) == 0){
      return(lm(as.formula(paste0(response, "~ 1")), data = diamond))
    }
    predictors <- paste(vec, collapse = "+")
    return(lm(as.formula(paste0(response, "~", predictors)), data = diamond))
  }
  compute_score <- function(lm){
    if(criteria == "AIC"){
      return(AIC(lm))
    } else if(criteria == "BIC"){
      return(BIC(lm))
    } else{
      return(mean((residuals(lm)/(1-hatvalues(lm)))^2))
    }
  }
  
  predictors_queue <- queue()
  score_queue <- queue()
  pushback(predictors_queue, init_vector)
  init.model <- get_model(init_vector)
  init.score <- compute_score(init.model)
  pushback(score_queue, init.score)
  best_models <- list(init_vector)
  for(i in 2:n.best){
    best_models[[i]] <- init_vector
  }
  best_scores <- rep(init.score, n.best)
  existing_vectors <- set(as.set(init_vector))
  j <- 0
  total <- 2^(length(all_vector) - length(init_vector))
  while(length(predictors_queue) > 0){
    j <- j + 1
    if(j %% 1 == 0){
      cat("\r", paste0("We have finished ", j, " trials out of up to ", 
                       total, " trials!"))
    }
    curr_vector <- pop(predictors_queue)
    curr_score <- pop(score_queue)
    additional_predictors <- all_vector[!all_vector %in% curr_vector]
    
    tmp_models <- NULL
    tmp_scores <- NULL
    i <- 0
    for(pred in additional_predictors){
      if(!set_contains_element(existing_vectors, as.set(c(curr_vector, pred)))){
        i <- i + 1
        model <- get_model(c(curr_vector, pred))
        score <- compute_score(model)
        existing_vectors <- set_union(existing_vectors, 
                                      set(as.set(c(curr_vector, pred))))
        if(score < curr_score){
          if(!greedy){
            pushback(predictors_queue, c(curr_vector, pred))
            pushback(score_queue, score)
          } else{
            tmp_models[[i]] <- c(curr_vector, pred)
            tmp_scores[i] <- score
          }
          if(score < max(best_scores)){
            idx <- which.max(best_scores)
            best_scores[idx] <- score
            best_models[[idx]] <- c(curr_vector, pred)
          }
        }
      }
    }
    if(greedy && !is.null(tmp_scores)){
      ord <- order(tmp_scores)
      for(i in 1:min(length(tmp_scores), n.best)){
        pushback(predictors_queue, tmp_models[[i]])
        pushback(score_queue, tmp_scores[i])
      }
    }
  }
  ordering <- order(best_scores)
  ret_scores <- c()
  ret_models <- c()
  for(o in ordering){
    ret_scores <- c(ret_scores, best_scores[o])
    ret_models <- c(ret_models, as.formula(paste0(response, "~", 
                                                  paste(best_models[[o]], 
                                                        collapse = "+"))))
    #ret_models <- c(ret_models, summary(get_model(best_models[[o]])))
  }
  return(list(best_models = ret_models, best_scores = ret_scores))
}

init_vector <- c("carat", "cut", "color", "clarity", "depth", "table", "x", 
                 "carat*cut", "carat*color", "carat*clarity", "cut*color",
                "cut*clarity", "color*clarity", "cut*depth", "carat*depth",
                "depth*table", "depth*x", "table*x", "carat*table", "carat*x")
additional <- c("cut*table", "cut*x", "color*depth", "color*table", "color*x",
                "clarity*depth", "clarity*table", "clarity*x")

methods <- c("AIC", "BIC", "LOOCV")
for(m in methods){
  print(paste0("According to ", m, ":"))
  non_greedy_model_selection(all_vector = c(init_vector, additional), response = "price", 
                             criteria = m, n.best = 3, greedy = F,
                             init_vector = init_vector)
}