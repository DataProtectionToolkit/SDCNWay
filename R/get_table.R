get_table <- function(x, vars, wgt=NULL) {
# Get (possibly weighted) cross tabulation for x[,vars].

  if(is.null(wgt))
    tmp <- x[,vars,drop=FALSE]
  else
    tmp <- cbind(x[,vars,drop=FALSE], wgt=x[, wgt])

  tmp <- na.omit(tmp)

  if(nrow(tmp) == 0)
    return("empty")
  
  results <- expand.grid(lapply(tmp[1:length(vars)], 
                                function(x) unique(x)),
                                stringsAsFactors = FALSE)

  if(!is.null(wgt)) {

    xtab_formula <- as.formula(eval(paste0("wgt~", paste(vars, collapse ="+"))))
    results <- merge(results, as.data.frame(xtabs(xtab_formula, tmp)), by=vars)

  }

  results <- merge(results, 
                   as.data.frame(xtabs(data=tmp[,1:length(vars),drop=FALSE])), 
                   by=vars)

  if(is.null(wgt))
    results <- dplyr::rename(results, .unweighted_freq=Freq)
  else
    results <- dplyr::rename(results, 
                             .weighted_freq=Freq.x,
                             .unweighted_freq=Freq.y)

  results <- results[results$.unweighted_freq > 0,]

  if(nrow(results) == 0)
    return("empty")

  return(results)

} # end 'get_table'