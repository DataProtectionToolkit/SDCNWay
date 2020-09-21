calculate_risk_measures <- function(data, wgt, varpool, tau1, tau2) {

  # Helper function ------------------------------------------------------------

  # Data should have mu-argus and related statistics attached,
  # i.e., should be first processed by 'add_mu_argus_and_statistics'
  el_emam_risk_measures <- function(data, tau1, tau2)
    return(list(pRa=mean(data$.cell_count < 1/tau1),
                pRb=1/min(data$.cell_count),
                pRc=length(unique(data$.cell_ID))/nrow(data),
                jRa=mean(data$.mu_argus > tau2),
                jRb=max(data$.mu_argus),
                jRc=mean(data$.mu_argus)))

  # Missing wgt not allowed. Check wgt. ----------------------------------------

  stopifnot(!missing(wgt) && 
            !is.null(wgt) && 
            inherits(wgt, "character") &&
            length(wgt) == 1 &&
            wgt %in% names(data))

  # Add mu-argus and other statistics ------------------------------------------

  data <- add_mu_argus_and_statistics(data, wgt, varpool)

  # Add el-emam risk measures and return ---------------------------------------

  return(list(data=data,
              mu_argus_summary=summarize_mu_argus(data, wgt),
              el_emam_measures=el_emam_risk_measures(data, tau1, tau2)))

} # end 'calculate_risk_measures'

# ------------------------------------------------------------------------------

calculate_risk_measures_unweighted <- function(data, varpool, tau1) {
# Same as calculate_risk_measures but excludes statistics that require a
# weighting variable.

  # Helper functions -----------------------------------------------------------

  # Data should have .cell_ID and .cell_count attached,
  # i.e., should be first processed by 'add_cell_ID', then 'add_cell_count'
  el_emam_risk_measures <- function(data)
    return(list(pRa=mean(data$.cell_count < 1/tau1),
                pRb=1/min(data$.cell_count),
                pRc=length(unique(data$.cell_ID))/nrow(data)))

  # ----------------------------------------------------------------------------
 
  add_cell_count <- function(data_with_cell_ID) {

    stopifnot(nrow(data_with_cell_ID) > 0)

    ids <- data_with_cell_ID$.cell_ID
    tmp <- data_with_cell_ID[,".cell_ID",drop=FALSE]
    tbl <- aggregate(tmp, by=list(.cell_ID=ids), FUN=length)
    names(tbl) <- c(".cell_ID", ".cell_count")

    return(merge(data, tbl))
    
  } # end 'add_statistics_per_cell'

  # Add cell IDs and cell counts -----------------------------------------------

  data <- add_cell_count(add_cell_ID(data, varpool))

  return(list(data=data,
              el_emam_measures=el_emam_risk_measures(data)))

} # end 'calculate_risk_measures_unweighted'

# ------------------------------------------------------------------------------

add_mu_argus_and_statistics <- function(data, wgt, varpool) {

  # Helper functions -----------------------------------------------------------

  # Term used in mu-argus score when f_k > 3
  get_mu_argus_term <- Vectorize(function(n, f_k, pi_k){
    return(factorial(n)*((1-pi_k)^n)/prod(f_k + 1:n))
  }, vectorize.args="n")

  # ----------------------------------------------------------------------------

  get_score <- function(f_k, pi_k) {
  # Return approximated mu-argus score for one cell given f_k and pi_k

    if(isTRUE(all.equal(pi_k, 1))) # pi_k approximately equal to 1
      return(1/f_k)

    if(f_k == 1)
      return(-log(pi_k)*pi_k/(1 - pi_k))
    else if(f_k == 2)
      return((pi_k*log(pi_k) + 1 - pi_k)*pi_k/((1 - pi_k)^2))
    else if(f_k == 3)
      return(((1 - pi_k)*(3*(1 - pi_k) - 2) - 2*pi_k*pi_k*log(pi_k))*pi_k/
             (2*((1-pi_k)^3)))
    else
      return((1 + sum(get_mu_argus_term(1:7, f_k, pi_k)))*pi_k/f_k)

  } # end 'get_score'

  # Add cell IDS and other needed statistics -----------------------------------

  data <- add_cell_ID(data, varpool)
  data <- add_statistics_per_cell(data, wgt)

  data$.mu_argus <- mapply(get_score, data$.cell_count, data$.pi_k)

  return(data)

}  # end 'add_mu_argus_and_statistics'

# ------------------------------------------------------------------------------

summarize_mu_argus <- function(data, wgt) {
# Create a table summarizing mu-argus scores by cell count

  max_counts <- c(1:3, Inf)

  n <- nrow(data)

  summarize_subset <- function(max_count) {

    dat <- data[data$.cell_count <= max_count,]

    sum_mu_argus <- sum(dat$.mu_argus)
    return(data.frame(num_cases=nrow(dat),
                      total_mu_argus=sum_mu_argus,
                      mean_mu_argus=mean(dat$.mu_argus),
                      total_mu_over_total_obs=sum_mu_argus/n,
                      total_mu_over_sum_wts=sum_mu_argus/sum(data[,wgt])))

  } # end 'summarize_subset'

  tbl <- do.call(rbind, lapply(max_counts, summarize_subset))
  tbl <- cbind(data.frame(max_count=max_counts), tbl)

  # Remove rows if no new information ------------------------------------------
  
  tbl <- tbl[1:(which(tbl$num_cases == tbl$num_cases[4])[1]),]

  return(tbl)

} # end 'summarize_mu_argus'

# ------------------------------------------------------------------------------

add_statistics_per_cell <- function(data_with_cell_ID, wgt) {
# Needed to calculate risk measures 
# Estimate and add pi_k using cell counts and weights. Also, add sum of weights.

  ids <- data_with_cell_ID$.cell_ID
  tmpdat <- data_with_cell_ID[,wgt,drop=FALSE]

  tbl <- aggregate(tmpdat, by=list(.cell_ID=ids), FUN=sum)
  names(tbl) <- c(".cell_ID", ".sum_weights")

  data <- merge(data_with_cell_ID, tbl)

  tbl <- aggregate(tmpdat, by=list(.cell_ID=ids), FUN=length)
  names(tbl) <- c(".cell_ID", ".cell_count")

  data <- merge(data, tbl)

  data$.pi_k <- data$.cell_count/data$.sum_weights

  return(data)

} # end 'add_statistics_per_cell'


