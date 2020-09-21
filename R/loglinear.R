#' Calculates file-level risk measures using a loglinear model.
#' 
#' @param data Data frame containing the data to be evaluated.
#' @param weight Column name for sampling weights.
#' @param varpool Vector of column names to be used in model.
#' @param degree Highest degree of interaction terms to be used in the model.
#' @param numiter Maximum number of iterations to run iterative proportional 
#'   fitting for the loglinear model.
#' @param epsilon Maximum deviation allowed between observed and fitted margins.
#' @param blanks_as_missing If TRUE, character and factor variables that are
#'   blank or pure whitespace are treated as missing values.
#' @param output_filename Name of the csv file to save the data set with 
#'   record-level risk measures, .tau1_rec and .tau2_rec, attached. NULL if no 
#'   output file is to be saved.
#' @details The data should not contain any missing values among \code{varpool}
#'   variables or the \code{weight} variable.
#' @return An object of type \code{sdc_loglinear} containing calculated risk
#' measures.
#' @examples
#' data(exampledata)
#' vars <- c("BORNUSA", "CENREG", "DAGE3", "DRACE3", "EDUC3", "GENDER")
#' wgt <- "WEIGHT"
#'
#' results <- sdc_loglinear(exampledata, wgt, vars, degree=3)
#' print(results)
#' plot(results, plotvar1="BORNUSA", plotvar2="WEIGHT")
#' @export
sdc_loglinear <- function(data, 
                          weight,
                          varpool,
                          degree=2,
                          numiter=40,
                          epsilon=0.001,
                          blanks_as_missing=TRUE,
                          output_filename=NULL) {

  # Input validation -----------------------------------------------------------

  sdc_loglinear_input_checks(data, 
                             weight,
                             varpool,
                             degree,
                             numiter,
                             epsilon,
                             blanks_as_missing,
                             output_filename)

  pi_results <- loglinear_workhorse(data, 
                                    weight, 
                                    varpool,
                                    degree,
                                    numiter,
                                    epsilon,
                                    fixed_pi=TRUE)

  pi_k_results <- loglinear_workhorse(data, 
                                      weight, 
                                      varpool,
                                      degree,
                                      numiter,
                                      epsilon,
                                      fixed_pi=FALSE)

  tables <- list(pi_results=pi_results$statistics,
                 pi_k_results=pi_k_results$statistics)

  return_value <- list(ll_tables=tables,
                       data_with_statistics=pi_k_results$data_with_statistics,
                       options=list(weight=weight, 
                                    varpool=varpool, 
                                    degree=degree,
                                    numiter=numiter,
                                    epsilon=epsilon,
                                    blanks_as_missing=blanks_as_missing,
                                    output_filename=output_filename))


  if(!is.null(output_filename)) {

    write.csv(pi_k_results$data_with_statistics, 
              output_filename,
              row.names=FALSE)

  } # end if

  class(return_value) <- "sdc_loglinear"
  return(return_value)

} # end 'loglinear_nway'

# ------------------------------------------------------------------------------

# This does most of the work within sdc_loglinear
loglinear_workhorse <- function(data, 
                                weight,
                                varpool,
                                degree=2,
                                numiter=40,
                                epsilon=0.01,
                                fixed_pi=TRUE) {

  # Get tables -----------------------------------------------------------------

  freq_tbl <- get_table(data, varpool, weight)

  tmp <- setNames(lapply(varpool, 
                         function(variable) unique((data[[variable]]))),
                  varpool)

  grid <- merge(expand.grid(tmp), freq_tbl, all.x=TRUE)
  grid[is.na(grid)] <- 0 

  # Check number of sample uniques ---------------------------------------------

  num_unique <- sum(grid$.unweighted_freq == 1)

  if(num_unique == 1)
    warning(paste("Only 1 sample unique. Some statistics will not",
                  "be available and appear as NA."))

  # Get average cell size ------------------------------------------------------

  avg_cell_size <- get_average_cell_size(data, varpool)
  
  # If no sample uniques, stop. ------------------------------------------------

  if(num_unique == 0) {

    values <- as.matrix(data.frame(sampsize=nrow(data),
                                   avg_cell_size=avg_cell_size,
                                   tau1Risk=0,
                                   tau2Risk=0,
                                   tau1=0,
                                   tau2=0,
                                   samp_unique=FALSE))

    row.names(values) <- NULL
    return(values)

  } # end if

  # Fit model ------------------------------------------------------------------

  varterms <- paste(varpool, collapse=" + ")
  formula <- as.formula(paste(".weighted_freq ~", varterms))

  if(degree == 1)
    formula2 <- as.formula(paste("~", varterms))
  else
    formula2 <- as.formula(paste0("~ (", varterms, ")^", degree))

  tab <- xtabs(formula, grid)

  model <- MASS::loglm(formula2, tab, fitted=TRUE, iter=numiter, eps=epsilon)

  # Get lambdas ----------------------------------------------------------------

  lambda_tbl <- as.data.frame.table(fitted(model))

  names(lambda_tbl)[names(lambda_tbl) == "Freq"] <- ".lambda"

  grid <- merge(grid, lambda_tbl, all.x=TRUE)
  lambda <- grid$.lambda

  # Get number of samples per cell, find sample uniques ------------------------

  nsamp <- grid$.unweighted_freq
  unique_sample <- nsamp == 1

  # Calculate statistics -------------------------------------------------------
  
  nsamp <- grid$.unweighted_freq
  unique_sample <- nsamp == 1

  pi <- rep(nrow(data)/sum(data[[weight]]), length(unique_sample))

  if(!fixed_pi) {

    pi_k <- grid$.unweighted_freq/grid$.weighted_freq
    pi <- ifelse(nsamp == 0, pi, pi_k)

  }

  wgt <- ifelse(grid$.unweighted_freq == 0, 
                0,
                grid$.weighted_freq/grid$.unweighted_freq)
  p <- ifelse(nsamp >= 1, nsamp/wgt, 0)

  arg1 <- ifelse(unique_sample & wgt <= 1, wgt,
          ifelse(unique_sample, p/(1 - p)*log(1/p), 0))
  arg2 <- ifelse(unique_sample, p, 0)

  exp1 <- ifelse(unique_sample, 
                 (1/((1 - pi)*lambda))*(1 - exp(-(1 - pi)*lambda)),
                 0)
  pusu <- ifelse(unique_sample, exp(-(1 - pi)*lambda), 0)

  r <- ifelse(lambda != 0,
              (1/((1 - pi)*(lambda)))*(1 - exp(-(1 - pi)*lambda)),
              0)

  a1 <- ifelse(lambda != 0, (1 - pi)*(lambda)*exp(-lambda), 0)
  b1 <- ifelse(pi != 0, ((1 - pi)^2)*(1/(2*pi))*(lambda)*exp(-lambda), 0)

  aa <- ifelse(lambda != 0, 
               a1*(nsamp - pi*lambda) + b1*(((nsamp - pi*lambda)^2) - nsamp),
               0)

  vaa <- a1*a1*pi*lambda + 2*b1*b1*pi*lambda*pi*lambda

  a2 <- ifelse(lambda != 0, exp(-pi*lambda)*r - (exp(-lambda)), 0)
  b2 <- ifelse(lambda != 0 & pi != 0, 
               ((1/(pi*lambda))*exp(-pi*lambda)*r) - 
                 ((1/(pi*lambda))*exp(-lambda)* (1 + .5*((1 - pi)*lambda))),
               0)
  bb2 <- a2*(nsamp - pi*lambda) + b2*(((nsamp - pi*lambda)^2) - nsamp)
  vbb <-a2*a2*pi*lambda+2*b2*b2*pi*lambda*pi*lambda

  stderr <- function(x) return(sd(x)/sqrt(length(x)))

  aanew <- mean(aa)/stderr(aa) #This is B1/sqrt(Vr)
  bbnew <- mean(bb2)/stderr(bb2) #This is B2/sqrt(Vr)

  scaa1 <- sum(aa)/sqrt(sum(vaa)) #This is B1/sqrt(V)
  scbb2 <- sum(bb2)/sqrt(sum(vbb)) #This is B2/sqrt(V)

  B_tau1_type1 <- mean(aa)/stderr(aa) #This is B1/sqrt(Vr)
  B_tau1_type2 <- sum(aa)/sqrt(sum(vaa)) #This is B1/sqrt(V)
  B_tau2_type1 <- mean(bb2)/stderr(bb2) #This is B2/sqrt(Vr)
  B_tau2_type2 <- sum(bb2)/sqrt(sum(vbb)) #This is B2/sqrt(V)

  tau2 <- sum(exp1)
  tau1 <- sum(pusu)
  tau1argus <- sum(arg1)
  tau2argus <- sum(arg2)

  values <- as.matrix(data.frame(sampsize=nrow(data),
                                 avg_cell_size=avg_cell_size,
                                 tau1Risk=tau1/sum(nsamp),
                                 tau2Risk=tau2/sum(nsamp),
                                 tau1=tau1,
                                 tau2=tau2,
                                 B_tau1_type1=B_tau1_type1,
                                 B_tau1_type2=B_tau1_type2,
                                 B_tau2_type1=B_tau2_type1,
                                 B_tau2_type2=B_tau2_type2,
                                 samp_unique=TRUE))

  if(degree > 1) {

    values <- rbind(loglinear_workhorse(data, 
                                        weight,
                                        varpool,
                                        1,
                                        numiter,
                                        epsilon,
                                        fixed_pi)$statistics,
                    values)

  } # end if

  row.names(values) <- NULL

  # Add risk measures to data set ----------------------------------------------

  if(!fixed_pi) { # Only need to calculate once

    keep_names <- !names(grid) %in% 
                  c(".weighted_freq", ".unweighted_freq", ".lambda")
    tmp <- grid[,keep_names,drop=FALSE]
    tmp$.tau1_rec <- pusu
    tmp$.tau2_rec <- exp1

    data <- merge(data, tmp, all.x=TRUE)

    if(any(is.na(data$.tau1_rec)))
      data$.tau1_rec[is.na(data$.tau1_rec)] <- 0

    if(any(is.na(data$.exp1)))
      data$.tau2_rec[is.na(data$.tau2_rec)] <- 0

  } # end if

  # Return statistics and modified data set ------------------------------------

  return(list(statistics=values, data_with_statistics=data))

} # end 'loglinear_workhorse'

# ------------------------------------------------------------------------------

#' @describeIn sdc_loglinear S3 print method for sdc_loglinear objects
#' 
#' Prints tables of file-level reidentification risk measures.
#' 
#' @param x Object of class sdc_loglinear, as returned by sdc_loglinear.
#' @param summary_outfile Name of summary output .txt file. If not NULL, console
#'   output is copied to the file. Default is NULL (no logging of output).
#'   Errors and warnings are not diverted (consider running in batch mode if
#'   logging is needed).
#' @param ... Currently unused. For NextMethod compatibility.
#' @export
print.sdc_loglinear <- function(x,
                                summary_outfile=NULL,
                                ...) {

  # Log output if requested ----------------------------------------------------

  if(!is.null(summary_outfile)) {

    stopifnot(is.character(summary_outfile))

    conn <- file(summary_outfile)
    sink(conn, append=FALSE, split=TRUE)

    on.exit(sink())
    on.exit(close(conn), add=TRUE)

  } # end if

  degree <- x$options$degree

  # Fixed pi -------------------------------------------------------------------

  cat("RESULTS - Uses overall average weight\n")
  loglinear_print_workhorse(x$ll_tables$pi_results, degree)

  # pi and pi_k ----------------------------------------------------------------
  
  cat("\nRESULTS - Uses cell average weights\n")
  loglinear_print_workhorse(x$ll_tables$pi_k_results, degree)

  return(invisible(NULL))

} # end 'print.sdc_loglinear'

# ------------------------------------------------------------------------------

loglinear_print_workhorse <- function(ll_table, degree) {

  if(degree == 1)
    interaction <- data.frame(interaction="none")
  else
    interaction <- data.frame(interaction=c("none", paste0(degree,"-way")))

  if(any(!ll_table[,"samp_unique"])) {

    interaction$interaction <- "No Sample Uniques"
    tbl1 <- cbind(as.data.frame(ll_table)[,c("sampsize", "avg_cell_size")],
                  interaction,
                  as.data.frame(ll_table)[,c("tau1Risk", "tau2Risk")])
    tbl2 <- as.data.frame(ll_table)[,5:6]

  } else {

    tbl1 <- cbind(as.data.frame(ll_table)[,c("sampsize", "avg_cell_size")],
                  interaction,
                  as.data.frame(ll_table)[,c("tau1Risk", "tau2Risk")])
    tbl2 <- as.data.frame(ll_table)[,5:10]

  } # end if...else

  old_scipen <- options("scipen")$scipen
  on.exit(options(scipen=old_scipen))
  options(scipen=9)
  print(tbl1, row.names=FALSE, digits=5)
  
  cat("\n")
  print(tbl2, row.names=FALSE, digits=5)

  return(invisible(NULL))

} # end 'loglinear_print_workhorse'

# ------------------------------------------------------------------------------

#' @describeIn sdc_loglinear S3 plot method for \code{sdc_loglinear} objects
#' 
#' Produces boxplots and scatterplots of record-level risk measures, tau1 and
#' tau2.
#' 
#' @param plotpath Directory to save plots. Plots are saved as \emph{jpeg} files 
#'   (quality = 100\%). If the directory does not exist, it is first created.
#'   If \code{plotpath} is NULL (default), plots are not saved.
#' @param plotvar1 A vector of names of discrete variables for boxplots. If 
#'   none, boxplots are not produced.
#' @param plotvar2 A vector of names of continuous variables for scatterplots. 
#'   If none, scatterplots are not produced.
#' @export
plot.sdc_loglinear <- function(x,
                               plotpath=NULL,
                               plotvar1=character(0),
                               plotvar2=character(0),
                               ...) {

  dat <- x$data_with_statistics

  # Input validation -----------------------------------------------------------

  stopifnot(is.null(plotpath) || 
            (is.character(plotpath) && length(plotpath) == 1))
  stopifnot(is.character(plotvar1) && all(plotvar1 %in% names(dat)))
  stopifnot(is.character(plotvar2) && all(plotvar2 %in% names(dat)))

  # Set save_plots flag --------------------------------------------------------

  save_plots <- !is.null(plotpath)

  # Run this again without saving, so that plots also appear in R console.
  if(save_plots)
    plot(x, plotvar1=plotvar1, plotvar2=plotvar2)

  # Helper functions -----------------------------------------------------------

  get_plot_filename <- function(plot_type, is_tau1_rec, variable) {

    fname <- paste0(plot_type,
                    "--", 
                    ifelse(is_tau1_rec, "tau1--", "tau2--"),
                    variable,
                    ".jpeg")

    return(file.path(plotpath, fname))

  } # end 'get_plot_filename'

  # ----------------------------------------------------------------------------

  # Call before making any plot
  plot_start <- function(plot_type, is_tau1_rec, variable) {

    if(save_plots)
      jpeg(get_plot_filename(plot_type, is_tau1_rec, variable), quality=100)
    else
      dev.new()

    return(invisible(NULL))

  } # end 'plot_start'

  # ----------------------------------------------------------------------------
  
  # Call after making any plot
  plot_end <- function() {

    if(save_plots)
      dev.off()

    return(invisible(NULL))

  } # end 'plot_end'

  # ----------------------------------------------------------------------------
  
  make_tau1_boxplot <- function(plotvar) {

    plot_start("box", is_tau1_rec=TRUE, plotvar)
    boxplot_formula <- as.formula(paste(".tau1_rec ~", plotvar))
    boxplot(boxplot_formula, 
            data=dat,
            ylab="tau1")
    plot_end()

    return(invisible(NULL))

  } # end 'make_tau1_boxplot'

  # ----------------------------------------------------------------------------
  
  make_tau2_boxplot <- function(plotvar) {

    plot_start("box", is_tau1_rec=FALSE, plotvar)
    boxplot_formula <- as.formula(paste(".tau2_rec ~", plotvar))
    boxplot(boxplot_formula, 
            data=dat,
            ylab="tau2")
    plot_end()

    return(invisible(NULL))

  } # end 'make_tau2_boxplot'

  # ----------------------------------------------------------------------------
  
  make_tau1_scatterplot <- function(plotvar) {

    plot_start("scatter", is_tau1_rec=TRUE, plotvar)
    scatrplt <- 
      ggplot2::ggplot(dat, 
                      mapping=ggplot2::aes_string(x=plotvar, 
                                                  y=".tau1_rec")) + 
      ggplot2::geom_point() +
      ggplot2::ylab("tau1")
    plot(scatrplt)
    plot_end()

    return(invisible(NULL))

  } # end 'make_tau1_scatterplot'

  # ----------------------------------------------------------------------------
  
  make_tau2_scatterplot <- function(plotvar) {

    plot_start("scatter", is_tau1_rec=FALSE, plotvar)
    scatrplt <- 
      ggplot2::ggplot(dat, 
                      mapping=ggplot2::aes_string(x=plotvar, 
                                                  y=".tau2_rec")) + 
      ggplot2::geom_point() +
      ggplot2::ylab("tau2")
    plot(scatrplt)
    plot_end()

    return(invisible(NULL))

  } # end 'make_tau2_scatterplot'

  # ----------------------------------------------------------------------------
  
  # Create output directory (if any) if it does not already exist --------------

  # Clean up directory name
  if(!is.null(plotpath))
    plotpath <- file.path(plotpath)

  if(!is.null(plotpath) && !dir.exists(plotpath))
    dir.create(plotpath, recursive=TRUE)

  # Produce tau1 boxplots ------------------------------------------------------

  lapply(plotvar1, make_tau1_boxplot)

  # Produce tau1 scatterplots --------------------------------------------------

  lapply(plotvar2, make_tau1_scatterplot)

  # Produce tau2 boxplots ------------------------------------------------------

  lapply(plotvar1, make_tau2_boxplot)

  # Produce tau2 scatterplots --------------------------------------------------

  lapply(plotvar2, make_tau2_scatterplot)

  return(invisible(NULL))

} # end 'plot.sdc_loglinear'

# ------------------------------------------------------------------------------

sdc_loglinear_input_checks <- function(data,
                                       weight,
                                       varpool,
                                       degree,
                                       numiter,
                                       epsilon,
                                       blanks_as_missing,
                                       output_filename) {

  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)

  stopifnot(is.character(weight))
  stopifnot(length(weight) == 1)
  stopifnot(weight %in% names(data))

  stopifnot(is.character(varpool))
  stopifnot(length(varpool) >= 1)
  stopifnot(all(varpool %in% names(data)))

  stopifnot(is.numeric(degree))
  stopifnot(length(degree) == 1)
  stopifnot(degree >= 1)
  stopifnot(degree == round(degree)) # is integer?
  stopifnot(degree <= length(varpool))

  stopifnot(is.numeric(numiter))
  stopifnot(length(numiter) == 1)
  stopifnot(numiter >= 1)
  stopifnot(numiter == round(numiter)) # is integer?

  stopifnot(is.numeric(epsilon))
  stopifnot(length(epsilon) == 1)
  stopifnot(epsilon > 0)

  stopifnot(is.logical(blanks_as_missing))

  stopifnot(is.null(output_filename) || 
            (is.character(output_filename) && 
             length(output_filename) == 1))

  # Check for missing values ---------------------------------------------------

  missing_variables <- sdc_loglinear_check_missing(data, 
                                                   varpool,
                                                   weight,
                                                   blanks_as_missing)

  # Throw exception if any missing values found --------------------------------
  
  if(length(missing_variables) > 0) {

    msg <- paste0("sdc_loglinear does not support missing values. ",
                  "Missing values found among: ", 
                  paste0(missing_variables, collapse=" "))

    stop(msg)

  }

  return(invisible(NULL))

} # end 'sdc_loglinear_input_checks'

# ------------------------------------------------------------------------------

sdc_loglinear_check_missing <- function(data, 
                                        varpool,
                                        weight,
                                        blanks_as_missing) {

  has_missing <- function(x)
    return(any(is.na(x)) || 
           (blanks_as_missing && any(trimws(x) == "")))

  variables <- c(varpool, weight)
  missing <- variables[sapply(c(varpool, weight),
                              function(var)
                                return(has_missing(data[,var])))]

  return(missing)

} # end 'sdc_loglinear_check_missing'

# ------------------------------------------------------------------------------

get_average_cell_size <- function(data, varpool) {

  num_of_cells <- 
    prod(sapply(varpool, 
                function(var) 
                  length(unique(as.vector(na.omit(data[[var]]))))
               )
        )

  return(nrow(data)/num_of_cells)

} # end 'get_average_cell_size'