#' Calculate risk measures through exhaustive tabulations, Mu-Argus, and other 
#'   methods.
#'
#' This function primarily uses the exhaustive tabulation method to quantify 
#' disclosure risk. It tabulates cell counts for different combinations of 
#' variables provided by the user. Using these counts, this function identifies
#' variable categories and records which are considered high risk for
#' disclosure. File-level re-identification risk measures are also provided, 
#' e.g., Mu-Argus (Polettini 2003) and the risk metrics promosed in El Emam 
#' (2011).
#' 
#' @param data Data frame containing the data for which we are to measure 
#'   disclosure risk. Unexpected behavior may result if any column name begins
#'   with a period.
#' @param ID Name of column which identifies records. If NULL (default), an
#'   ID column named .ROW_NUMBER is created and used in reports.
#' @param weight Column name for sampling weights. NULL or empty if none.
#' @param varpool Vector of column names over which to form tables.
#' @param forcelist Vector of variable names. Some are included in all 
#'   tabulations. Optional. 
#' @param forcenum Number of variables in \code{forcelist} that are mandatory for 
#'   all tabulations. That is, all tabulations will have a number of variables
#'   from forcelist exactly equal to \code{forcenum}.
#' @param missingdef A named list specifying missing values. The names
#'   correspond to column names in \code{data}. 
#' @param mindim Integer specifying the minimum number of \code{varpool} variables 
#'   (including \code{forcelist} variables) that can be used to form tables. 
#' @param maxdim Integer specifying the maximum number of \code{varpool} variables
#' (including ]code{forcelist} variables) that can be used to form tables. 
#' @param threshold Threshold to determine the number of violations in terms of 
#'   cell counts. If the number of cases in a cell is less than \code{threshold}, 
#'   the cell is flagged as a violation. If threshold is NULL and wgthreshold is
#'   not NULL, then only a weighted threshold will be used. If both are NULL,
#'   threshold will be set to 3 and the weighted threshold will not be used.
#' @param wgtthreshold Threshold to determine violations in terms of weighted
#'   cell counts. If NULL, a weighted threshold will not be used.
#' @param condition Character string describing how weighted and unweighted
#'   thresholds are combined when both are used. If used, it must be "and" or 
#'   "or" (case insensitive). This parameter is ignored if \code{weight} is NULL.
#' @param output_filename Name of the csv file to save the data set with 
#'   violation counts and Mu-Argus scores attached. NULL if no output file is to
#'   be saved.
#' @param tau1 A threshold to compute the risk measure, pRa. See User Manual for
#'   more details.
#' @param tau2 A threshold to compute the risk measure, jRa. This parameter is 
#' ignored if \code{weight} is NULL. See User Manual for more details.
#' @param include_mu_argus Flag indicating whether Mu-Argus and El-Emam 
#'   metrics should be calculated.
#' @return An object of type \code{sdc_extabs}. Internally, a named list of 
#'   statistics.
#' \describe{
#'   \item{\emph{tabulation}}{Cell counts and violation flags. Represented as
#'      a list with each element corresponding to a \code{varpool} combination.}
#'   \item{\emph{data_with_statistics}}{The original data with new 
#'         columns showing statistics such as violation counts and Mu-Argus 
#'         score for each record.}
#'   \item{\emph{recoded_data_with_statistics}}{Same as 
#'      \code{data_with_statistics} but with missing value recodes.}
#'   \item{\emph{mu_argus_summary}}{Summary table of Mu-Argus by cell count. For
#'      this summary, all variables in \code{varpool} are used to define a cell. If
#'      weight is NULL, then this summary is omitted.}
#'   \item{\emph{el_emam_measures}}{List of file-level re-identification risk 
#'     measures.}
#'   \item{\emph{percent_violations_by_var_and_level}}{Table with percent of
#'      records that are in violation for each variable/category.}
#'   \item{\emph{percent_violations_by_dim_var_and_level}}{Table with percent 
#'      of cells that are in violation for each dimension/variable/category.}
#'   \item{\emph{options}}{Options provided to \code{sdc_extabs} by the user, 
#'      such as \code{missingdef}, \code{mindim}, etc.}
#' } 
#' @details If a specified missing value contains
#'   only whitespace, it will match any element with only whitespace. NA values
#'   in data are treated as missing regardless of \code{missingdef}. If you do not
#'   want NA values to be treated as missing, please recode them before
#'   passing the data to this function.
#' 
#'   Note that if a weight variable is not provided, the number of statistics
#'   and plots that are produced is significantly reduced.
#' @examples 
#' data(exampledata)
#' vars <- c("BIB1201", "BIC0501", "BID0101", "BIE0601", "BORNUSA", "CENREG",
#'           "DAGE3", "DRACE3", "EDUC3", "GENDER")
#' results <- sdc_extabs(exampledata, 
#'                       ID="CASEID",
#'                       weight="WEIGHT", 
#'                       varpool=vars,
#'                       mindim=2,
#'                       maxdim=3,
#'                       missingdef=list(BIE0601=5),
#'                       wgtthreshold=3000,
#'                       condition="or")
#' print(results, cutoff=15)
#' plot(results, plotvar1="BORNUSA", plotvar2="WEIGHT")
#' @references 
#' \insertRef{el2011methods}{SDCNway}
#' 
#' \insertRef{polettini2003some}{SDCNway}
#' @importFrom graphics boxplot plot
#' @importFrom methods as
#' @importFrom stats aggregate as.formula fitted na.omit sd setNames xtabs
#' @importFrom utils combn write.csv
#' @importFrom grDevices dev.new dev.off jpeg
#' @importFrom Rdpack reprompt
#' @export
sdc_extabs <- function(data, 
                       ID=NULL,
                       weight=NULL,
                       varpool=names(data),
                       forcelist=character(0),
                       forcenum=1,
                       missingdef=list(),
                       mindim=1,
                       maxdim=2,
                       threshold=NULL,
                       wgtthreshold=NULL,
                       condition=NULL,
                       output_filename=NULL,
                       tau1=0.2,
                       tau2=0.2,
                       include_mu_argus=TRUE) {

  # Input validation -----------------------------------------------------------

  sdc_extabs_input_checks(data,
                          ID,
                          weight,
                          varpool,
                          forcelist,
                          forcenum,
                          missingdef,
                          mindim,
                          maxdim,
                          threshold,
                          wgtthreshold,
                          condition,
                          output_filename)

  # Adjust condition if needed -------------------------------------------------

  if(!is.null(condition))
    condition <- tolower(condition)

  # Save original data, before recoding missing values -------------------------

  original_data <- data

  # Make ID if none provided ---------------------------------------------------

  if(is.null(ID)) {

    ID <- ".ROW_NUMBER"
    data$.ROW_NUMBER <- 1:nrow(data)
    original_data$.ROW_NUMBER 
    warning("ID variable not provided. ID .ROW_NUMBER created.")

  } # end if 

  # Provide warning if no weight provided, set wgtthreshold to NULL ------------

  if(is.null(weight)) {

    msg <- paste("No weight variable provided.",
                 "Some risk measures cannot be calculated.")
    warning(msg)

    wgtthreshold <- NULL

  } # end if

  # If both threshold and wgtthreshold are NULL, set threshold to 3 ------------

  if(is.null(threshold) && is.null(wgtthreshold))
    threshold <- 3

  # Handle missingdef ----------------------------------------------------------

  cols <- names(missingdef)

  for(i in seq_along(cols)) {

    column <- cols[i]
    missing_values <- missingdef[[i]]

    check_condition <- !is.na(data[,column]) & data[,column] %in% missing_values
    if(any(check_condition))
      data[check_condition, column] <- NA

    if(any(trimws(missing_values) == "")) {

      check_condition <- !is.na(data[,column]) & trimws(data[,column]) == ""
      
      if(any(check_condition))
        data[check_condition, column] <- NA

    } # end if

  } # end for

  # Get tabulations ------------------------------------------------------------

  combinations <- get_all_variable_combinations(data,
                                                varpool,
                                                mindim,
                                                maxdim,
                                                forcelist,
                                                forcenum)

  tabulation <- tabulate_data(data,
                              weight,
                              combinations)

  # Flag violations in tabulation ----------------------------------------------

  violations <- flag_violations(tabulation,
                                threshold,
                                wgtthreshold,
                                condition)

  # Find which records are violations ------------------------------------------

  data_with_statistics <- violations_by_record(data, violations)

  # Calculate risk measures such as mu-argus -----------------------------------

  if(!is.null(weight)) {

    tmp <- calculate_risk_measures(data_with_statistics, 
                                   weight,
                                   varpool,
                                   tau1,
                                   tau2)
    data_with_statistics <- tmp$data

    if(include_mu_argus) {

      mu_argus_summary <- tmp$mu_argus_summary
      el_emam_measures <- tmp$el_emam_measures

    }

  } else {

    tmp <- calculate_risk_measures_unweighted(data_with_statistics,
                                              varpool,
                                              tau1)

      data_with_statistics <- tmp$data

      if(include_mu_argus)
        el_emam_measures <- tmp$el_emam_measures

  } # end if...else

  # Merge statistics with original data ----------------------------------------
    
  filter_condition <- grepl("^\\.", names(data_with_statistics))
  new_columns <- names(data_with_statistics)[filter_condition]
    
  original_data_with_statistics <- 
    merge(original_data,
          data_with_statistics[,c(ID, new_columns)],
          all.x=TRUE)

  # Save data if requested -----------------------------------------------------

  if(!is.null(output_filename)) {

    exclude <- c(".cell_count", ".pi_k", ".sum_weights", ".cell_ID")
    include_condition <- !names(original_data_with_statistics) %in% exclude
    include <- names(original_data_with_statistics)[include_condition]

    write.csv(original_data_with_statistics[,include,drop=FALSE], 
              output_filename,
              row.names=FALSE)

  } # end if

  # Get percent cell violations by dim, variable, and category -----------------

  percent_violations_by_dim_var_and_level <- 
    get_percent_cell_violations(data,
                                violations,
                                varpool,
                                mindim,
                                maxdim) 

  # Get percent violations by variable and category ----------------------------

  percent_violations_by_var_and_level <- 
    get_percent_violations_by_var_and_level(original_data_with_statistics,
                                            varpool)

  # Save options ---------------------------------------------------------------

  options <- list() 
  option_names <- c("ID", 
                    "weight", 
                    "varpool",
                    "forcelist",
                    "forcenum",
                    "missingdef",
                    "mindim",
                    "maxdim",
                    "threshold",
                    "wgtthreshold",
                    "condition",
                    "output_filename")

  for(option_name in option_names) {

    if(exists(option_name)) {

      tmp <- list(get(option_name))
      names(tmp) <- option_name

      options <- append(options, tmp)

    } # end if

  } # end for

  # Sort data sets -------------------------------------------------------------
  
  original_data_with_statistics <-
    original_data_with_statistics[order(original_data_with_statistics[[ID]]),]
  data_with_statistics <-
    data_with_statistics[order(data_with_statistics[[ID]]),]

  # Prepare return value and return --------------------------------------------

  return_value <- 
    list(tabulation=violations,
         data_with_statistics=original_data_with_statistics,
         recoded_data_with_statistics=data_with_statistics,
         percent_violations_by_var_and_level=
           percent_violations_by_var_and_level,
         percent_violations_by_dim_var_and_level=
           percent_violations_by_dim_var_and_level,
         options=options)

  if(exists("mu_argus_summary", envir=environment(), inherits=FALSE))
    return_value <- append(return_value, 
                           list(mu_argus_summary=mu_argus_summary))

  if(exists("el_emam_measures", envir=environment(), inherits=FALSE))
    return_value <- append(return_value, 
                           list(el_emam_measures=el_emam_measures))

  class(return_value) <- "sdc_extabs"

  return(return_value)

} # end 'sdc_extabs'

# ------------------------------------------------------------------------------

#' @describeIn sdc_extabs S3 print method for \code{sdc_extabs} objects
#' 
#' Prints a nicely formatted version of the percent record violations by 
#' variable/category and percent cell violations by dimension/variable/category
#' 
#' @param x An object of class \code{sdc_extabs}, as returned by the 
#'   \code{\link{sdc_extabs}} function.
#' @param cutoff The number of variable categories with the highest percentage 
#'   of cell violations for each table dimension. Default is
#'   50.
#' @param summary_outfile Name of summary output .txt file. If not NULL, console
#'   output is copied to the file. Default is NULL (no logging of output).
#'   Errors and warnings are not diverted (consider running in batch mode if
#'   logging of errors and warnings is needed).
#' @param ... Currently unused. For NextMethod compatibility.
#' @export
print.sdc_extabs <- function(x,
                             cutoff=50,
                             summary_outfile=NULL,
                             ...) {

  # Helper functions for output formatting -------------------------------------

  # Truncate string
  trunc_str <- Vectorize(function(x, ncharacters) {

    x <- as.character(x)

    if(nchar(x) > ncharacters)
      x <- substr(x, 1, ncharacters)

    return(x)

  }, vectorize.args="x")

  # ----------------------------------------------------------------------------

  format_percent <- function(x)
    return(paste0(trunc_str(formatC(100*x, mode="double", format="f", digits=2),
                            5), 
                  "%"))
  # ----------------------------------------------------------------------------

  # Return a version of x$percent_violations_by_dim_var_and_level
  # with some rows removed if number of rows in dimension exceeds cutoff
  apply_cutoff <- function() {

    tmp <- x$percent_violations_by_dim_var_and_level
    dimensions <- unique(tmp$dimension)

    new_tbl <- lapply(dimensions,
      function(dimension) {

        tbl <- tmp[tmp$dimension == dimension,]

        if(nrow(tbl) > cutoff)
          tbl <- tbl[1:cutoff,]

        return(tbl)

      }
    )

    return(do.call(rbind, new_tbl))

  } # end 'apply_cutoff'

  # ----------------------------------------------------------------------------

  make_space_str <- function(n) {

    n <- n[n > 0]

    if(length(n) > 0)
      return(paste0(rep(" ", n), collapse=""))
    else
      return("")

  }

  # Log output if requested ----------------------------------------------------

  if(!is.null(summary_outfile)) {

    stopifnot(is.character(summary_outfile))

    conn <- file(summary_outfile)
    sink(conn, append=FALSE, split=TRUE)

    on.exit(sink())
    on.exit(close(conn), add=TRUE)

  } # end if

  # Print tabulations of pre- vs. post-recoded variables.

  cat("Cross-tabulation of original vs. recoded variables.\n\n")

  for(var in x$options$varpool) {

    cat("\n", var, ":\n\n")

    old <- addNA(x$data_with_statistics[[var]])
    new <- addNA(x$recoded_data_with_statistics[[var]])

    tbl <- as.data.frame(table(old, new))
    names(tbl) <- c("Original", "Recoded", "Frequency")
    tbl <- tbl[tbl$Frequency > 0,]

    print(tbl, row.names=FALSE)

  } # end 'for'

  # Print number of records ----------------------------------------------------

  cat("\nNumber of records: ", nrow(x$data_with_statistics), "\n\n")

  # Print summary of violation counts ------------------------------------------

  vio_count_summary <- 
    c(summary(x$data_with_statistics$.violation_count),
      sum(x$data_with_statistics$.violation_count))
  names(vio_count_summary)[length(vio_count_summary)] <- "Sum"

  cat("       Summary of violation counts\n\n")
  print(vio_count_summary,digits=4)
  cat("\n\n")

  # format percentages for both tables -----------------------------------------

  x$percent_violations_by_dim_var_and_level$perc <-
    format_percent(x$percent_violations_by_dim_var_and_level$perc)

  x$percent_violations_by_var_and_level$perc <-
    format_percent(x$percent_violations_by_var_and_level$perc)

  # Get x$percent_violations_by_dim_var_and_level w/ cutoff applied 

  tmp <- apply_cutoff() 

  # Do some calculations needed for formatting ---------------------------------

  length_of_str <- 25 

  if(!is.null(x$options$threshold))
    length_of_str <- length_of_str + 
                     nchar(trunc(x$options$threshold))

  if(!is.null(x$options$wgtthreshold))
    length_of_str <- length_of_str + 6 + 
                     nchar(trunc(x$options$wgtthreshold))

  num_spaces <- floor(37 - length_of_str/2)
  space_str <- make_space_str(num_spaces)

  max_var_name <- max(nchar(as.vector(na.omit(tmp$variable))))
  num_spaces2 <- max(7, max_var_name)
  space_str2 <- make_space_str(num_spaces2)

  max_cat <- max(nchar(as.character(as.vector(na.omit(tmp$category)))))
  num_spaces3 <- max(5, max_cat)
  space_str3 <- make_space_str(num_spaces3)

  space_str4 <- make_space_str(num_spaces2 - 3)

  # Print header for 1st table -------------------------------------------------

  violation_str <- paste0("Top ",
                          cutoff,
                          " violations by variable and variable categories ",
                          "for each table dimension")

  cat(center_string(violation_str, 37), "\n")

  if(!is.null(x$options$threshold)) {

    count_string <- paste("Unweighted_count <", 
                          x$options$threshold)

    if(!is.null(x$options$wgtthreshold))
      count_string <- paste(count_string, 
                            "or Weighted count <",
                            x$options$wgtthreshold)

  } else {

    count_string <- paste("Weighted_count <", 
                          x$options$wgtthreshold)

  } # end if...else

  cat(center_string(count_string, 37), "\n\n")

  cat("Number of variables", space_str2, space_str3, "    Category of")
  cat("                Percent of Cells\n")
  cat("involved in tables", space_str3, " Variable      variable      ")
  cat("with violations based on threshold rules\n")

  names(tmp) <-
    sapply(list(rep(" ", 8), 
                rep(" ", num_spaces2 + 17), 
                rep(" ", num_spaces3 + 3), 
                rep(" ", 30)), 
           function(vec) paste(vec,collapse=""))

  # Print 1st table ------------------------------------------------------------

  print(tmp, row.names=FALSE)

  # Print header for 2nd table -------------------------------------------------

  cat("\n\n", center_string("Percent records with violations by variable and category", 37), "\n\n")

  cat(space_str4, "  Variable  Category of variable  Percent\n")

  tmp <- x$percent_violations_by_var_and_level
  names(tmp) <-
    sapply(list(rep(" ", num_spaces2 + 7), 
                rep(" ", num_spaces3 - 3), 
                rep(" ", 26)), 
           function(vec) paste(vec,collapse=""))

  # Print 2nd table ------------------------------------------------------------

  print(tmp, row.names=FALSE)

  # Handle mu-argus table if information available -----------------------------

  if("mu_argus_summary" %in% names(x)) {

    # Prepare mu-argus summary for printing ------------------------------------

    mu_argus_summary <- x$mu_argus_summary
    mu_argus_summary$total_mu_argus <- 
      formatC(mu_argus_summary$total_mu_argus, format="f", digits=4)
    mu_argus_summary$mean_mu_argus <- 
      formatC(mu_argus_summary$mean_mu_argus, format="f", digits=4)
    mu_argus_summary$total_mu_over_sum_wts <-
      formatC(mu_argus_summary$total_mu_over_sum_wts, format="f", digits=9)
    mu_argus_summary$total_mu_over_total_obs <-
      formatC(mu_argus_summary$total_mu_over_total_obs, format="f", digits=9)
    mu_argus_summary$max_count <- 
      c(" =1", "<=2", "<=3", "ALL")[1:nrow(mu_argus_summary)]
    names(mu_argus_summary) <- c("      ",
                                 "         ",
                                 "                ",
                                 "             ",
                                 "             ",
                                 "             ")

    # Print mu-argus table -----------------------------------------------------

    cat("\n\n                       MU-ARGUS summaries\n\n")
    cat("                                                    Total Argus  Total Argus\n")
    cat("                                                    score/total score/overall\n")
    cat("   Cell    Number      Total Argus    Mean Argus     number of      sum of\n")
    cat("   count  of cases           score         score   observations    weights\n")
    print(mu_argus_summary, row.names=FALSE)

  } # end if

  # Handle el-emam table -------------------------------------------------------

  if("el_emam_measures" %in% names(x)) {

    el_emam <- as.data.frame(x$el_emam_measures)

    for(i in 1:ncol(el_emam))
      el_emam[,i] <- formatC(signif(el_emam[,i], digits=7), 
                             format="fg",
                             flag="#")

    cat("\n\n")

    if(ncol(el_emam) > 3)
      cat("   ")

    cat("   Re-identification Risk Metrics (El-Emam)\n\n")

    print(as.data.frame(el_emam), row.names=FALSE)

  } # end if
  
  # Print top 10 violations ----------------------------------------------------

  varpool <- x$options$varpool
  dat <- x$data_with_statistics
  id_var <- x$options$ID
  rows_to_print <- order(dat$.violation_count, decreasing=TRUE)[1:10]

  if("mu_argus_summary" %in% names(x)) {

    weight <- x$options$weight

    orig_columns <- c(id_var, varpool, weight, ".mu_argus", ".violation_count")
    new_columns <- c(id_var, 
                     varpool, 
                     weight, 
                     "Mu-Argus Score", 
                     "Violation Count")
    top_violations <- dat[rows_to_print, orig_columns]
    colnames(top_violations) <- new_columns
      
  } else {

    orig_columns <- c(id_var, varpool, ".violation_count")
    new_columns <- c(id_var, varpool, "Violation Count")
    top_violations <- dat[rows_to_print, orig_columns]
    colnames(top_violations) <- new_columns
      
  } # end if...else

  cat("\n\nTop 10 Records with Most Violations\n\n")
  print(top_violations, row.names=FALSE)

  cat("\n") # Extra space, in case of warning messages

  return(invisible(NULL))

} # end 'print.sdc_extabs'

# ------------------------------------------------------------------------------

#' @describeIn sdc_extabs S3 plot method for \code{sdc_extabs} objects
#' 
#' Produces boxplots and scatterplots of violation counts and mu-argus scores.
#' 
#' @param plotpath Directory to save plots. Plots are saved as \emph{jpeg} files 
#'   (quality = 100\%). If the directory does not exist, it is first created.
#'   If \code{plotpath} is NULL (default), plots are not saved.
#' @param plotvar1 A vector of names of discrete variables for boxplots. If 
#'   none, boxplots are not produced.
#' @param plotvar2 A vector of names of continuous variables for scatterplots. 
#'   If none, scatterplots are not produced.
#' @export
plot.sdc_extabs <- function(x,
                            plotpath=NULL,
                            plotvar1=character(0),
                            plotvar2=character(0),
                            ...) {

  dat <- x$data_with_statistics
  has_mu_argus <- ".mu_argus" %in% names(dat)

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

  get_plot_filename <- function(plot_type, is_vio_count, variable) {

    fname <- paste0(plot_type,
                    "--", 
                    ifelse(is_vio_count, "violation_count--", "mu_argus--"),
                    variable,
                    ".jpeg")

    return(file.path(plotpath, fname))

  } # end 'get_plot_filename'

  # ----------------------------------------------------------------------------

  # Call before making any plot
  plot_start <- function(plot_type, is_vio_count, variable) {

    if(save_plots)
      jpeg(get_plot_filename(plot_type, is_vio_count, variable), quality=100)
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
  
  make_violation_count_boxplot <- function(plotvar) {

    plot_start("box", is_vio_count=TRUE, plotvar)
    boxplot_formula <- as.formula(paste(".violation_count ~", plotvar))
    boxplot(boxplot_formula, 
            data=dat,
            ylab="Violation count")
    plot_end()

    return(invisible(NULL))

  } # end 'make_violation_count_boxplot'

  # ----------------------------------------------------------------------------
  
  make_mu_argus_boxplot <- function(plotvar) {

    if(has_mu_argus) {

      plot_start("box", is_vio_count=FALSE, plotvar)
      boxplot_formula <- as.formula(paste(".mu_argus ~", plotvar))
      boxplot(boxplot_formula, 
              data=dat,
              ylab="Mu-Argus Score")
      plot_end()

    }

    return(invisible(NULL))

  } # end 'make_mu_argus_boxplot'

  # ----------------------------------------------------------------------------
  
  make_violation_count_scatterplot <- function(plotvar) {

    plot_start("scatter", is_vio_count=TRUE, plotvar)
    scatrplt <- 
      ggplot2::ggplot(dat, 
                      mapping=ggplot2::aes_string(x=plotvar, 
                                                  y=".violation_count")) + 
      ggplot2::geom_point() +
      ggplot2::ylab("Violation Count")
    plot(scatrplt)
    plot_end()

    return(invisible(NULL))

  } # end 'make_violation_count_scatterplot'

  # ----------------------------------------------------------------------------
  
  make_mu_argus_scatterplot <- function(plotvar) {

    if(has_mu_argus) {

      plot_start("scatter", is_vio_count=FALSE, plotvar)
      scatrplt <- 
        ggplot2::ggplot(dat, 
                        mapping=ggplot2::aes_string(x=plotvar, 
                                                    y=".mu_argus")) + 
        ggplot2::geom_point() +
        ggplot2::ylab("Mu-Argus Score")
      plot(scatrplt)
      plot_end()

    }

    return(invisible(NULL))

  } # end 'make_mu_argus_scatterplot'

  # ----------------------------------------------------------------------------
  
  # Create output directory (if any) if it does not already exist --------------

  # Clean up directory name
  if(!is.null(plotpath))
    plotpath <- file.path(plotpath)

  if(!is.null(plotpath) && !dir.exists(plotpath))
    dir.create(plotpath, recursive=TRUE)

  # Produce violation count boxplots -------------------------------------------

  lapply(plotvar1, make_violation_count_boxplot)

  # Produce violation count scatterplots ---------------------------------------

  lapply(plotvar2, make_violation_count_scatterplot)

  # Produce Mu-Argus boxplots --------------------------------------------------

  lapply(plotvar1, make_mu_argus_boxplot)

  # Produce mu_argus scatterplots ----------------------------------------------

  lapply(plotvar2, make_mu_argus_scatterplot)

  return(invisible(NULL))

} # end 'plot.sdc_extabs'

# ------------------------------------------------------------------------------

# Input validation
sdc_extabs_input_checks <- function(data, 
                                    ID,
                                    weight,
                                    varpool,
                                    forcelist,
                                    forcenum,
                                    missingdef,
                                    mindim,
                                    maxdim,
                                    threshold,
                                    wgtthreshold,
                                    condition,
                                    output_filename) { 

  stopifnot(inherits(data, "data.frame"))
  stopifnot(nrow(data) > 0)

  stopifnot(is.null(ID) ||
            (is.character(ID) && all(ID %in% names(data))))


  stopifnot(is.null(weight) || 
            (is.character(weight) && 
             length(weight) == 1 && 
             weight %in% names(data)))

  stopifnot(is.character(varpool))
  stopifnot(length(varpool) > 0)
  stopifnot(all(varpool %in% names(data)))

  stopifnot(inherits(forcelist, "character"))
  stopifnot(all(forcelist %in% varpool))

  stopifnot(length(forcelist) == 0 || 
            (inherits(forcenum, c("integer", "numeric")) &&
             forcenum == as.integer(forcenum) &&
             forcenum >= 0L &&
             forcenum <= length(forcelist)))

  stopifnot(all(forcelist %in% varpool))

  stopifnot(inherits(missingdef, "list"))
  stopifnot(length(missingdef) == 0 || !is.null(names(missingdef)))
  stopifnot(is.null(names(missingdef)) || 
            all(names(missingdef) %in% names(data)))

  stopifnot(inherits(mindim, c("integer", "numeric")))
  stopifnot(length(mindim) == 1)
  stopifnot(mindim == as.integer(mindim))
  stopifnot(mindim > 0)
  stopifnot(mindim <= length(varpool))

  stopifnot(length(forcelist) == 0 || 
            forcenum <= mindim)

  stopifnot(inherits(maxdim, c("integer", "numeric")))
  stopifnot(length(maxdim) == 1)
  stopifnot(maxdim == as.integer(maxdim))
  stopifnot(maxdim > 0)
  stopifnot(maxdim <= length(varpool))
  stopifnot(maxdim >= mindim)

  stopifnot(is.null(threshold) || 
            inherits(threshold, c("integer", "numeric")))

  stopifnot(is.null(wgtthreshold) ||
            inherits(wgtthreshold, c("integer", "numeric")))
  
  stopifnot(is.null(condition) || 
            (is.character(condition) &&
             length(condition) == 1 && 
             tolower(condition) %in% c("and", "or")))

  stopifnot(is.null(output_filename) || 
            (is.character(output_filename) && 
             length(output_filename) == 1))

} # end 'sdc_extabs_input_checks'

# ------------------------------------------------------------------------------

center_string <- function(str, midpoint, print_str=FALSE) {

  num_spaces <- midpoint - floor(nchar(str)/2)

  if(num_spaces > 0)
    new_str <- paste(paste(rep(" ", num_spaces), collapse=""), str)
  else
    new_str <- str

  if(print_str)
    print(new_str)

  return(new_str)

} # end 'center_string'