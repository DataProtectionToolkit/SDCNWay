tabulate_data <- function(data,
                          weight=NULL,
                          combinations) {
# Calculates cell counts.
#
# Calculates 1way, 2way, 3way, etc. tabulations of variables.
#
# data: Data frame containing the data to tabulation
# weight: Column name for weighting variable. NULL  if none.
# combinations: All allowed combinations of variables, as returned by
#   get_all_variable_combinations.
# return value: A list of tabulations
  tabulations <- lapply(combinations, 
                        function(combo) get_table(data, combo, weight))
  tabulations <- tabulations[sapply(tabulations, 
                                    function(tab) !inherits(tab, "character"))]

  return(tabulations)

} # end 'tabulate_data'

# ------------------------------------------------------------------------------

get_all_variable_combinations <- function(data,
                                          varpool,
                                          mindim=1,
                                          maxdim=2,
                                          forcelist=character(0),
                                          forcenum=1) {
# Calculate all allowed combinations of variables in varpool
# 
# data: Data set. Needed to check missing values.
# varpool: Vector of column names over which to form tables.
# mindim: Integer specifies the minimum number of varpool variables 
#   (including forcelist variables) that can be used to form tables. 
# maxdim: Integer specifies the maximum number of varpool variables
#   (including forcelist variables) that can be used to form tables. 
# forcelist Vector of variables that are mandatory for all tabulations. 
# forcenum Number of variables in forcelist that are mandatory for all
#   tabulations.
  
  # Remove variables from varpool with all missing values ----------------------

  any_not_missing <- function(variable)
    return(any(!is.na(data[,variable])))

  varpool <- varpool[sapply(varpool, any_not_missing)]

  # Find all combinations of the appropriate dimensions ------------------------

  combinations <- lapply(mindim:maxdim, 
    function(dimension) 
      as.list(as.data.frame(combn(varpool, dimension), 
                            stringsAsFactors=FALSE)))
  combinations <- unlist(combinations, recursive=FALSE)

  if(length(forcelist) > 1) {

    is_okay <- sapply(combinations,
                      function(combo)
                       return(sum(forcelist %in% combo) == forcenum))
    combinations <- combinations[is_okay]

  } # end if

  return(combinations)
  
} # end 'get_all_variable_combinations'

# ------------------------------------------------------------------------------

# Add ID variable, .cell_ID, to data. Two records will have the same cell ID if
# and only if they are identical across all variables in varpool.
add_cell_ID <- function(data, varpool=names(data)) {

  udata <- unique(data[,varpool,drop=FALSE])
  udata$.cell_ID <- 1:nrow(udata)

  return(merge(data, udata))

} # end 'add_cell_ID'

# ------------------------------------------------------------------------------

flag_violations <- function(tabulation,
                            threshold,
                            wgtthreshold=NULL,
                            condition=NULL) {
# Identify cells with violations
#
# tabulation: List returned by tabulate_data function.
# threshold: If number of cases in cell less than threshold, flag cell
#   as a violation
# wgtthreshold: Threshold to determine violations in terms of weighted
#   counts.
# condition: Character variable describing how weighted and unweighted
#   thresholds are combined when both are used. If used, must be "and" or "or".
# combinations: Variable combinations as returned by 
#   get_all_variable_combinations
# return value: A list of containing values with violations

  # Helper function ------------------------------------------------------------

  flag_violations_one_table <- function(tbl) {

    if(!is.null(threshold) && ".unweighted_freq" %in% names(tbl))
      is_violation <- tbl$.unweighted_freq > 0 &
                      tbl$.unweighted_freq < threshold

    if(!is.null(wgtthreshold) && ".weighted_freq" %in% names(tbl)) {

      tmp_violation <- (tbl$.weighted_freq > 0 & 
                        tbl$.weighted_freq < wgtthreshold)

      if(is.null(threshold))
        is_violation <- tmp_violation
      else if(condition == "and")
        is_violation <- is_violation & tmp_violation
      else
        is_violation <- is_violation | tmp_violation         

    } # end if

    tbl$.violation <- is_violation

    return(tbl)

  } # end 'flag_violations_one_table'

  # Apply helper to each table -------------------------------------------------

  return(lapply(tabulation, 
                flag_violations_one_table))

} # end 'flag_violations'

# ------------------------------------------------------------------------------

get_percent_cell_violations <- function(data,
                                        tabs_with_flags,
                                        varpool,
                                        mindim,
                                        maxdim) {
# Get percent cell violations for each dimension/variable/category
# 
# data: Data frame containing the original data
# tabs_with_flags: Tabulations with violation flag attached, as returned
#   by flag_violations
# varpool: Variable list of variables used for tabulations. Same list
#   as provided to sdc_extabs.
# mindim: Minimum dimension. Same as provided to sdc_extabs. 
# maxdim: Maximum dimension. Same as provided to sdc_extabs. 

  # Determine if data is weighted ----------------------------------------------

  has_weighted <- ".weighted_freq" %in% names(tabs_with_flags[[1]])

  # Helper functions -----------------------------------------------------------

  initialize_perc_table <- function() {

    dimensions <- mindim:maxdim 

    perc_tbl <- data.frame(dimension=integer(0), 
                           variable=character(0),
                           category=character(0))

    tmp <- lapply(varpool,
      function(variable) {

      values <- as.character(na.omit(unique(data[,variable])))

      tmp_tbl <- expand.grid(dimension=dimensions, category=values)
      tmp_tbl$variable <- variable

      return(tmp_tbl)

    }) # end lapply

    perc_tbl <- do.call(rbind, tmp)

    return(perc_tbl)

  } # end 'initialize_perc_table'

  # ----------------------------------------------------------------------------

  get_perc_violation <- function(dimension, variable, category) {

    ncols_desired <- dimension + 2 + has_weighted

    # Filter tables of the right dimension -------------------------------------

    condition <- sapply(tabs_with_flags, 
                        function(tab) return(ncol(tab) == ncols_desired))
    tabs <- tabs_with_flags[condition]

    # Filter tables containing this variable -----------------------------------

    condition <- sapply(tabs, function(tab) return(variable %in% names(tab)))
    tabs <- tabs[condition]

    # Initialize tallies -------------------------------------------------------

    violation_tally <- 0
    nonviolation_tally <- 0

    # Loop through tables, adding to tallies -----------------------------------

    for(tab in tabs) {

      # Filter rows with the correct category ----------------------------------

      tmp <- tab[as.character(tab[[variable]]) == as.character(category),]

      # Tally violation and non-violation cells --------------------------------

      if(nrow(tmp) > 0) {

        violation_tally <- violation_tally + sum(tmp$.violation)
        nonviolation_tally <- nonviolation_tally + sum(!tmp$.violation)

      }

    } # end for

    return(violation_tally/(violation_tally + nonviolation_tally))

  } # end 'get_perc_violation'

  # ----------------------------------------------------------------------------

  sort_perc_table <- function(tbl)
    return(dplyr::arrange(tbl, 
                          dimension,
                          dplyr::desc(perc), 
                          variable,
                          as.character(category)))

  # Initialize table -----------------------------------------------------------

  perc_table <- initialize_perc_table()
  
  # Convert varpool variables to factors. Make NAs explicit --------------------

  for(j in length(varpool))
    data[[varpool[j]]] <- addNA(data[[varpool[j]]], ifany=TRUE)

  # Acquire violation percentages ----------------------------------------------

  perc_table$perc <- with(perc_table, mapply(get_perc_violation, 
                                             dimension,
                                             variable,
                                             category))

  # Filter out rows where perc is NaN (can happen when forcelist is used). -----

  perc_table <- perc_table[!is.nan(perc_table$perc),]

  # Filter rows with violaton percentages > 0 ----------------------------------

  perc_table <- perc_table[perc_table$perc > 0,]
  
  # Sort table by dimension, percentage (descending), variable, and category, --
  # in that order --------------------------------------------------------------

  perc_table <- sort_perc_table(perc_table)

  # Reorder columns and return -------------------------------------------------

  return(perc_table[,c("dimension", "variable", "category", "perc")])

} # end 'get_percent_cell_violations'

# ------------------------------------------------------------------------------

get_percent_violations_by_var_and_level <- function(data_with_vio_counts,
                                                    varpool) {
# Get percent record violations by variable/category
# 
# data_with_vio_counts: Data with a .violation_count column. As returned
#   by violations_by_record.
# varpool: Variable list of variables used for tabulations. Same list
#   as provided to sdc_extabs.

  # Helper function ------------------------------------------------------------

  get_perc_violation <- function(x) {
  
    variable <- x["variable"]
    category <- x["category"]

    is_corresponding_record <- !is.na(data_with_vio_counts[[variable]]) &
                               as.character(data_with_vio_counts[[variable]]) ==
                                 category

    return(mean(data_with_vio_counts$.violation[is_corresponding_record]))

  } # end 'get_perc_violation'

  # Set violation flag ---------------------------------------------------------

  data_with_vio_counts$.violation <- data_with_vio_counts$.violation_count > 0

  # Get list of variables and categories ---------------------------------------

  vars <- character(0)
  levs <- character(0)

  for(var in varpool) {

    columntype <- class(data_with_vio_counts[[var]])
    unique_levels <- as(unique(na.omit(data_with_vio_counts[[var]])), 
                        columntype)

    vars <- c(vars, rep(var, length(unique_levels)))
    levs <- c(levs, as.character(unique_levels))

  } # end for

  perc_vio_tbl <- data.frame(variable=vars, 
                             category=levs,
                             stringsAsFactors=FALSE)

  # Calculate percent of violations in cells -----------------------------------

  perc_vio_tbl$perc <- apply(perc_vio_tbl, 1, get_perc_violation)
  perc_vio_tbl <- perc_vio_tbl[perc_vio_tbl$perc > 0,]

  # Reformat table (sort, remove row names) and return -------------------------

  order_of_rows <- with(perc_vio_tbl, order(variable, 
                                            as.character(category), 
                                            perc))
  perc_vio_tbl <- perc_vio_tbl[order_of_rows,]
  row.names(perc_vio_tbl) <- NULL

  return(perc_vio_tbl)

} # end 'get_percent_violations_by_var_and_level'

# ------------------------------------------------------------------------------

violations_by_record <- function(data, tbls_with_violation_flag) {
# Adds violation counts to data.
# 
# data: Data frame containing the original data
# tbls_with_violation_flag: Tabulations with violation flag attached, as 
#   returned by flag_violations
  
  results <- data
  results$.violation_count <- 0

  for(i in seq_along(tbls_with_violation_flag)) {

    tbl <- tbls_with_violation_flag[[i]]

    is_var_name <- !names(tbl) %in% c(paste0(".", c("", "un"), "weighted_freq"), 
                                      ".violation")
    var_names <- names(tbl)[is_var_name]

    # 'merge' function messes up row order, which causes problems in next line.
    # So, we use plyr::join instead.
    tmp <- plyr::join(data, tbl, by=var_names)
  
    if(any(is.na(tmp$.violation)))
      tmp$.violation[is.na(tmp$.violation)] <- FALSE

    results$.violation_count <- results$.violation_count + tmp$.violation

  } # end for

  return(results)

} # end 'violations_by_record'
