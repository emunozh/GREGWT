# Created by Esteban Munoz (emunozh@gmail.com).
#
# 04.02.2015
# last edit: 20.04.2015
# Mon 04 Apr 2016 02:16:15 PM CEST
#


#' @title prepareData
#'
#' @description
#' Prepares the data for simulation
#'
#' @param census Census data of small areas.
#' @param survey A survey of individual records (microdata).
#' @param census_area_id (optional, default=1) row name or row index with
#' area id in the census data. Define as `FALSE` if area code should be
#' generated.
#' @param census_categories (optional, default=FALSE) row names or row index
#' of with categories to be used in the simulation.
#' @param survey_weights (optional, default=FALSE) row name or row index of
#' initial weights in the survey data. `FALSE` will use the last column.
#' @param survey_id (optional, default=1) individual records id's. Define
#' as `FALSE` to generate an id.
#' @param survey_categories (optional, default=FALSE) survey categories to be
#' used in the simulation.
#' @param reference_col (optional, default=FALSE) Category used as reference.
#' @param group (optional, default=FALSE) Used variable to run an integrated
#' re-weighting simulation.
#' @param convert (optional, default=TRUE) Converts data to binary format.
#' @param na.rm (optional, default=FALSE) remove records with nan values.
#' @param use_base (optional, default=TRUE) use the model.matrix function form
#' base R.
#' @param breaks (optional, default=FALSE) define the beaks to calculate
#' population totals, if FALSE population totals won't be computed
#' @param pop_benchmark (optional, default=FALSE) define the benchmark to be
#' used for the computation of the total population, pass as a vector/
#' containing the breaks of the benchmark (e.g. \code{pop_benchmark=c(1,5)}).
#' If FALSE the function will compute total population as the mean of the all
#' benchmarks.
#' @param pop_du (optional, default=FALSE) define the benchmark to be used for
#' the computation of total dwelling units. Analog to
#' \code{pop_benchmark}
#' @param pop_building (optional, default=FALSE) define the benchmark to be
#' used for the computation of total building units. Analog to
#' \code{pop_benchmark}
#' @param pop_total_col (optional, default=FALSE) col containing the population
#' totals
#' @param align (optional, default=FALSE) align values to population totals
#' @param verbose (optional, default=FALSE) be verbose
#' @return X Prepared survey matrix.
#' @return Tx Marginal totals for simulation area.
#' @return dx Survey design weights.
#' @return area_id Small area ID.
#' @return total_pop mean population totals for each area
#' @return X_complete binary formatted survey with all all categories
#' @return Tx_complete marginal sums with all categories
#' @examples
#' data("GREGWT.census")
#' data("GREGWT.survey")
#'
#' simulation_data <- prepareData(GREGWT.census, GREGWT.survey,
#'                                census_categories=seq(2,24),
#'                                survey_categories=seq(1,3))
#'
#' simulation_data1 <- prepareData(GREGWT.census, GREGWT.survey,
#'                                 census_categories=seq(2,24),
#'                                 survey_categories=seq(1,3),
#'                                 pop_benchmark=c(2,12),
#'                                 verbose=TRUE)
#'
#' # compute the total population as the mean of all benchmarks. Breaks parameters
#' # needs to be defined. In this case the breaks are displaced by one because the
#' # area code is on the first column.
#' simulation_data2 <- prepareData(GREGWT.census, GREGWT.survey,
#'                                 census_categories=seq(2,24),
#'                                 survey_categories=seq(1,3),
#'                                 breaks=c(11, 17),
#'                                 verbose=TRUE)
#'
#' total_pop1 <- simulation_data1$total_pop
#' plot(total_pop1$pop)
#' total_pop2 <- simulation_data2$total_pop
#' points(total_pop2$pop, col="red", pch="+")
#'
#' @author M. Esteban Munoz H.
prepareData <- function(census, survey,            # require input data
    census_area_id     = 1,     # census area id
    survey_id          = 1,     # survey id
    convert            = TRUE,  # convert data to binary
    use_base           = TRUE,  # use model.matrix
    census_categories  = FALSE, # census categories
    survey_weights     = FALSE, # survey weights
    survey_categories  = FALSE, # survey categories
    reference_col      = FALSE, # define a reference column
    group              = FALSE, # if iterated reweighting
    na.rm              = FALSE, # remove nan
    breaks             = FALSE, # compute population totals
    pop_benchmark      = FALSE, # benchmark to be use for pop
    du_benchmark       = FALSE, # benchmark to be use for du
    building_benchmark = FALSE, # benchmark to be use for building
    align              = FALSE, # align to total population
    pop_total_col      = FALSE, # population totals
    verbose            = FALSE  # be verbose
                        ){
    if (verbose){
        cat("preparing data --> ")
        cat("\ndim(survey):", dim(survey))
        cat("\ndim(census):", dim(census))
    }

    if (class(group) != "logical"){
        group_id <- survey[group]
        survey <- survey[!(names(survey) %in% group)]
    }

    ### Format data census
    # get id pos pop
    area <- getAreaid(census,
                      census_area_id     = census_area_id,
                      breaks             = breaks,
                      pop_benchmark      = pop_benchmark,
                      du_benchmark       = du_benchmark,
                      building_benchmark = building_benchmark,
                      pop_total_col      = pop_total_col,
                      verbose            = verbose
                      )

    ### Filter data
    filtered <- filterData(census, census_categories, census_area_id,
                           survey, survey_categories, survey_id,
                           survey_weights, na.rm,
                           verbose=verbose)
    Tx <- filtered$Tx
    dx <- filtered$dx
    survey_id <- filtered$survey_id
    X <- filtered$X

    ### Align data
    Tx <- alignData(Tx, align, breaks, area$pop, verbose=verbose)

    ### Transform categorical data to binary columns
    if (!(is.logical(reference_col))) use_base <- FALSE
    if (convert){
        if(verbose) cat("\n\t|--> categories --> ")
        X_old_temp <- X
        X_complete <- toBinary(X, reference_col=T, use_base=use_base,
                               verbose=verbose)
        if(verbose) cat(" X(complete) ")
        X <- toBinary(X, reference_col=reference_col, use_base=use_base,
                      verbose=verbose)
        if(verbose) cat(" X ")
        keep_cols <- rownames(X_old_temp) %in% rownames(X)
        dx <- dx[keep_cols]
        survey_id <- survey_id[keep_cols, ]
        if(verbose) cat(" ok")
    } else {
        X_complete <- X
    }

    ### Add a suffix to data names
    colnames(X) <- paste("G", colnames(X), sep = ".")
    colnames(X_complete) <- paste("G", colnames(X_complete), sep = ".")
    # Add a suffix to data names
    colnames(Tx) <- paste("G", colnames(Tx), sep = ".")

    ### Convert matrix to data frame
    Tx <- as.data.frame(Tx)
    Tx_complete <- Tx[names(Tx) %in% colnames(X_complete)]
    Tx <- Tx[names(Tx) %in% colnames(X)]

    ### Check data dimensions (1) number of constrains
    if(verbose){
        cat("\n\t|--> dimensions and col names -->")
        cat(" dim(X):", dim(X))
        cat(" dim(Tx)[2]:", dim(Tx)[2])
        cat(" length(dx):", length(dx))
    }
    if(dim(X)[2] != dim(Tx)[2]){
        namesdX <- namesDiff(X, Tx)
        namesdTx <- namesDiff(Tx, X)
        stop(paste(
            "\ndimension error",
            "\ndim Tx:", dim(Tx)[2],
            "\n\t Tx names not in X : ",
            namesdTx,
            "\ndim  X:", dim(X)[2],
            "\n",
            # paste(names(X), collapse="; "),
            "\n\t X  names not in Tx: ",
            namesdX))
    }else{ if(verbose) cat(" dim1:ok")}
    # Check data dimensions (2) weights - survey
    if(dim(X)[1] != length(dx)){
        stop(paste(
            "\ndimension error",
            "\nweights length:", length(dx),
            "\nsurvey  length:", dim(X)[1]))
    }else{ if(verbose) cat(" dim2:ok ")}
    # Check data dimensions (3) Tx/X are not empty
    if(dim(Tx)[1] == 0){stop("Tx is empty")}
    if(dim(X)[1] == 0){stop("X is empty")}
    if(verbose) cat(" ok")

    # Add group variable if defined
    Xt <- X
    if(class(group) != "logical"){
        X <- cbind(X, group_id)
        X_complete <- cbind(X_complete, group_id)
    }

    ### Sort variables by name
    X <- X[, order(colnames(X))]
    Tx <- Tx[, order(colnames(Tx))]

    # Save original data set
    Xin <- X

    # remove columns of 1 and 0
    Xo  <-  X[colSums(X,na.rm=TRUE)  != 0 & colSums(X,na.rm=TRUE)  != dim(X)[1]]
    Txo <- Tx[colSums(Xt,na.rm=TRUE) != 0 & colSums(Xt,na.rm=TRUE) != dim(Xt)[1]]
    removed.var <- dim(X)[2] - dim(Xo)[2]
    if (removed.var>=1) { if (verbose){cat("removed:", removed.var) } }

    # Add survey id
    Xo  <- cbind(survey_id, Xo)

    # Convert data to matrix
    if (class(Xo) == "data.frame"){ Xo <- data.matrix(Xo) }
    Xo  <- as.matrix(Xo)
    #
    if (class(Txo) == "data.frame"){ Txo <- data.matrix(Txo) }
    Txo <- as.matrix(Txo)

    if (verbose) cat("data:OK\n")
    prepared.data <- list(X           = Xo,
                          Tx          = Txo,
                          dx          = dx,
                          breaks      = breaks,
                          X_complete  = X_complete,
                          Tx_complete = Tx_complete,
                          area_id     = area$id,
                          total_pop   = area$pop,
                          survey      = survey)
    class(prepared.data) <- "prepareData"
    return(prepared.data)
}


alignData <- function(Tx, align, breaks, area, verbose=F){
    if(!(class(align) == "logical")){
        if(verbose) cat("\n\t|--> align data -->")
        if(class(breaks) == "logical"){
            stop("in order to align the data I need to know the category breaks")
        }else{
            breaks <- c(1, breaks, dim(Tx)[2])
            if (verbose) printBreaks(names(Tx), breaks)
            for(i in seq(dim(align)[2])){
                vara <- names(align)[i]
                sums <- area[vara]
                if (length(sums) >= 1 && sums != 0){
                    breaks_sum = breaks[align[1,i]:align[2,i]]
                    if (verbose) printBreaks(names(Tx), breaks_sum, name=vara)
                    pos = 0
                    if(i > 1) pos = 1
                    Tx <- alignWithTotals(Tx, breaks_sum, sums, pos,
                                          name=vara,
                                          verbose=verbose)
                }
            }
        }
    }
    return(Tx)
}


printBreaks <- function(names_tx, breaks, name="all"){
    cat("\n\t\tBreaks: ", name, "\n")
    for (br in breaks){
        cat("\n\t\t", br, "-->", names_tx[br])
    }
    cat("\n")
}


#' @title getAreaid
#'
#' @description
#' Gets the area id from the survey data or generates a new id for the survey.
#' This function also calculates the population size of each area.
#'
#' @param census binary census survey data
#' @param census.area.if (optional, default=False) position or name of the area
#' id in the survey. If False a new id vector will be generated for the survey.
#' @param breaks (optional, default=FALSE) breaks of the categories, needed to
#' estimate the total population of the areas. If FALSE the function will not
#' compute the total population for the simulation areas.
#' @param pop_benchmark (optional, default=FALSE) define the benchmark to be
#' used for the computation of the total population, pass as a vector
#' containing the breaks of the benchmark (e.g. \code{pop_benchmark=c(1,5)}).
#' If FALSE the function will compute total population as the mean of the all
#' benchmarks.
#' @param pop_total_col (optional, default=FALSE) col containing the population
#' totals
#' @param verbose (optional, default=FALSE) be verbose
#' @return list(id, pos, pop)
#' @author M. Esteban Munoz H.
#TODO: examples
getAreaid <- function(census,
                      census_area_id     = FALSE,
                      breaks             = FALSE,
                      pop_benchmark      = FALSE,
                      du_benchmark       = FALSE,
                      building_benchmark = FALSE,
                      pop_total_col      = FALSE,
                      verbose            = FALSE
                      ){
    if(verbose) cat("\n\t|--> get area id --> ")
    if(class(census_area_id) %in% c("numeric", "integer")){
        if (verbose) cat("by position (", names(census)[census_area_id], ")")
        area_id     <- census[census_area_id]
        area_id_pos <- census_area_id
    }else if(is.logical(census_area_id)){
        if (verbose) cat("created new")
        area_id     <- as.data.frame(seq(dim(census)[1]))
        area_id_pos <- census_area_id
    }else{
        if (verbose) cat("by name (", census_area_id, ")")
        area_id     <- census[names(census) %in% census_area_id]
        area_id_pos <- which(names(census) == census_area_id)
    }
    # Compute population totals
    if(is.logical(breaks) & is.logical(pop_benchmark) & is.logical(pop_total_col)){
        total_pop <- FALSE
    }else{
        total_pop <- getPopulation(census, breaks,
                                   area_id            = area_id,
                                   area_id_pos        = area_id_pos,
                                   pop_benchmark      = pop_benchmark,
                                   du_benchmark       = du_benchmark,
                                   building_benchmark = building_benchmark,
                                   pop_total_col      = pop_total_col,
                                   verbose            = verbose
                                   )}
    if(verbose) cat("ok ")
    list(id=area_id, pos=area_id_pos, pop=total_pop)
}


#' @title getPopulation
#'
#' @description
#' Computes the total population of simulation areas using the constrain
#' variables.
#'
#' @param Tx census data
#' @param breaks positions of the variable breaks in the census
#' @param area_id_pos (optional, default=FALSE) position of the area id on the
#' census. If FALSE the function will try to use \code{area_id} as area codes,
#' if \code{area_id} is also FALSE the function will issue an error.
#' @param area_id (optional, default=FALSE) vector of area codes to be used. If
#' If FALSE the function will try to use \code{area_id} to identify the column
#' in the census containing the area codes, if \code{area_id_pos} is also FALSE
#' the function will issue an error.
#' @param pop_benchmark (optional, default=FALSE) define the benchmark to be
#' used for the computation of the total population, pass as a vector
#' containing the breaks of the benchmark (e.g. \code{pop_benchmark=c(1,5)}).
#' If FALSE the function will compute total population as the mean of the all
#' benchmarks.
#' @param pop_total_col (optional, default=FALSE) col containing the population
#' totals
#' @return total_pop population total for each simulation area
#' @examples
#' data("GREGWT.census")
#' total_pop <- getPopulation(GREGWT.census, c(11,17), area_id_pos=1)
#' @author M. Esteban Munoz H.
getPopulation <- function(Tx, breaks,
                          area_id_pos        = FALSE,
                          area_id            = FALSE,
                          pop_benchmark      = FALSE,
                          du_benchmark       = FALSE,
                          building_benchmark = FALSE,
                          pop_total_col      = FALSE,
                          verbose            = FALSE
                          ){
    Tx_for_pop <- Tx
    if (verbose) cat("\n\t\t|--> get population totals --> ")
    if (!(is.logical(area_id_pos))) {
        if (verbose) cat("\n\t\t\t|--> got area id position --> ")
        area_id_dat <- Tx[ area_id_pos]
        Tx          <- Tx[-area_id_pos]
        if (verbose) cat("ok ")
    } else if (is.logical(area_id_pos) & is.logical(area_id)){
        stop("To compute the population total this function needs either:
    (a) the position of the area id area_id_pos=int; or
    (b) a vector with the area id area_id=vector()")
    } else if (!(is.logical(area_id))){
        area_id_dat <- area_id
        if (verbose) cat(" ok ")
    }

    pop_sums_t <- FALSE
    if (class(pop_total_col) != "logical"){
        if (verbose) cat("\n\t\t\t|--> got total population column --> ")
        pop_sums_t <- rowSums(Tx[pop_total_col])
        if (verbose) cat("ok ")
    }

    if (class(building_benchmark) != "logical"){
        building_sums <- computeTotal(Tx_for_pop,
                                      building_benchmark,
                                      name    = "bui",
                                      verbose = verbose)
    } else {
        building_sums <- FALSE
    }

    if (class(du_benchmark) != "logical"){
        du_sums <- computeTotal(Tx_for_pop,
                                du_benchmark,
                                name    = "hhd",
                                verbose = verbose)
    } else {
        du_sums <- FALSE
    }

    if (sum(is.na(pop_sums_t)) >= 1 | is.logical(pop_sums_t)){
        if(class(pop_benchmark) != "logical"){
            pop_sums <- computeTotal(Tx_for_pop,
                                     pop_benchmark,
                                     name    = "ind",
                                     verbose = verbose)
        } else if (is.logical(du_benchmark) & is.logical(building_benchmark)){
            if (verbose) cat("\n\t\t\t|--> got breaks --> ")
            breaks <- c(breaks, dim(Tx)[2])
            pop_sums <- vector(length=dim(Tx)[1])
            i = 1; j = 0
            for (b in breaks){
                pop_sums <- pop_sums + rowSums(Tx[names(Tx)[i:b]], na.rm=TRUE)
                i = b+1; j = j+1}
                pop_sums = pop_sums / j
            if (verbose) cat(" ok ")
        } else {
            cat("\nWARNING! no total population defined\n")
            pop_sums <- FALSE
        }
        if (sum(is.na(pop_sums_t))) {
            if (verbose) cat("\n\t\t\t|--> fill missing --> ")
            pop_sums_t[is.na(pop_sums_t)] <- pop_sums[is.nan(pop_sums_t)]
            pop_sums <- pop_sums_t
            if (verbose) cat("ok ")
        }
    } else {
        pop_sums <- pop_sums_t
    }

    total_pop <- matrix(0, nrow=dim(area_id_dat)[1], ncol=2)
    total_pop <- data.frame(
        code     <- area_id_dat,
        pop      <- pop_sums,
        du       <- du_sums,
        building <- building_sums
        )
    names(total_pop) <- c("area_id", "pop", "du", "building")
    if(verbose) cat("ok ")
    return(total_pop)
}


computeTotal <- function(Tx, benchmark, verbose=FALSE, name="None"){
    if (verbose) {
        cat("\n\t\t\t|--> got benchmarks for: --> ", name, " at: ",
            benchmark[1], benchmark[2], "\t<",
            names(Tx)[benchmark[1]], " ... ",
            names(Tx)[benchmark[2]], ">\t"
                    )}
    benchmark <- as.numeric(benchmark)
    sums <- rowSums(Tx[names(Tx)[benchmark[1]:benchmark[2]]], na.rm=T)
    if (verbose) cat(sum(sums, na.rm=T), " ok ")
    return(sums)
}


filterData <- function(
               census, census_categories, census_area_id,
               survey, survey_categories, survey_id,
               survey_weights, na.rm,
               verbose=verbose){
    if(verbose) cat("\n\t|--> filter data -->")
    # Tx categories
    if(class(census_categories) %in% c("numeric", "integer")){
        if(verbose) cat("numeric categories position -->")
        Tx <- census[census_categories]
    }else if(class(census_categories)=="logical"){
        if(verbose) cat("using all categories -->")
        Tx <- census
    }else{
        if(verbose) cat("named categories position -->")
        Tx <- census[names(census) %in% census_categories]}
    if(verbose) cat(" dim(Tx):", dim(Tx))

    ### Remove area_id
    if(class(census_categories)=="logical"){
        if(class(census_area_id) %in% c("numeric", "integer")){
            if(verbose) cat("\n\t\tnumeric area ID")
            Tx <- Tx[-census_area_id]
        }else{
            if(verbose) cat("\n\t\tnamed area ID")
            Tx <- Tx[!(names(Tx) %in% census_area_id)]
        }
        if(verbose) cat(" dim(Tx):", dim(Tx))
    }

    ### Format data survey
    # dx weights
    if(class(survey_weights) %in% c("numeric", "integer")){
        dx <- survey[survey_weights]
    }else if(class(survey_weights)=="logical"){
        dx <- survey[dim(survey)[2]]
    }else{
        dx <- survey[names(survey) %in% survey_weights]}
    # Transform weights to vector
    weight.colnames <- names(dx)
    dx <- as.numeric(as.matrix(dx))

    # survey id
    if(class(survey_id) %in% c("numeric", "integer")){
        survey_id <- survey[survey_id]
    }else if(class(survey_id) == "logical"){
        survey_id <- as.data.frame(seq(dim(survey)[1]))
    }else{
        survey_id <- survey[names(survey) %in% survey_id]}
    # rename id
    names(survey_id) <- "survey_id"

    # survey categories
    if(verbose) cat("\n\t\tdim(X) --> ", dim(survey))
    if(class(survey_categories) %in% c("numeric", "integer")){
        if(verbose) cat("\n\t\t|--> using categories by index\n\t\t")
        X <- survey[survey_categories]
    }else if(class(survey_categories)=="logical"){
        if(verbose) cat("\n\t\t|--> using all categories\n\t\t")
        X <- survey[!(names(survey) %in% c(weight.colnames, "survey_id"))]
    }else{
        if(verbose) cat("\n\t\t|--> using all categories by name\n\t\t")
        X <- survey[names(survey) %in% survey_categories]}

    # Remove nan
    if(na.rm){
        remove = rowSums(is.na(X)) != 1
        X <- X[remove, ]
        dx <- dx[remove]
        survey_id <- survey_id[remove]}
    if(verbose) cat("dim(Tx):", dim(Tx))
    if(verbose) cat(" length(dx):", length(dx))
    if(verbose) cat(" dim(X):", dim(X))
    if(verbose) cat(" ok")
    list(Tx=Tx, dx=dx, survey_id=survey_id, X=X)
}


#' @title toBinary
#'
#' @description
#' Transforms categorical data to binary rows, if the data type of a column is
#' numerical the function will ignore this column and maintain it as numerical.
#'
#' The function \code{\link{prepareData}} will use this function to transform
#' the inputed survey data.
#'
#' @param X data matrix containing categorical data
#' @param reference_col (optional, default=FALSE) variable defining which
#' columns to use as reference columns. If FALSE the algorithm will
#' automatically select which categories to use as reference, if TRUE all
#' categories will be maintain.
#' @param use_base (optional, default=TRUE) use base R function
#' @param verbose (optional, default=FALSE) be verbose
#' @return X data matrix containing binary data
#' @author M. Esteban Munoz H.
#TODO: examples
toBinary <- function(X, reference_col=FALSE, use_base=TRUE, verbose=FALSE){
    if(verbose) cat("\n\t\t|--> to binary -->")
    # Check if variables are categorical data
    var.type <- sapply(X, class)
    X_numeric <- X[var.type == "numeric"]
    X <- X[var.type == "factor"]

    if(is.logical(reference_col)){
        if(reference_col){
            if(verbose) cat(" (a) Using all columns")
            X_out <- manualModelMatrix(X, type="all",
                                       use_base=use_base, verbose=verbose)
            if(verbose) cat(" with", dim(X_out)[2]-1, "categories- ")
        }else{
            if(verbose) cat(" (b) Automatic selection of reference col")
            X_out <- manualModelMatrix(X, type="normal",
                                       use_base=use_base, verbose=verbose)
        }
    }else{
        if(verbose) cat(" (c) Using reference columns <")
        if(verbose) cat(reference_col)
        X_out <- manualModelMatrix(X, type="all",
                                     use_base=use_base, verbose=verbose)
    }

    for(cn in colnames(X)){colnames(X_out) <- gsub(cn,"",colnames(X_out))}
    # Convert matrix to data frame
    X_out <- as.data.frame(X_out)
    # We don't need the intercept vector
    X_out <- X_out[!(names(X_out) %in% "(Intercept)")]

    # Delete the reference columns
    if(is.list(reference_col) | is.character(reference_col)){
        if(verbose) cat("--> ")
        if(verbose) cat(length(reference_col))
        X_out <- X_out[!(names(X_out) %in% reference_col)]}

    #TODO: I can't use model.matrix and numerical categories
    #X_numeric_bind <- X_numeric[rownames(X_numeric) %in% rownames(X_out), ]
    X_numeric_bind <- X_numeric
    X_out_bind <- cbind(X_out, X_numeric_bind)
    colnames(X_out_bind) <- c(colnames(X_out), colnames(X_numeric))
    if(verbose) cat(" ok")

    return(X_out_bind)
}


#' @title manualModelMatrix
#'
#' @description
#' function to manually create a design matrix
#'
#' @param X input matrix
#' @param type (optional, default="all"). Maintain all categories or create
#' reference categories
#' @param use_base (optional, default=TRUE). Use the R base function
#' `model.matrix`
#' @param verbose (optional, default=FALSE) be verbose
#' @return X_out. Design matrix
#' @author M. Esteban Munoz H.
#TODO: examples
manualModelMatrix <- function(X, type="all",
                              use_base=TRUE, verbose=verbose){
    if(verbose) cat("\n\t\t\t|--> model matrix -->")
    if(use_base){
        if(type=="normal"){
            if(verbose) cat(" normal ")
            X_out <- model.matrix(~ ., data=X)
        }else if(type=="all"){
            if(verbose) cat(" all ")
            X_out <- model.matrix(~ ., data=X,
                contrasts.arg = lapply(X, contrasts, contrasts=FALSE))}
    }else{
        if(verbose) cat(" dummy ")
        X_out <- makeDummy(X, verbose=verbose)}
    if(verbose) cat(" ok")
    return(X_out)}


#' @title makeDummy
#'
#' @description
#' Make dummy variables out of categorical variables.
#'
#' @param X data frame of categorical variables.
#' @param verbose (optional, default=FALSE) be verbose
#' @return X_out data frame of dummy variables
#' @examples
#' data("GREGWT.survey")
#' dummy.survey <- makeDummy(GREGWT.survey)
#' @author M. Esteban Munoz H.
makeDummy <- function(X, verbose=FALSE){

    if(verbose) cat("\n\t\t\t\t|--> make dummy variables -->")
    if(verbose) cat("dim(X):", dim(X))

    lev <- lapply(X, levels)
    lev_names <- unique(unlist(lev))
    lev_num <- sum(sapply(lev,  length))

    temp <- matrix(0, nrow=dim(X)[1], ncol=lev_num)
    X_out <- as.data.frame(temp)
    colnames(X_out) <- lev_names

    for(i in seq(length(lev))){
        X_tmp <- X[names(lev[i])]
        for(l in lev[[i]]){
            X_out[l][X_tmp == l] <- 1}}
    if(verbose) cat("dim(X):", dim(X_out))
    if(verbose) cat(" ok")
    return(X_out)}


namesDiff <- function(d_1, d_2){
    namesd <- names(d_1)[!(names(d_1) %in% names(d_2))]
    return(paste(namesd, collapse="\t"))}


#' @title plotPrepareData
#'
#' @description
#' Visualize the correlation between categories of all constrains.
#'
#' @param prepareData object
#' @return NULL
#' @author M. Esteban Munoz H.
plotPrepareData <- function(data.in,...){
    # Visualize correlation using a correlogram
    require(corrgram)
    X <- data.in$X
    print(dim(X))
    corrgram(X, lower.panel=panel.shade, upper.panel=panel.pie,
         main="Correlation between categories of the benchmarks.")}


#' @title fillMissing
#'
#' @description
#' fill the missing in a data set based on the distribution of the variable
#'
#' @param Tx, data with missing values (dummy values)
#' @param name_breaks, breaks defining the different variables
#' @param omit (optional, default=FALSE), omit a specific variable
#' @return Tx_out, data with estimated values
#' @author M. Esteban Munoz H.
#TODO: make example fillMissing
fillMissing <- function(
        Tx,          # marginal totals with missing values
        name_breaks, # position of the categories breaks
        omit=FALSE   # optional to omit a certain category
                        ){
    name_breaks <- c(1, name_breaks, dim(Tx)[2]+1)

    Tx_out <- data.frame(temp=vector(length=dim(Tx)[1]))

    for(i in seq(length(name_breaks)-1)){
        n <- names(Tx)[name_breaks[i]:(name_breaks[i+1]-1)]
        subTx <- subset(Tx, select=n)
        if(i %in% omit){
            Tx_out <- cbind(Tx_out, subTx)
        }else{
            freq_rel <- colSums(subTx, na.rm=T)/sum(subTx, na.rm=T)
            tempVar <- apply(
                subTx[(rowSums(is.na(subTx))>0), ], 1 , fillNA, freq_rel)
            subTx[(rowSums(is.na(subTx))>0), ] <- t(tempVar)
            Tx_out <- cbind(Tx_out, subTx)}}

    Tx_out <- subset(Tx_out, select=-temp)
    return(Tx_out)}


fillNA <- function(subTx_na, freq_rel){
    col_sel <- names(subTx_na)[is.na(subTx_na)]
    marginal_tot <- sum(subTx_na, na.rm=T)
    subTx_na[is.na(subTx_na)] <- round(
        freq_rel[col_sel] * sum(subTx_na, na.rm=T))
    return(subTx_na)}


#' @title alignWithTotals
#'
#' @description
#' align the census values to total population, the sum of all categories from
#' all variables has to be equal to the population total.
#'
#' @param Tx, un-aligned census data
#' @param pop_totals, population totals vector
#' @param verbose (optional, default=FALSE) be verbose
#' @return Tx_align, align census data
#' @author M. Esteban Munoz H.
#TODO: make example alignWithTotals
alignWithTotals <- function(Tx, breaks, pop_totals, pos,
                            name="standard", verbose=FALSE){
    i = breaks[1] + pos
    for(j in seq(length(breaks)-1)){
        if(length(breaks) == 1){
            b <- breaks[1]
        }else{
            b <- breaks[j+1]
        }
        if(verbose) cat("\n\t\tfrom: ", i, "\tto: ", b)
        this_Tx <- Tx[i:b]
        if(verbose) cat("\n\t\tpop_totals: ", sum(pop_totals), "\t")
        if(verbose) cat("sum Tx: ", sum(this_Tx, na.rm=T), "\n")
        if(sum(this_Tx, na.rm=T) != sum(pop_totals)){
            if(verbose){
                cat("\n\t\t-----------\n")
                cat("\t\t---", name ,"---\n")
                cat("\t\t-----------\n")
                cat("\t\tbreak ", names(Tx)[i], names(Tx)[b], " num: ", i, b)
                cat("\n\t\tdiff: ", sum(Tx[i:b], na.rm=T) - sum(pop_totals))
            }
            if(verbose) cat("\n\t\t|--> inflate -->")
            for(p in seq(1, dim(Tx)[1])){
                Tx[p, i:b] <- inflate(Tx, pop_totals, p, i, b, verbose=verbose)
            }
            if(verbose) cat(" ok")
        }
        i = b+1
    }
    return(Tx)
}


inflate <- function(Tx, pop_totals, p, i, b, verbose=FALSE){
    Tx_row = Tx[p, i:b]
    Tx_row[is.na(Tx_row)] <- 0
    Tx_prop = colMeans(Tx[,i:b], na.rm=T)
    dif <- pop_totals[p,1] - sum(Tx_row, na.rm=T)
    dif_start <- dif
    samp <- sample(seq(i, b), abs(dif), replace=TRUE, prob=Tx_prop)
    k = 0
    for(j in seq(i,b)){
        k = k + 1
        if(dif > 0){
            Tx_row[k] <- Tx[p, j] + length(samp[samp==j])
        }else if(dif < 0){
            Tx_row[k] <- Tx[p, j] - length(samp[samp==j])
        }
    }
    Tx_row[Tx_row<0] <- 0
    dif <- pop_totals[p, 1] - sum(Tx_row, na.rm=T)
    return(Tx_row)
}


#' @title sumFactors
#'
#' @description
#' multiply factors by a given weight
#'
#' @param fact vector containing categorical data
#' @param multiplier weights for each record of fact
#' @return total categories count
#' @examples
#' fact <- as.factor(c(
#'     "age.class.1",
#'     "age.class.2",
#'     "age.class.1",
#'     "age.class.3"))
#' levels(fact)
#' # [1] "age.class.1" "age.class.2" "age.class.3"
#' multiplier <- c(3, 1, 4, 10)
#' sumFactors(fact, multiplier)
#' #   age.class.1 age.class.2 age.class.3
#' # 1           7           1          10
#' @author M. Esteban Munoz H.
sumFactors <- function(fact, multiplier){
    input_variable <- as.factor(as.matrix(fact))
    input_weight <- as.numeric(as.matrix(multiplier))
    output_data <- lapply(levels(input_variable), function(x){
        sum(input_weight[input_variable == x])
                            })
    output_data <- as.data.frame(output_data)
    names(output_data) <- levels(input_variable)
    return(output_data)
}
