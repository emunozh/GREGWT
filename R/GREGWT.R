# Created by Esteban Munoz (emunozh@gmail.com).
#
# 13.10.2014
# last edit:
# Mi 26 Aug 2015 11:46:05 AEST
#


#TODO: document GREGWT
#TODO: implement weighted benchmarks

#' @title GREGWT
#'
#' @description
#' Implementation of the GREGWT algorithm in the R language
#'
#' @param data_in (default=FALSE) input data
#' @param area_code (default = 1) area code to use for the reweighting, defines
#' the area code name. If area_code is a list the GREGWT function will loop
#' through that list and reweight every area on the list. If area_code is of
#' class logical the GREGWT function will loop through all simulation area in
#' the data set.
#' @param use_ginv (default = FALSE) use the ginv function of the MASS library
#' to compute the inverse matrix.
#' @param Tx (default = FALSE) manually specify benchmark vector
#' @param X_input (default = FALSE) manually specify survey matrix as dummy
#' variables without reference categories.
#' @param X_complete_input (default = FALSE) manually specify complete survey
#' matrix as dummy variables with all categories.
#' @param dx_input (default = FALSE) manually specify initial weights.
#' @param area_pop (default = FALSE) manually specify the total population. The
#' area_pop is a vector of size == dim(Tx)[2]. If Tx is a vector area_pop needs
#' to be define as a single value representing the total population of a single
#' area.
#' @param group (default = FALSE) position of the grouping variable, implies an
#' integrated reweighting. 
#' @param bounds (default = c(0,Inf)) defines the bounds of the estimated
#' new weights.
#' @param epsilon (default = 0.001) defines the desire precision of the model
#' to achieve convergence. 
#' @param max_iter (default = 10) maximum number of iterations to run
#' @param cat_weights (default = FALSE, under development) weights for
#' individual benchmarks. If not FALSE the function will compute individual
#' weights for each benchmark (or benchmark group) and compute a weighted mean
#' with the weights define in variable 'cat_weights'.
#' @param align_pop (default = TRUE) align weights to total population.
#' @param error (default = TRUE) compute the error in the simulation.
#' @param output_log (default = FALSE) create en output log of the simulation.
#' This will create a file called GREGWT.out on the current working directory.
#' This option will suppress the output to the command line. You can process
#' the output by calling the function logtocsv() provided by this package.
#' @param verbose (optional, default = FALSE) be verbose
#' @return GREGWT
#' @examples
#'
#' library('GREGWT')
#' data("GREGWT.census")
#' data("GREGWT.survey")
#' 
#' simulation_data <- prepareData(
#'   GREGWT.census, GREGWT.survey,
#'   survey_id=FALSE,
#'   pop_benchmark=c(2,12),
#'   census_categories=seq(2,24),
#'   survey_categories=seq(1,3)
#' )
#' 
#' # reweight for 4 areas
#' areas <- c("02", "11", "04011", "04012")
#' weights_GREGWT = GREGWT(data_in=simulation_data, area_code=areas)
#' plot(weights_GREGWT)
#' 
#' # reweight and plot the results for a single area
#' acode = "02"
#' weights_GREGWT_02 <- GREGWT(data_in= simulation_data, area_code=acode)
#' 
#' plot(weights_GREGWT_02)
#' print(weights_GREGWT_02)
#' summary(weights_GREGWT_02)
#'
#' @author M. Estebna Munoz H.
GREGWT <- function(x, ...) UseMethod("GREGWT")

GREGWT.default <- function(data_in=FALSE,
                           X_input=FALSE,
                           X_complete_input=FALSE,
                           Tx=FALSE,
                           dx_input=FALSE,
                           area_pop=FALSE,     # area population
                           area_code=1,
                           use_ginv=FALSE,
                           group=FALSE,
                           bounds=c(0,Inf),
                           epsilon=0.001,
                           max_iter=10,
                           benchmark_weights=FALSE,
                           align_pop=TRUE,
                           error=TRUE,
                           output_log=FALSE,
                           verbose=FALSE){

    if (output_log){
        log_con <- file("GREGWT.log")
        sink(log_con, append=TRUE)
        sink(log_con, append=TRUE, type="message")
    }

    if (!(is.logical(data_in))){
        if (class(data_in) != "prepareData"){
            stop("data_in has to be of class 'prepareData'")
        }
    }

    if (verbose) cat("\nstart GREGWT, read variables")
    if (is.logical(data_in) &
       !(is.logical(X_input)) &
       !(is.logical(dx_input))
       ){
        if (verbose) cat("\nusing X_complete, Tx and dx")
        X_complete <- X_complete_input
        X_temp <- X_input
        dx <- dx_input
        if (!(is.null(dim(Tx)))){
            if (is.numeric(area_code)){
                total_pop <- area_pop[area_code]
            } else if (is.character(area_code)){
                total_pop <- area_pop[rownames(Tx) %in% area_code, ]
            }
        } else {
            total_pop <- area_pop
        }
        survey <- X_complete_input
    } else if (!(is.logical(data_in))){
        if (verbose) cat("\nusing data_in")
        X_complete  <- data_in$X_complete
        X_temp <- data_in$X
        dx <- data_in$dx
        total_pop  <- data_in$total_pop
        survey <- data_in$survey
    } else {
        stop("Either data_in or X_input and dx_input have to be define")
    }

    survey_id <- X_temp[colnames(X_temp) %in% "survey_id"]
    X <- X_temp[,!(colnames(X_temp) %in% "survey_id")]
    X <- as.matrix(X)
    if (verbose) cat(" ... ok")

    # Get area benchmarks
    if (verbose) cat("\nget area benchmarks")
    if (is.logical(Tx) & !(is.logical(data_in))){
        if(verbose) cat("\n\tget Tx ")
        if (length(area_code) == 1){
            if(verbose) cat("(for single area):")
            area_index = which(data_in$area_id==area_code)
        } else {
            if(verbose) cat("(for", length(area_code), "areas):")
            area_index = data_in$area_id[,1] %in% area_code
        }
        Tx <- data_in$Tx[area_index, ]
        Tx_complete <- data_in$Tx_complete[area_index,]
    } else if (!(is.logical(Tx))){
        if(verbose) cat("\n\tgiven Tx")
        if (is.numeric(area_code) & !(is.null(dim(Tx)))){
            Tx <- Tx[area_code, ]
        } else if (is.character(area_code) & !(is.null(dim(Tx)))){
            Tx <- Tx[rownames(Tx) %in% area_code, ]
        }
        Tx_complete <- Tx
    } else {
        stop("Either data_in or Tx have to be define as input")
    }
    if (verbose){
        cat("\n\tX_complete: ", dim(X_complete))
        if (!(is.null(dim(Tx)))){
            cat("\n\tTx: ", dim(Tx))
            cat("\n\ttotal_pop: ", dim(total_pop)[1], dim(total_pop)[2])
        } else {
            cat("\n\tTx: ", length(Tx))
            cat("\n\ttotal_pop: ", dim(total_pop)[1], dim(total_pop)[2])
        }
        cat("\n\tdx: ", length(dx))
    }
    if(verbose) cat(" ... ok")

    ## Define breaks and weighted benchmarks
    #if (!(is.logical(data_in)) & is.logical(breaks)){
    #    breaks <- data_in$breaks
    #} else if (is.logical(data_in) & is.logical(breaks)){
    #    stop("Either data_in or breaks have to be defined as input")
    #}
    #breaks <- c(breaks, dim(Tx)[2])
    #if((length(breaks) != length(benchmark_weights)) & !(is.logical(benchmark_weights))){
    #    #TODO: implement weighted benchmarks
    #    stop("Expected ", length(breaks), "weights, ",
    #         length(benchmark_weights), "given\n")
    #}

    # Get total area population
    if(verbose) cat("\nget area population")
    if(is.logical(area_pop)){
        if(is.logical(total_pop)){
            pop <- FALSE
            align_pop <- FALSE
            if(verbose) cat("\nunknown total population: ", pop)
        }else{
            if (length(area_code) == 1){
                pop <- total_pop[which(total_pop$area_id==area_code), ]
                pop <- pop$pop
            } else {
                pop <- total_pop[total_pop$area_id %in% area_code, ]
                pop <- pop$pop
            }
            if(verbose) cat("\nusing total population by area_id: ",
                            area_code, pop)
        }
    }else{
        pop <- area_pop
        if(verbose) cat("\nusing variable <area_pop>: ", pop)
    }
    if(verbose) cat("\npopulation: ", pop)
    if(verbose) cat(" ... ok")

    if(!(is.logical(pop))){
        # Add population vector
        X <- cbind(X, vector(length=dim(X)[1]) + 1)
        colnames(X)[length(colnames(X))] <- "pop" 
        X_complete <- cbind(X_complete, vector(length=dim(X)[1]) + 1)
        colnames(X_complete)[length(colnames(X_complete))] <- "pop" 

        # And the corresponding population benchmark to Tx
        if (is.null(dim(Tx))){
            Tx[length(Tx)+1] <- pop
            Tx_complete[length(Tx_complete)+1] <- pop
        } else {
            Tx <- cbind(Tx, pop)
            Tx_complete <- cbind(Tx_complete, pop)
        }
        names(Tx)[length(names(Tx))] <- "pop"
        names(Tx_complete)[length(names(Tx_complete))] <- "pop"
        
    }
    if(verbose) cat(" ... ok")

    if(verbose) cat("\narrange data and data formats")
    constrains_complete <- names(Tx_complete)
    constrains_names <- names(Tx)
    constrains_names <- gsub("G.", "", constrains_names)
    Tx <- as.matrix(Tx)
    Tx_complete <- as.matrix(Tx_complete)
    bounds <- as.numeric(bounds)
    epsilon <- as.numeric(epsilon)
    max_iter <- as.numeric(max_iter)
    if(verbose) cat(" ... ok")

    if(verbose) cat("\nconvert NA to 0")
    Tx[is.na(Tx)] <- 0
    Tx_complete[is.na(Tx_complete)] <- 0

    # get number of simulation areas
    if (verbose) cat("\ndim(Tx) --> ", dim(Tx), "\n")
    if (is.null(dim(Tx))){
        area_numbers <- 1
    } else {
        if (dim(Tx)[2] == 1){
            if (verbose) cat("\nndim 2 TX", "\n")
            area_numbers <- dim(Tx)[2]
        } else {
            if (verbose) cat("\nndim 1 TX", "\n")
            area_numbers <- dim(Tx)[1]
        }
    }
   
    # loop throgh all simulation areas
    final_weights = as.numeric()
    TAE <- as.numeric()
    SAE <- as.numeric()
    PSAE <- as.numeric()
    TAD <- as.numeric()
    DW <- as.numeric()
    DChi2 <- as.numeric()
    TDChi2 <- as.numeric()
    Z <- as.numeric()
    EM <- as.numeric()
    ED <- as.numeric()
    for (i in seq(area_numbers)){
        if (verbose) cat("\n", "area_numbers -->", area_numbers, "\n")
        if (verbose) cat("\n", "loop --> ", i, "\n")
        if (area_numbers == 1){
            area_code_i <- area_code
            Tx_i <- Tx
            Tx_complete_i <- Tx_complete
            pop_i <- pop
        } else {
            area_code_i <- area_code[i]
            Tx_i <- Tx[i, ]
            Tx_complete_i <- Tx_complete[i, ]
            pop_i <- pop[i]
        }

    # Divide initial weights by area population
    if(!(is.logical(pop_i))){
        if(verbose) cat("\ndivide initial weights by area population")
        if(verbose) cat(" length(dx)", length(dx), "\n")
        if(verbose) cat(" pop: ", pop_i, "\n")
        dx <- as.numeric(dx)
        dx <- dx*pop_i/sum(dx)
    }
    if(verbose) cat(" length(dx)", length(dx))
    if(verbose) cat(" ... ok")

    # Group initial weights for integrated reweight
    dx_oinput <- dx
    if(is.character(group)){
        if(verbose) cat("\ngroup weights if integrated reweighting")
        dx <- groupDx(X, dx, group, pop_i)
        if(verbose) cat(" length(dx-group): ", length(dx))
    }else{
        dx_oinput=FALSE
    }
    if(verbose) cat(" ... ok")

    if(verbose) cat("\nGREGWT...")
    model_iter <- GREGWTest(X, dx, Tx_i, X_complete, Tx_complete, pop_i,
                       survey_id=survey_id,
                       use_ginv=use_ginv, 
                       group=group,
                       bounds=bounds,
                       epsilon=epsilon,
                       max_iter=max_iter,
                       X_input=X_input,
                       dx_oinput=dx_oinput,
                       area_code=area_code_i,
                       align_pop=align_pop,
                       verbose=verbose)
    if (verbose) cat(" ... ok")

    model_iter$X <- X
    model_iter$X_complete <- X_complete
    model_iter$constrains <- constrains_names
    model_iter$constrains_complete <- constrains_complete
    model_iter$survey <- survey
    # loop
    model_iter$Tx <- Tx_i
    model_iter$Tx_complete <- Tx_complete_i
    model_iter$pop <- pop_i

    if (error){
        model_iter <- computeError(model_iter, group, verbose=verbose)
        TAE <- cbind(TAE, model_iter$TAE)
        SAE <- cbind(SAE, model_iter$SAE)
        PSAE <- cbind(PSAE, model_iter$PSAE)
        TAD <- cbind(TAD, model_iter$TAD)
        DW <- cbind(DW, model_iter$DW)
        DChi2 <- cbind(DChi2, model_iter$DChi2)
        TDChi2 <- cbind(TDChi2, model_iter$TDChi2)
        Z <- cbind(Z, model_iter$Z)
        EM <- cbind(EM, model_iter$EM)
        ED <- cbind(ED, model_iter$ED)
    } else {
        cat("\t|\n")
    }

    final_weights <- cbind(final_weights, model_iter$final_weights)
    } # end loop areas
   
    model <- model_iter #TODO: is error implemented

    model$TAE <- rowMeans(TAE)
    model$SAE <- rowMeans(SAE)
    model$PSAE <- rowMeans(PSAE)
    model$TAD <- mean(TAD)
    model$DW <- mean(DW)
    model$DChi2 <- mean(DChi2)
    model$TDChi2 <- mean(TDChi2)
    model$Z <- rowMeans(Z)

    model$final_weights <- final_weights
    model$input_weights <- model_iter$input_weights
    model$Tx <- Tx
    model$Tx_complete <- Tx_complete
    model$pop <- pop
    model$X <- X
    model$X_complete <- X_complete
    model$constrains <- constrains_names
    model$constrains_complete <- constrains_complete
    model$survey <- survey

    if (verbose){
        cat("\nprepare output")
        cat("\n\t|--> dim(X)               = ", dim(X))
        cat("\n\t|--> dim(X_complete)      = ", dim(X_complete))
        if (is.null(dim(Tx))){
            cat("\n\t|--> length(Tx)           = ", length(Tx))
            cat("\n\t|--> length(Tx_complete)  = ", length(Tx_complete))
        } else {
            cat("\n\t|--> dim(Tx)           = ", dim(Tx))
            cat("\n\t|--> dim(Tx_complete)  = ", dim(Tx_complete))
        }
        cat("\n\t|--> length(constrains)   = ", length(constrains_names))
        cat("\n\t|--> length(constrains_c) = ", length(constrains_complete))
        if (is.null(dim(Tx))){
            cat("\n\t|--> pop                  = ", pop)
        } else {
            cat("\n\t|--> dim(pop)             = ", dim(pop))
        }
        cat("\n\t|--> dim(survey)          = ", dim(survey))
    }

    class(model) <- "GREGWT"

    if (output_log){
        sink() 
        sink(type="message")
    }
    return(model)
}


GREGWTest <- function(X, dx, Tx, X_complete, Tx_complete, pop,
                      # Optional variables
                      survey_id="survey_id",
                      use_ginv=FALSE,
                      group=FALSE,
                      bounds=c(-Inf,Inf),
                      epsilon=0.001,
                      max_iter=10,
                      X_input=FALSE,
                      dx_oinput=FALSE,
                      align_pop=TRUE,
                      area_code="/",
                      verbose=FALSE){

    # Save the aggregation ids
    X_g <- X

    if(verbose) cat("\n\tgroup...")
    if(is.character(group)){
        X <- X[,colnames(X)[colnames(X)!="Group.1"]]
        X_complete <- X_complete[,colnames(X_complete)[colnames(X_complete)!="Group.1"]]
    }
    if(verbose) cat("ok")

    # Number of attributes in the input sample
    att_num <- dim(X)[2]

    if(verbose){
        cat("\n\tget first lambda")
        cat("\n\tdim(X)", dim(X))
        #cat("\n\tcolnames(X)", colnames(X))
        cat("\n\tlength(dx)", length(dx))
    }
    # Lambda
    hTx <- colSums(X * dx, na.rm=T)  # Sample totals
    lambda <- getLambda(X, dx, Tx, hTx, att_num, use_ginv, verbose=verbose)

    if(is.na(sum(lambda))){
        cat("Warning, nans in lambda vector")
        return(list(input_weights=F, final_weights=F))
    }

    if(verbose) cat("ok")

    if(verbose) cat("\n\tmain loop")
    # Truncate weights
    convergence = FALSE
    number_iter = 0
    while(!convergence){
        number_iter = number_iter + 1

        # Get new weights
        wx = dx * (1 + X %*% lambda)  # New weights

        # Truncate weights
        wx[wx<bounds[1]] <- bounds[1]
        wx[wx>bounds[2]] <- bounds[2]
        
        # Truncate initial wights 
        dx[wx<bounds[1]] <- 0
        dx[wx>bounds[2]] <- 0
       
        # Recompute lambda
        if(verbose) cat("\n\t\tNew Lambda")
        hTx <- colSums(X * as.numeric(wx), na.rm=T)  # Sample totals
        lambdaS <- getLambda(X, dx, Tx, hTx, att_num, use_ginv, verbose=verbose)
        # Save lambda m-1
        lambdaO <- lambda
        # Compute new lambda
        lambda = lambda + lambdaS

        # Compute values for convergence
        if(all(Tx > epsilon)){
            #dlta_tx <- abs(Tx_complete - hTx_complete)
            dlta_tx <- abs(Tx - hTx)
            dlta_l  <- abs(lambdaS-lambdaO)
        }else{
            dlta_tx <- 0
            dlta_l <- 0
        }

        convergence <- (
            all(dlta_tx < epsilon) |
            all(dlta_l <= epsilon) |
            number_iter >= max_iter)

        if(convergence){
            cat("| dTx | ",
                format(max(dlta_tx), digits=2, scientific=T, width=7), "< ",
                format(epsilon, digits=2, scientific=T) ,"-> ")
            cat(format(all(dlta_tx < epsilon), width=5), " ")
            cat("| dAm | ",
                format(max(dlta_l), digits=2, scientific=T, width=7), "<=",
                format(epsilon, digits=2, scientific=T) ,"-> ")
            cat(format(all(dlta_l <= epsilon), width=5), " ") 
            cat("| itr: ", format(number_iter, digits=0, width=4), "")}
        }
    if(verbose) cat("ok")
    if(number_iter >= max_iter){
        cat("| NC | ")
    }else{
        cat("| OK | ")}

    cat("AC: ", format(area_code, scientific=FALSE))

    if(is.character(group)){
        if (verbose) cat("Expanding group ")
        weights <- expandGroup(X_g, wx, X, group)
        dx_output=dx_oinput
        wx_output=weights
    }else{
        dx_output=dx
        wx_output=wx}

    wx_output <- cbind(survey_id, wx_output)
    if (verbose) cat("bind index")
    if(align_pop){
        if (verbose) cat("population alignment")
        wx_output <- alignPop(wx_output, pop)
    }

    return(list(input_weights=dx_output, final_weights=wx_output))
} 

#TODO: document alignPop 

#' @title alignPop
#'
#' @description #' Aligns the survey weights with the total population
#'
#' @param w new computed weights
#' @param p population total
#' @return wo align weights
#' @author M. Estebna Munoz H.  
#TODO: make example alignPop
alignPop <- function(w, p){
    #wo <- w + w / sum(w) * (p - sum(w))
    wo <- p * w / sum(w)
    return(wo)
} 


getLambda <- function(X, dx, Tx, hTx, att_num, use_ginv, verbose=FALSE){
    if(verbose) cat("\n\t\tget A")
    A <- crossprod(as.matrix(dx * X), as.matrix(X))
    if(verbose){
        cat("... ok")
        if(sum(is.na(A))) cat("!!A contains NaN's!!")
        if(sum(is.na(Tx))) cat("!!Tx contains NaN's!!")
        if(sum(is.na(hTx))) cat("!!hTx contains NaN's!!")
    }
    # A solution for 'A %*% lambda = (Tx - hTx)'
    if(use_ginv){
        if(verbose) cat("\n\t\t\tusing ginv ")
        if(verbose) cat("\n\t\tlength(Tx): ", length(Tx))
        if(verbose) cat("\t\tlength(hTx): ", length(hTx))
        require('MASS')
        At <- ginv(A)
        lambda = At %*% (Tx - hTx)
        if(verbose) cat("ok")
    }else{
        if(verbose) cat("\n\t\t\tfind a solution for A ")
        #A <- solve(A)
        #lambda = A %*% (Tx - hTx)
        lambda <- solve(A, (Tx - hTx))   # Lagrange multipliers
        if(verbose) cat("ok")
    }
    if(verbose) cat("ok")
    if(verbose){
        cat("\n\t\tLambda:")
        cat(lambda)
    }
    return(lambda)}


groupX <- function(X, group){
    X_g <- aggregate(X, by=list(X[,group]), FUN=sum)
    X_g <- X_g[names(X_g)!=group]
    #X_g <- as.matrix(X_g)
    return(X_g)}


#TODO: divide initial weights by area population
groupDx <- function(X_input, dx, group, pop){
    dx_g <- aggregate(dx, by=list(X_input[,group]), FUN=sum)
    dx_g <- as.matrix(dx_g)
    dx_g <- as.numeric(dx_g[,"x"])
    return(dx_g)}


#TODO: make example expandGroup
expandGroup <- function(X_g, wx, X_input, group){
    Origin <- X_g[,"Group.1"]
    Grouped <- X_input[,group]
    index <- match(Grouped, Origin, nomatch=0) 
    weights <- wx[index]
    return(weights)}


computeError <- function(model, group, verbose=FALSE){
    if(verbose) cat("\n compute error")
    model$call <- match.call()
    # Weights
    d <- model$input_weights
    w <- model$final_weights
    # Marginal sums 
    Tx <- model$Tx_complete
    #cat("\n\n", dim(model$X_complete))
    #cat("\n\n", length(w))
    #print(head(model$X_complete * w))
    hTx <- colSums(model$X_complete * w)[model$constrains_complete]
    #cat("\nsum(X_c) = ", sum(model$X_complete))
    #cat("\nsum(X) = ", sum(model$X))
    #cat("\nsum(Tx) = ", sum(Tx))
    #cat("\nsum(hTx) = ", sum(hTx))
    #cat("\nsum(w) = ", sum(w), "\n")

    #if(verbose){
    #    cat("\n Tx:")
    #    print(class(Tx))
    #}
    #Tx <- Tx[, order(colnames(Tx))]

    # Compute the internal error of the computation
    # Total absolute error (TAE)
    model$TAE <- abs(Tx-hTx)
    cat(" | TAE:", format(sum(model$TAE),digits=2,scientific=T,width=7))
    # Standardized absolute error (SAE)
    model$SAE <- abs(Tx-hTx) / model$pop
    # Percentage error (PSAE)
    if(verbose) cat("\n using total population: ", model$pop)
    model$PSAE <- model$SAE * 100
    cat(" | PSAE:", format(sum(model$PSAE),digits=2,scientific=T,width=7))
    ## Correlation Coefficient (Pearson Correlation)
    #model$pearson <- cor(cbind(Tx, hTx), use="complete.obs", method="pearson")
    ## Independent samples t-Test
    #model$ttest <- t.test(Tx, hTx)
    ## Coefficient of determination
    #lm.X <- lm(Tx~hTx)
    #model$r2 <- summary(lm.X)$r.squared
    #model$r2.adj <- summary(lm.X)$adj.r.squared

    # Total absolute distance (TAD)
    model$TAD <- sum(abs(w-d))
    # Weights Distance
    model$DW <- abs(w - d)
    # Chi-squared distance
    model$DChi2 <- 1/2 * (d * w)^2 / d
    # Total Chi-squared distance
    model$TDChi2 <- sum(1/2 * (d * w)^2 / d)

    # Z-statistic
    r = hTx/sum(Tx, na.rm=T)
    p = Tx/sum(Tx, na.rm=T)
    model$Z <- (r-p)/sqrt(p*(1-p)/sum(Tx, na.rm=T))
    # Error in Margin (EM)
    model$EM <- (sum(d) - sum(w)) / sum(d)
    # Error in Distribution (ED)
    model$ED <- abs(sum(d) - sum(w)) / sum(d)

    cat(" |\n")

    return(model)}


print.GREGWT <- function(x, ...){
    #cat("Call:\n")
    #print(x$call)
    cat("\nMean Weight Distance [mean(abs(w - d))]:\n")
    print(mean(x$DW))
    cat("\nMean Chi-squared Distance [mean(1/2 * (w * d)^2 / d)]:\n")
    print(mean(x$DChi2))
    cat("\nTotal Chi-squared Distance [sum(1/2 * (w * d)^2 / d)]:\n")
    print(x$TDChi2)
    cat("\nTotal absolute distance [sum(abs(w - d))]:\n")
    print(x$TAD)
    cat("\nTotal absolute error [abs(Tx - hTx)]:\n")
    print(x$TAE)
    cat("\nStandardized absolute error [TAE / n]:\n")
    print(x$SAE)
    cat("\nPercentage absolute error [SAE * 100]:\n")
    print(x$PSAE)
    #cat(Correlation Coefficient (Pearson Correlation))
    #print(x$pearson)
    #cat("\nIndependent samples t-Test [t.test(Tx, hTx)]:\n")
    #print(x$ttest)
    #cat("\nR squared (Tx~hTx):\n")
    #print(x$r2)
    #cat("\nAdjusted R squared:\n")
    #print(x$r2.adj)
    cat("\nModified z-statistic [Z = (r-p)/sqrt(p*(1-p)/sum(Tx))]:\n")
    print(x$Z)
    cat("\nError in margin (EM) [(sum(d)-sum(w))/sum(d)]:\n")
    print(x$EM)
    cat("\nError in distribution (ED) [abs(sum(d)-sum(w))/sum(d)]:\n")
    print(x$ED)
    cat("\nNew Weights. Access via: X$final_weights:\n")
    print(paste(length(x$final_weights),"New weights"))}


summary.GREGWT <- function(x, ...){
    ttest = x$ttest
    res <- list(
        call=x$call,
        # 1 ######
        TAD = sum(x$TAD),
        M.Weight.D = mean(x$DW),
        M.Chi2.D = mean(x$DChi2),
        T.Chi2.D = x$TDChi2,
        # 2 ######
        TAE = sum(x$TAE),
        SAE = sum(x$SAE),
        PSAE = mean(x$PSAE),
        # 3 ######
        #p.val.ttest = ttest$p.value,
        #pearson = mean(x$pearson[x$pearson==1]),
        #r2 = x$r2,
        #r2.adj = x$r2.adj,
        # 4 ######
        Z = mean(x$Z), 
        EM = x$EM,
        ED = x$ED
        )
    class(res) <- "summary.GREGWT"
    res}


print.summary.GREGWT <- function(x, ...){
    TWidth = 10
    #cat("Call:\n")
    #print(x$call)
    cat("\n================================================\n")
    cat(format("TAD", width=TWidth), "|")
    cat(format("M.Weight.D", width=TWidth), "|")
    cat(format("M.Chi2.D", width=TWidth), "|")
    cat(format("T.Chi2.D", width=TWidth), "|")
    cat("\n------------------------------------------------\n")
    cat(format(x$TAD, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$M.Weight.D, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$M.Chi2.D, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$T.Chi2.D, digits=2, scientific=T, width=TWidth), "|")
    cat("\n------------------------------------------------\n")
    cat(format("TAE", width=TWidth), "|")
    cat(format("SAE", width=TWidth), "|")
    cat(format("PSAE", width=TWidth), "|")
    cat(format("", width=TWidth), "|")
    cat("\n------------------------------------------------\n")
    cat(format(x$TAE, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$SAE, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$PSAE, digits=2, scientific=T, width=TWidth), "|")
    cat(format("", width=TWidth), "|")
    cat("\n------------------------------------------------\n")
    #cat(format("t-Test", width=TWidth), "|")
    #cat(format("M.Pearson", width=TWidth), "|")
    #cat(format("R2", width=TWidth), "|")
    #cat(format("adj R2", width=TWidth), "|")
    #cat("\n------------------------------------------------\n")
    #cat(format(x$p.val.ttest, digits=2, scientific=T, width=TWidth), "|")
    #cat(format(x$pearson, digits=2, scientific=T, width=TWidth), "|")
    #cat(format(x$r2, digits=2, scientific=T, width=TWidth), "|")
    #cat(format(x$r2.adj, digits=2, scientific=T, width=TWidth), "|")
    #cat("\n------------------------------------------------\n")
    cat(format("M.Z", width=TWidth), "|")
    cat(format("EM", width=TWidth), "|")
    cat(format("ED", width=TWidth), "|")
    cat(format("", width=TWidth), "|")
    cat("\n------------------------------------------------\n")
    cat(format(x$Z, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$EM, digits=2, scientific=T, width=TWidth), "|")
    cat(format(x$ED, digits=2, scientific=T, width=TWidth), "|")
    cat(format("", width=TWidth), "|")
    cat("\n================================================\n")
    #printCoefmat(x$Distance)
}


plot.GREGWT <- function(x, ...){
    layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
    names_y <- x$constrains_complete
    names_y <- gsub("G.", "", names_y)
    hTx <- colSums(
        x$final_weights * x$X_complete, na.rm=T)
    Tx <- x$Tx_complete
    #PSAE <- as.numeric(abs(Tx-hTx)/length(x$final_weights)*100)
    #PSAE <- x$PSAE
    PSAE <- abs(Tx-hTx) / x$pop * 100
    barplot(PSAE,
            names.arg=names_y,
            main="Percentage Error of model constrains",
            ylab="PSAE = (Tx - hTx)/n*100", xlab="")
    # Z-statistic
    r = hTx/sum(Tx, na.rm=T)
    p = Tx/sum(Tx, na.rm=T)
    Z <- (r-p)/sqrt(p*(1-p)/sum(Tx, na.rm=T))
    barplot(Z,
            names.arg=names_y,
            main="Z-Statistic of model constrains",
            ylab="Z = (r-p)/sqrt(p*(1-p)/sum(Tx))", xlab="")

    if (dim(x$final_weights)[2] == 1){
        area_numbers <- 1
        plot(as.numeric(x$final_weights), x$input_weights, 
             pch=".",
             ylab="Input Weights [d]",
             xlab="New Weights [w]",
             main="Initial and estimated new weights")
    } else {
        area_numbers <- dim(x$final_weights)[2]
        plot(as.numeric(x$final_weights[,1]), x$input_weights, 
             pch=".",
             ylab="Input Weights [d]",
             xlab="New Weights [w]",
             main="Initial and estimated new weights")
        for (i in seq(2, area_numbers)){
            points(as.numeric(x$final_weights[, i]),
                x$input_weights, col=i, pch=".")
    }}
    abline(0,1,col="red")

    plot(sort(as.numeric(
         x$final_weights) - x$input_weights),
         ylab="Weight distance [w - d]",
         xlab="",
         main="Weight distance")
    abline(h=0,col="red")
    #ED <- abs(sum(x$input_weights)-sum(x$final_weights))/sum(x$input_weights)
    #EM <- sum(x$input_weights)-sum(x$final_weights)/sum(x$input_weights)
    #plot(ED, EM)
}

#' @title logtocsv
#'
#' @description
#' This function will try to read the log output created by function GREGWT
#' with the output_log variable set to TRUE and generate a csv file on the same
#' path called GREGWT_log.csv.
#'
#' @param file_name (default = FALSE) alternative file name to process, if
#' FALSE the function will look for a file name GREGWT.log in the root folder.
#'
#' @author M. Estebna Munoz H.
logtocsv <- function(file_name=FALSE){
    file_name <- "GREGWT.log" 
    data_table <- read.table(file_name)
    if (dim(data_table)[2] == 25){
        col_index <- c(24,8,16,19,21)
        col_names <- c("area_code", "dTx_convergence", "dAm_convergence",
                       "iterations", "convergence")
    } else {
        col_index <- c(24,8,16,19,21,27,30)
        col_names <- c("area_code", "dTx_convergence", "dAm_convergence",
                       "iterations", "convergence", "TAE", "PSAE")
    }
    data_table <- data_table[col_index]
    names(data_table) <- col_names
    data_table$area_code <- as.character(data_table$area_code)
    data_table$convergence <- as.character(data_table$convergence)
    data_table$convergence[data_table$convergence == "OK"] <- "TRUE"
    data_table$convergence[data_table$convergence == "NC"] <- "FALSE"
    data_table$convergence <- as.logical(data_table$convergence)
    write.csv(data_table, file="GREGWT_log.csv", row.names=FALSE)
}
