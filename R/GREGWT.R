# Created by Esteban Munoz (emunozh@gmail.com).
#
# 13.10.2014
# last edit:
# Fr 24 Jul 2015 16:41:52 AEST
#


#TODO: document GREGWT
#' @title GREGWT
#'
#' @description
#' Implementation of the GREGWT algorithm in the R language
#'
#' @param data_in (default=FALSE) input data
#' @param area_code (default = 1) position of the area code in the benchmark data
#' @param use_ginv (default = FALSE) use the ginv function of the MASS library
#' to compute the inverse matrix.
#' @param Tx (default = FALSE) manually specify benchmark vector
#' @param X_input (default = FALSE) manually specify survey matrix as dummy
#' variables without reference categories.
#' @param X_complete_input (default = FALSE) manually specify complete survey
#' matrix as dummy variables with all categories.
#' @param dx_input (default = FALSE) manually specify initial weights.
#' @param area_pop (default = FALSE) manually specify the total population 
#' @param group (default = FALSE) position of the grouping variable, implies an
#' integrated reweighting. 
#' @param bounds (default = c(0,Inf)) defines the bounds of the estimated
#' new weights.
#' @param epsilon (default = 0.001) defines the desire precision of the model
#' to achieve convergence. 
#' @param max_iter (default = 10)
#' @param cat_weights (default = FALSE) weights for individual benchmarks. If
#' not FALSE the function will compute individual weights for each benchmark
#' (or benchmark group) and compute a weighted mean with the weights define in
#' variable 'cat <- weights'.
#' @param align_pop (default = TRUE) align weights to total population.
#' @param error (default = TRUE) compute the error in the simulation.
#' @param verbose (optional, default = 10) be verbose
#' @return GREGWT
#' @examples
#'
#' library('GREGWT')
#' data("GREGWT.census")
#' data("GREGWT.survey")
#' 
#' Simulation.Data <- prepareData(
#'   GREGWT.census, GREGWT.survey,
#'   survey_id=FALSE,
#'   pop_benchmark=c(1,11),
#'   census_categories=seq(2,24),
#'   survey_categories=seq(1,3)
#' )
#' 
#' # reweight for 4 areas
#' areas <- c("02", "11", "04011", "04012")
#' for(acode in areas){
#'     Weights.GREGWT = GREGWT(
#'       data_in=Simulation.Data,
#'       area_code=acode)
#'     print(mean(Weights.GREGWT$final_weights))
#' }
#' 
#' # reweight and plot the results for a single area
#' Tx = Simulation.Data$Tx[which(Simulation.Data$area_id==acode), ]
#' X = as.data.frame(Simulation.Data$X)
#' X_complete = as.data.frame(Simulation.Data$X_complete)
#' Pop = Simulation.Data$total_pop[which(Simulation.Data$area_id==acode), ]
#' Pop = Pop$pop
#' Weights.GREGWT.02 <- GREGWT(
#'   Tx = Tx,
#'   dx_input = Simulation.Data$dx,
#'   X_input = X,
#'   X_complete_input = X_complete,
#'   area_pop = Pop
#'   )
#' 
#' plot(Weights.GREGWT.02)
#' print(Weights.GREGWT.02)
#' summary(Weights.GREGWT.02)
#'
#' @author M. Estebna Munoz H.
#TODO: make example GREGWT
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
                           verbose=FALSE){

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
        #X_complete <- X_complete[,!(colnames(X_complete) %in% "survey_id")]
        X_temp <- X_input
        dx <- dx_input
        total_pop <- area_pop
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
        if(verbose) cat("\n\tget Tx")
        area_index = which(data_in$area_id==area_code)
        Tx <- data_in$Tx[area_index, ]
        Tx_complete <- data_in$Tx_complete[which(data_in$area_id==area_code),]
    } else if (!(is.logical(Tx))){
        if(verbose) cat("\n\tgiven Tx")
        Tx <- Tx
        Tx_complete <- Tx
    } else {
        stop("Either data_in or Tx have to be define as input")
    }
    if (verbose){
        cat("\n\t X_complete: ", dim(X_complete))
        cat("\n\t Tx: ", length(Tx))
        cat("\n\t dx: ", length(dx))
        cat("\n\t total_pop: ", dim(total_pop))
    }
    #if(verbose) cat("Tx:\n")
    #if(verbose) print(Tx)
    #if(verbose) cat("Tx_complete:\n")
    #if(verbose) print(Tx_complete)
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
            if(verbose) cat("\nunknown total population: ", pop)
        }else{
            pop <- total_pop[which(total_pop$area_id==area_code), ]
            pop <- pop$pop
            if(verbose) cat("\nusing total population by area_id: ",
                            area_code, pop)
        }
    }else{
        pop <- area_pop
        if(verbose) cat("\nusing variable <area_pop>: ", pop)
    }
    if(verbose) cat("\npopulation: ", pop)
    if(verbose) cat(" ... ok")

    # Divide initial weights by area population
    if(verbose) cat("\ndivide initial weights by area population")
    if(verbose) cat(" length(dx)", length(dx), "\n")
    if(verbose) cat(" pop: ", pop, "\n")
    if(!(is.logical(pop))){
        dx <- as.numeric(dx)
        dx <- dx*pop/sum(dx)

        # Add population vector
        X <- cbind(X, vector(length=dim(X)[1]) + 1)
        colnames(X)[length(colnames(X))] <- "pop" 
        X_complete <- cbind(X_complete, vector(length=dim(X)[1]) + 1)
        colnames(X_complete)[length(colnames(X_complete))] <- "pop" 

        # And the corresponding population benchmark to Tx
        Tx[length(Tx)+1] <- pop
        names(Tx)[length(names(Tx))] <- "pop"

        Tx_complete[length(Tx_complete)+1] <- pop
        names(Tx_complete)[length(names(Tx_complete))] <- "pop"
        #if(verbose) cat(" consrtains:\n\t", names(Tx_complete), "\n")
    }
    if(verbose) cat(" length(dx)", length(dx))
    if(verbose) cat(" ... ok")

    # Group initial weights for integrated reweight
    if(verbose) cat(" length(dx): ", length(dx))
    dx_oinput <- dx
    if(is.character(group)){
        if(verbose) cat("\ngroup weights if integrated re-weighting")
        dx <- groupDx(X, dx, group, pop)
    }else{
        dx_oinput=FALSE
    }
    if(verbose) cat(" ... ok")

    if(verbose) cat("\narrange data and data formats")
    constrains_complete <- names(Tx_complete)
    constrains_names <- names(Tx)
    constrains_names <- gsub("G.", "", constrains_names)
    Tx <- as.numeric(Tx)
    Tx_complete <- as.numeric(Tx_complete)
    bounds <- as.numeric(bounds)
    epsilon <- as.numeric(epsilon)
    max_iter <- as.numeric(max_iter)
    if(verbose) cat(" ... ok")

    if(verbose) cat("\nconvert NA to 0")
    Tx[is.na(Tx)] <- 0
    Tx_complete[is.na(Tx_complete)] <- 0

    if(verbose) cat("\nGREGWT...")
    model <- GREGWTest(X, dx, Tx, X_complete, Tx_complete, pop,
                       # TODO: do I need the survey_id in the result
                       # survey_id=survey_id,
                       use_ginv=use_ginv, 
                       group=group,
                       bounds=bounds,
                       epsilon=epsilon,
                       max_iter=max_iter,
                       X_input=X_input,
                       dx_oinput=dx_oinput,
                       area_code=area_code,
                       align_pop=align_pop,
                       verbose=verbose)
    if (verbose) cat(" ... ok")

    if (verbose){
        cat("\nprepare output")
        cat("\n\t|--> dim(X)               = ", dim(X))
        cat("\n\t|--> dim(X_complete)      = ", dim(X_complete))
        cat("\n\t|--> length(Tx)           = ", length(Tx))
        cat("\n\t|--> length(Tx_complete)  = ", length(Tx_complete))
        cat("\n\t|--> length(constrains)   = ", length(constrains_names))
        cat("\n\t|--> length(constrains_c) = ", length(constrains_complete))
        cat("\n\t|--> pop                  = ", pop)
        cat("\n\t|--> dim(survey)          = ", dim(survey))
    }
    model$Tx <- Tx
    model$Tx_complete <- Tx_complete
    model$X <- X
    model$X_complete <- X_complete
    model$constrains <- constrains_names
    model$constrains_complete <- constrains_complete
    model$pop <- pop
    model$survey <- survey
    if(verbose) cat(" ... ok")
    cat("\t")

    if (error){
        model <- ComputeError(model, group, verbose=verbose)}
    class(model) <- "GREGWT"
    model}


GREGWTest <- function(X, dx, Tx, X_complete, Tx_complete, pop,
                      # Optional variables
                      # survey_id="survey_id",
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

    cat("AC: ", area_code)

    if(is.character(group)){
        if (verbose) cat("Expanding group ")
        weights <- expandGroup(X_g, wx, X, group)
        dx_output=dx_oinput
        wx_output=weights
    }else{
        dx_output=dx
        wx_output=wx}

    #wx_output <- cbind(survey_id, wx_output)
    if(align_pop){
        if (verbose) cat("population alignment")
        wx_output <- alignPop(wx_output, pop)
    }

    return(list(input_weights=dx_output, final_weights=wx_output))
} 

#TODO: document alignPop #' @title alignPop
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


ComputeError <- function(model, group, verbose=FALSE){
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
    model$Distance.Weights <- abs(w - d)
    # Chi-squared distance
    model$Distance.Chi2 <- 1/2 * (d * w)^2 / d
    # Total Chi-squared distance
    model$TDistance.Chi2 <- sum(1/2 * (d * w)^2 / d)

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
    print(mean(x$Distance.Weights))
    cat("\nMean Chi-squared Distance [mean(1/2 * (w * d)^2 / d)]:\n")
    print(mean(x$Distance.Chi2))
    cat("\nTotal Chi-squared Distance [sum(1/2 * (w * d)^2 / d)]:\n")
    print(x$TDistance.Chi2)
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
        M.Weight.D = mean(x$Distance.Weights),
        M.Chi2.D = mean(x$Distance.Chi2),
        T.Chi2.D = x$TDistance.Chi2,
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
    hTx <- colSums(
        Weights.GREGWT$final_weights * Weights.GREGWT$X_complete, na.rm=T)
    #Tx <- colSums(Weights.GREGWT$input_weights * x$X, na.rm=T)
    Tx <- Weights.GREGWT$Tx_complete
    #PSAE <- as.numeric(abs(Tx-hTx)/length(Weights.GREGWT$final_weights)*100)
    PSAE <- Weights.GREGWT$PSAE
    barplot(PSAE,
            main="Percentage Error of model constrains",
            ylab="PSAE [(Tx - hTx)/n*100]", xlab="")
    # Z-statistic
    r = hTx/sum(Tx)
    p = Tx/sum(Tx)
    Z <- (r-p)/sqrt(p*(1-p)/sum(Tx))
    barplot(as.numeric(Z),
            #names.arg=x$constrains,
            names.arg=names(Weights.GREGWT$X_complete),
            main="Z-Statistic of model constrains",
            ylab="Z-stat [(r-p)/sqrt(p*(1-p)/sum(Tx))]", xlab="")
    plot(as.numeric(Weights.GREGWT$final_weights), Weights.GREGWT$input_weights, 
         ylab="Input Weights [d]",
         xlab="New Weights [w]",
         main="Initial and estimated new weights")
    abline(0,1,col="red")
    plot(sort(as.numeric(
         Weights.GREGWT$final_weights) - Weights.GREGWT$input_weights),
         ylab="Weight distance [w - d]",
         xlab="",
         main="Weight distance")
    abline(h=0,col="red")
    #ED <- abs(sum(x$input_weights)-sum(x$final_weights))/sum(x$input_weights)
    #EM <- sum(x$input_weights)-sum(x$final_weights)/sum(x$input_weights)
    #plot(ED, EM)
}
