# Created by Esteban Munoz (emunozh@gmail.com).
#
# 04.02.2015
# last edit:
# Do 13 Aug 2015 14:17:42 AEST
#

#TODO: document Synthetize
#' @title Synthetize
#'
#' @description
#' Create a synthetic population based on survey weights.
#'
#' @param data_in (default = FALSE) an object of class gregwt, compute with the
#' function gregwt
#' @param X (default = FALSE) original survey, must be define if data_in is not
#' data_in object will overwrite this parameter
#' @param w (default = FALSE) output weights used to sample the survey.
#' data_in object will overwrite this parameter. If data_in = FALSE and w =
#' FALSE all weight will be set to 1.
#' @param benchmarks (default = NULL) use extra benchmarks for the creation of
#' synthetic individuals.
#' @param pop_size_input (default = FALSE) total population size. This
#' parameter will overwrite the parameter on data_in object.
#' @param method (default = c('random', 1)) define the method for creating the
#' synthetic population. (1) c('random', n) will output an array of size
#' [pop,x,n] where pop is the population size, x are the number of parameters on
#' the original survey and n are the number of random samples; (2) c('best', n)
#' will output an array of size [pop,x], this array is the 'best' fit to
#' benchmarks; (3) c('bestpop', n) similar to best bust gives more attention
#' to population totals, important if running an integrated reweight. 
#' @param HHsize_mean_fit (default = 1.8) mean household size.
#' @param max_iter (default = 100) maximal iterations
#' @param fit_tolerance (default = 0.1) tolerance to achieve a fit
#' @param random_seed (default = 12345) seed for random number generator,
#' needed to reproduce results.
#' @param output (default = FALSE) file to attached to the input survey.
#' @param verbose (default = FALSE) be verbose.
#' @return result matrix of synthetic population.
#' @author M. Estebna Munoz H.
#TODO: make example Synthetize
Synthetize <- function(x, ...) UseMethod("Synthetize")

Synthetize.default <- function(
                          data_in = FALSE,
                          X = FALSE,
                          w = FALSE,
                          pop_size_input = FALSE,
                          benchmarks = NULL,
                          method = c("random", 1),
                          group = FALSE,
                          HHsize_mean_fit = 1.8,
                          max_iter = 100,
                          fit_tolerance = 0.1,
                          random_seed = 12345,
                          verbose = FALSE,
                          output = FALSE){

    if (verbose) cat("Initiated\n")
    # set the random seed to ensure reproducibility
    set.seed(random_seed)

    # define number of samples to take from the survey
    n_samples = as.numeric(method[2])

    if (is.logical(data_in)){
        if (is.logical(X)){
            survey <- X
        }
        if (is.logical(w)){
            if (is.logical(pop_size_input)){
                pop_size_input <- dim(X)[1]
                w <- vector(length=dim(X)[1])
                w <- w + 1
            }
        } else {
            if (is.logical(pop_size_input)){
                pop_size_input <- sum(w)
            }
        }
    } else {
        survey <- data_in$X_complete
        w <- data_in$final_weights
    }

    # prepare array for synthetic population
    if (!(is.logical(output))){survey <- cbind(survey, output)}
    survey <- as.data.frame(survey)

    if (!(is.logical(pop_size_input))){
        pop_size_input <- pop_size_input
    }else if (!(is.logical(data_in$pop))){
        pop_size_input <- data_in$pop
    }else{
        stop("missing population total")
    }

    if (is.null(benchmarks)){
        benchmarks <- data_in$Tx_complete
        names(benchmarks) <- data_in$constrains_complete
    }

    # test method
    if ((method[1] == "best" | method[1] == "bestpop") & is.null(benchmarks)){
        stop("benchmarks have to be define to use method 'best' or methoth 'bestpop'")
    }

    if (verbose){
        cat("--Input Data--")
        cat("\n==============================")
        cat("\nrandom seed = ", random_seed)
        cat("\nn samples   = ", n_samples)
        cat("\nsurvey size = ", dim(survey))
        cat("\nweights l   = ", length(w))
        cat("\ntotal pop   = ", pop_size_input)
        cat("\n==============================\n")

    }

    model <- Synthetizeget(survey, w, n_samples, pop_size_input,
                           method = method[1],
                           group = group,
                           HHsize_mean_fit = HHsize_mean_fit,
                           max_iter = max_iter,
                           fit_tolerance = fit_tolerance,
                           random_seed = random_seed,
                           benchmarks = benchmarks,
                           verbose = verbose,
                           output = output)

    weight_index <- unlist(
        lapply(colnames(model), function(x) grepl("weight", tolower(x))))

    if (sum(weight_index) > 0){
        if (verbose) print(weight_index)
        if (verbose) cat("\n dim(mode)", dim(model), "\n")
        model[, weight_index, ] <- 1
    }

    # model <- as.data.frame(model)
    model
}


Synthetizeget <- function(survey, w, n_samples, pop_size_input,
                          method="random",
                          group = FALSE,
                          HHsize_mean_fit = 1.8,
                          max_iter = 100,
                          fit_tolerance = 0.1,
                          random_seed = 12345,
                          benchmarks = NULL,
                          verbose = FALSE,
                          output = FALSE){

    # create the result array
    result=array(NaN, dim=c(pop_size_input, dim(survey)[2], n_samples))
    dimnames(result) <- list(NULL,colnames(survey),NULL)

    if(class(group) == "character"){
        if (verbose) cat("\ngroupped data\n")
        results <- unGroupData(survey, pop_size_input,
                               n_samples, HHsize_mean_fit)
    }else{ if (verbose) cat("\nungroupded data")}

    if (method == "fbs"){
        result <- fbs(n_samples)
    } else {
        for(i in seq(1, n_samples)){
            if (verbose) cat("\n\t|-->sample ", i, "/", n_samples)
            index_survey <- sample(
                nrow(survey), pop_size_input, replace=TRUE, prob=w)
            synthetic_pop <- survey[index_survey, ]
            result[,,i] <- as.matrix(synthetic_pop)
        }
    }

    if (verbose) cat("\nusing method: ", method)

    switch(method,
        "random"={result=result},
        "fbs"={result=result},
        "best"={result=findBest(result, n_samples,
                                benchmarks, verbose=verbose)},
        "bestpop"={
            result=findBest(result, n_samples, benchmarks,
                            pop_size_input, verbose=verbose)
        },
        result=NULL)

    if (verbose) cat("\nOK\n")
    return(result)
}


unGroupData <- function(survey, pop_size_input, n_samples, HHsize_mean_fit){
    # create the result array
    result=array(NaN, dim=c(pop_size_input, dim(survey)[2], n_samples))
    dimnames(result) <- list(NULL,colnames(survey),NULL)

    X.g <- data.frame(HHid=survey[,group],
                    HHsize=(vector(length=dim(survey)[1])+1)) 
    X.g <- aggregate(X.g, by=list(survey[,group]), FUN=sum)
    X.g <- X.g[names(X.g)!=group]

    wx.g <- aggregate(w, by=list(survey[,group]), FUN=mean)
    wx.g <- as.matrix(wx.g)
    wx.g <- as.numeric(wx.g[,"x"])
    mean.w <- mean(wx.g)

    pop_size.sel = pop_size_input / HHsize_mean_fit

    for(i in seq(1,n_samples)){
        HHsize_mean_sel = 0 
        HHsize_sum_sel = 0
        HHsize_delta = Inf
        Popsize_delta = Inf
        iter.num = 0
        while(!(HHsize_mean_fit < HHsize_mean_sel+fit_tolerance  & 
                HHsize_mean_fit > HHsize_mean_sel-fit_tolerance) | 
                pop_size_input != HHsize_sum_sel
                ){
            if(iter.num == max_iter) break
            iter.num = iter.num + 1
            pop.index <- sample(
                nrow(X.g), pop_size.sel, replace=TRUE, prob=wx.g)
            pop_sel_temp <- X.g[pop.index, ]
            HHsize_mean_sel.temp <- mean(pop_sel_temp$HHsize)
            #cat("household size: ", HHsize_mean_sel.temp, "\n")
            HHsize_sum_sel.temp <- sum(pop_sel_temp$HHsize)
            if((abs(HHsize_mean_sel.temp-HHsize_mean_fit)<HHsize_delta) &&
            (abs(HHsize_sum_sel.temp- pop_size_input) <Popsize_delta)){
                HHsize_mean_sel <- HHsize_mean_sel.temp 
                HHsize_delta <- abs(HHsize_mean_sel-HHsize_mean_fit)
                HHsize_sum_sel <- HHsize_sum_sel_temp
                Popsize_delta <- abs(HHsize_sum_sel-pop_size_input)
                pop.sel <- pop_sel_temp
            }else{
                #HHsize_delta.2 <- HHsize_mean_sel-HHsize_mean_fit
                if(HHsize_mean_sel > HHsize_mean_fit){
                    this_index = X.g$HHsize <= HHsize_mean_fit
                }else{
                    this_index = X.g$HHsize >= HHsize_mean_fit}
                wx.g[this_index] <- wx.g[this_index] + mean_w}}

        synthetic_pop <- data.frame(
            matrix(NA, nrow=0, ncol=dim(survey)[2]))
        names(synthetic_pop) <- names(survey)

        p <- pop_sel$Group.1 
        for(j in seq(length(p))){
            spp <- survey[which(X[,group] == p[j]), ]
            spp[group] <- j
            synthetic_pop <- rbind(synthetic_pop, spp)}

        if(dim(result)[1] != dim(synthetic_pop)[1]){
            #cat("people!\t got:", dim(synthetic.pop)[1])
            #cat(" want:", dim(result)[1],"\n")
            synthetic_pop <- fillPop(synthetic_pop, dim(result)[1])}

        result[,,i] <- as.matrix(synthetic_pop)
    }
}


fillPop <- function(synthetic.pop, expected.dim){
    gap <- expected.dim - dim(synthetic.pop)[1]
    #cat("I have a gap of:", gap)
    if(gap>=1){
        #cat(" adding",gap,"individuals\n")
        fill <- vector(length=(dim(synthetic.pop)[2]))
        fill[] <- NaN
        for(i in seq(1,abs(gap))){
            synthetic.pop <- rbind(synthetic.pop, fill)}
    }else{
        to.index = dim(synthetic.pop)[1]+gap
        #cat(" to index:", to.index,"\n")
        synthetic.pop <- synthetic.pop[1:to.index,]}
    return(synthetic.pop)
}


findBest <- function(result, n_samples, benchmarks,
                     pop_size=FALSE, verbose=FALSE){
    # select only columns to benchmak to
    bench_names <- names(benchmarks)
    result_name <- dimnames(result)[[2]]
    index <- sapply(result_name,"%in%",bench_names)
    c_sums <- result[,index,]

    # get the sample marginal totals for the computation of TAE
    if(dim(result)[1]==1){
        # If there is a single person in the sample
        TAE_sample <- sum(colSums(c_sums, na.rm=TRUE))
    }else{
        TAE_sample <- colSums(colSums(c_sums, na.rm=TRUE))
    }

    # get the marginal totals from the benchmarks for the computation of the
    # TAE
    TAE_benchm <- sum(benchmarks)
    # compute the TAE
    TAE_m <- abs(TAE_sample - TAE_benchm)
    # get the min TAE index
    TAE_index <- which(TAE_m == min(TAE_m))

    if(pop_size){
        # Resulting sample, count people
        if(dim(result)[1]==1){
            # If there is a single person in the sample
            rs = sum(!is.nan(result[,1,1]))
        }else{
            rs <- rowSums(!is.nan(aperm(result[,1,])))}
        # population difference
        diff_pop <- abs(rs - pop_size) 
        # select the best sample based on TAE difference and total population
        # difference
        diff_both = TAE_m + diff_pop
        result_index <- which(diff_both == min(diff_both))
    }else{
        result_index <- TAE_index
    }
    return(result[,,result_index[1]])
}
