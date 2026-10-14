
#' @name cellMeans
#'
#' @title Cell Means
#'
#' @description Extract cell means and other information from a linear model. 
#' 
#' @param object \code{\link[stats]{lm}} or \code{\link[stats]{aov}} model object
#' @param factors character vector providing the model terms (factors/categorical predictors) to obtain the means for each level or level combinations
#' @param var.equal currently ignored
#' @param ... currently ignored
#'
#' @details The purpose of this function is to summarize results from a linear model in preperation to test contrast(s) with a call to the \code{\link[linearModel]{contrast}} function.
#'
#' The weights of the linear model object are currently being ignored. 
#'
#' 
#' @return Object of class anova and data.frame
#'
#' @seealso \code{\link[linearModel]{contrast}}
#'
#' @export
#'
#' @import stats
#' @import utils
#' @import methods
#' 
#'
#' @examples
#' data(genericData)
#'
#' mod <- lm(Y~A+B+A*B,data=genericData)
#' cellMeans(mod, factors='A')
#' cellMeans(mod)
#' 




cellMeans <- function(object, factors=NULL, var.equal=TRUE, ...){


    ## ########################################
    ## argument checking
    ## ########################################

    if(missing(var.equal) || !is.vector(var.equal) || !is.logical(var.equal) || length(var.equal) !=1){
        stop('Argument var.equal must be a single logical value')
    }#end

    
    if(missing(object) || !methods::is(object, 'lm')){
        stop('Argument object must be a fitted model from a call to lm() or aov()')
    }#end if

    m <- object
    ## data from model
    df <- m$model
    (responseName <- names(df)[attr(m$terms,'response')])
    (allPredNames <- names(df)[names(df)%in%attr(m$terms, 'term.labels')])


    if(!is.null(factors)){

        if(!is.vector(factors) | !is.character(factors)){

            stop('Argument factors must be a character vector')
        }#end if


        predInModel <- factors%in%allPredNames


        if(all(!predInModel)){

            stop('Argument factors must be a character vector \n',
                 'These values must match the terms in the model argument')
        }#end if
        
        if(any(!predInModel)){
            warning('The following terms do not exist in the model and will be ignored: \n',
            paste0(factors[!predInModel], collapse=', '))
        }#end if


        predNames <- factors[predInModel]

                
    }else{

        predNames <- allPredNames

    }#end else if
    
    isFactor <- (sapply(df[,predNames], is.factor) | sapply(df[,predNames], is.character))

    if(all(!isFactor)){
        stop('All terms in the model are not factors/characters, cell means cannot be calculated')
    }#end if


    if(any(!isFactor)){


        warning('The follow terms in the model are not factors/characters and will be ignored: \n',
                paste0(predNames[!isFactor], collapse=', '))

                predNames <- predNames[isFactor]
    }#end if
    

    checkNames <- predNames%in%c('mean','n','sd')
    if(any(checkNames)){
        names(df)[names(df)%in%predNames[checkNames]] <- base::toupper(names(df)[names(df)%in%predNames[checkNames]])
        predNames[checkNames] <- toupper(predNames[checkNames])
    }#end if 


    

    ## ########################################
    ## calculate means and extract info from the model 
    ## ########################################

    

    aggFun <- function(x, rn){
        out <- unique(x[,names(x)!=rn,drop=FALSE])
        out$mean <- mean(x[,rn])
        out$sd <- sd(x[,rn])
        out$n <- length(x[,rn])
        return(out)
    }#end aggFun

    outMeans <- do.call(rbind,lapply(split(x=df[,c(responseName,predNames)],f=df[,predNames]), FUN=aggFun, rn=responseName))
    row.names(outMeans) <- NULL

    outMeans

    ## extract from lm object
    AT <- linearModel::anovaTable(m)
    (dfM <- AT['Model','df'])
    (dfE <- AT['Residuals','df'])
    (mse <- AT['Residuals','MS'])


    out <- list(cellMeans=outMeans, mse=mse, dfModel=dfM, dfError=dfE, var.equal=var.equal)

    class(out) <- 'cellMeans'

    return(out)

}#end cellMeans



#' @title Print \code{cellMeans} object
#'
#' @description Print \code{cellMeans} object.
#'
#' @param x A \code{cellMeans} object.
#' @param ... Currently ignored.
#'
#' @details see \code{\link{cellMeans}}
#'
#' @export
#'
#' @return Print for cellMeans object
#'
#'


print.cellMeans <- function(x, ...){

    print(x$cellMeans)
    
    cat('\n MSE =', x$mse, 'with df', x$dfModel, 'and', x$dfError, '\n')
    invisible(x)
} #end function


