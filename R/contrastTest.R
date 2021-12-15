
#' @name contrastTest
#'
#' @title Test Contrasts
#'
#' @description Contrast testing function. Designed to test contrasts of parameter estimates from a linear model. 
#' 
#' @param estVec numeric vector of parameter estimates for comparison
#' @param n numeric vector indicating the sample size for the parameter estimates, if a single value is given it is assumed to apply to all estiamtes
#' @param dfModel numeric value for the model degrees of freedom
#' @param dfError numeric value for the error or residual degrees of freedom
#' @param mse numeric value for the mean squared error from the model
#' @param C numeric matrix, each row is a contrast that should sum to zero, see details
#' @param test character, indicating which testing method should be used, see details
#' @param ... currently ignored
#'
#' @details The test argument can be one of the following: 'scheffe','bonferroni','hsd', or 'lsd'. 'hsd' is the Tukey HSD test. 'lsd' is th Fisher LSD test. The other two are the Scheffe test and Bonferroni adjustment.
#'
#' The matrix C is the contrast matrix. Each row is a separate contrast. The number of columns of C must be equal to the \code{length(estVec)}. Row names for C are retained in the output, but they are not required. 
#'
#' @return Object of class anova and data.frame
#'
#' @export
#'
#' @import stats
#' @import utils
#' 
#'
#' @examples
#' data(genericData)
#'
#' mod <- lm(Y~A,data=genericData)
#' dfModel <- anovaTable(mod)['Model','df']
#' dfError <- anovaTable(mod)['Residual','df']
#' mse <- anovaTable(mod)['Residual','MS']
#' meanVec <- aggregate(Y~A,FUN=mean,data=genericData)$Y
#' n <- aggregate(Y~A,FUN=length,data=genericData)$Y
#'
#' ## can add names for ease of interpretation with the output
#' names(meanVec) <- c('group 1','group 2','group 3')
#' contrastTest(estVec=meanVec,n=n,dfModel=dfModel,dfError=dfError,mse=mse,test='hsd')
#'
#' ## each group vs the mean of the other two
#' C <- rbind(c(1,-0.5,-0.5),c(-0.5,1,-0.5),c(-0.5,-0.5,1))
#' ## row names are not required but are helpful
#' row.names(C) <- c('1 vs 2+3','2 vs 1+3','3 vs 1+2')
#' contrastTest(estVec=meanVec,n=n,dfModel=dfModel,dfError=dfError,mse=mse,C=C,test='scheffe')
#' 


contrastTest <- function(estVec,n,dfModel,dfError,mse,C=NULL,test=c('scheffe','bonferroni','hsd','lsd'),...){


    ## ===================================================
    ## argument checking
    if(!is.numeric(estVec) || any(is.na(estVec)) || !is.vector(estVec) || length(estVec)<2){
        stop('The argument estVec must be a numeric vector with length >=2.')
    }#end if

    if(!is.numeric(n) || any(is.na(estVec)) || !is.vector(n) || any(n<=0)){
        stop('The argument n must be a postive numeric vector.')
    }#end if

    if(length(n)==1){
        nVec <- rep(n,length(estVec))
    }else{
        nVec <- n
    }#end else if

    if(length(nVec) != length(estVec)){
        stop('The argument n must be length 1 or the same length as estVec.')
    }#end if


    modList <- list(dfModel=dfModel,dfError=dfError,mse=mse)

    for(i in seq_along(modList)){
        thisVal <- modList[[i]]
        thisName <- names(modList)[i]

        if(!is.numeric(thisVal) || !is.vector(thisVal) || length(thisVal)!=1 || thisVal<=0){
            stop('The argument ',thisName,' must be a single positive numeric value.')
        }#end if
    }#end for i



    testMethod <- match.arg(test)
    ##testMethod <- 'hsd'


    if(is.null(names(estVec))){
        names(estVec) <- paste0('est',1:length(estVec))
    }#end if


    ## ===================================================
    ## create constract matrix if none provide
    ## does all pairwise comparisons 
    if(is.null(C)){
        
        combResult <- utils::combn(length(estVec) ,2)

        contrastNames <- NULL
        for(i in 1:ncol(combResult)){
            tmpContrast <- rep(0,length(estVec))
            tmpContrast[combResult[1,i]] <- 1
            tmpContrast[combResult[2,i]] <- -1
            C <- rbind(C,tmpContrast)

            contrastNames <- c(contrastNames,paste0(names(estVec)[combResult[1,i]],' vs ',names(estVec)[combResult[2,i]]))
         
        }#end for i

        row.names(C) <- contrastNames
    }#end if
    

    if(!is.matrix(C) || !is.numeric(C) || !all(rowSums(C)==0) || ncol(C) != length(estVec)){
        stop('The argument C must either be NULL or a numeric matrix. Each row corresponds to a contrast to be tested. The number of columns must be equal to the number of total number of parameters considered.')
    }#end if


    if(is.null(row.names(C))){
        row.names(C) <- paste0('Contrast ',1:nrow(C))
    }#end if

    ## contrast estimates
    est <- C%*%estVec

    
    ## generic constrast variance
    estVar <- (C^2)%*%(1/nVec) *mse 


    if(testMethod=='scheffe'){
        
        out <- data.frame(estimate=est,df1=dfModel,
                          df2=dfError,F=est^2/(estVar*dfModel))
        out[,'Pr(>F)'] <- stats::pf(q=out$F,df1=dfModel,df2=dfError,lower.tail=FALSE)

    }else if(testMethod=='hsd'){

        out <- data.frame(estimate=est,df1=length(estVec),
                          df2=dfError,q=abs(est)/sqrt(estVar/2))
        out[,'Pr(>q)'] <- stats::ptukey(q=out$q,nmeans=length(estVec),df=dfError,
                                 lower.tail=FALSE)
    
    }else{

        out <- data.frame(estimate=est,df=dfError,t=est/sqrt(estVar))
        out[,'Pr(>|t|)'] <- stats::pt(q=abs(out$t),df=dfError,lower.tail=FALSE)*2

        if(testMethod=='bonferroni'){
            out[,'Pr(>|t|)'] <- out[,'Pr(>|t|)']*nrow(out)
            if(any(out[,'Pr(>|t|)']>1)){
                message('Bonferroni adjusted p-value that are greater than 1 are set to 1')
                out[out[,'Pr(>|t|)']>1,'Pr(>|t|)'] <- 1
                }#end if
        }#end if

    }#end else if


    class(out) <- c('anova','data.frame')
    return(out)

}# contrastTest function
