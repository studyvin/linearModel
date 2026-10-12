

#' @name contrast
#'
#' @title Contrasts
#'
#' @description Contrast testing function. Designed to test contrasts of parameter estimates from a linear model. 
#' 
#' @param cm \code{\link[linearModel]{cellMeans}} object
#' @param C numeric matrix, each row is a contrast that should sum to zero, see details
#' @param test character, indicating which testing method should be used, see details
#' @param ... currently ignored
#'
#' @details The test argument can be one of the following: 'scheffe','bonferroni','hsd', or 'lsd'. 'hsd' is the Tukey HSD test. 'lsd' is th Fisher LSD test. The other two are the Scheffe test and Bonferroni adjustment.
#'
#' The matrix C is the contrast matrix. Each row is a separate contrast. The number of columns of C must be equal to the \code{length(estVec)}. Row names for C are retained in the output, but they are not required.
#'
#' This is a wrapper function and the workhorse is the \code{\link[linearModel]{contrastTest}} function. 
#'
#' @return Object of class anova and data.frame
#'
#' @seealso \code{\link[linearModel]{cellMeans}}
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





contrast <- function(cm, C=NULL,test = c("hsd", "scheffe", 
    "bonferroni", "lsd"),...){


    if(!methods::is(cm, 'cellMeans')){
        stop('Argument cm must be the result from a call to linearModel::cellMeans')
    }#end



    df <- cm$cellMeans

    estVec <- df$mean

    names(estVec) <- apply(df[,!names(df)%in%c('mean','var','r'), drop=FALSE], 1, paste0, collapse='_')
    estVec
    
 
    out <- linearModel::contrastTest(estVec=estVec,
                                     n=cm$cellMeans$r,
                                     dfModel=cm$dfModel,
                                     dfError=cm$dfError,
                                     mse=cm$mse,                                     
                                    C=C, test=test)

    return(out)

}#end function
