##anova table calculation function
#' @name anovaTable
#'
#' @title ANOVA Table
#'
#' @description Produces the overall ANOVA table where the model sum of squares are not partioned into their parts. 
#' 
#' @param object lm or aov model object
#' @param ... currently ignored
#' 
#'
#' @return Object of class anova and data.frame
#'
#' @export
#'
#' @import stats
#'
#' @examples
#' data(depression)
#' 
#' ## MLR model
#' modMLR <- lm(depress~trauma+control,data=depression)
#' anovaTable(modMLR)
#'
#' ## ANOVA model
#' depression$gender <- factor(depression$gender)
#' depression$history <- factor(depression$history)
#' modAOV  <- lm(depress~-1+gender+history+gender:history,data=depression)
#' anovaTable(modAOV)
#' 



anovaTable <- function(object,...){
    ##object <- modMLR
    ##object <- modAOV

    if(!any(class(object) %in% c('lm','aov'))){
        stop('The argument object must be an lm or aov class.')
    }#end if

    ## extract results
    yHat <- object$fitted.values
    eHat <- object$residuals
    w <- object$weights
    p  <- object$rank

    ## DF
    dfError <- object$df.residual
    dfModel <- p-attr(object$terms, "intercept")

    
    ## SS calc
    if (is.null(w)) {
        SSM <- if (attr(object$terms, "intercept")){ 
                   sum((yHat - mean(yHat))^2)
               } else { sum(yHat^2) }
        
        SSE <- sum(eHat^2)
        
    } else {
        
        SSM <- if (attr(object$terms, "intercept")) {
                   m <- sum(w * yHat/sum(w))
                   sum(w * (yHat - m)^2)
               } else { sum(w * f^2) }
        
        SSE <- sum(w * eHat^2)
        
    }#end else if

    f <- (SSM/dfModel)/(SSE/dfError)

    out <- data.frame(df=c(dfModel,dfError),SS=c(SSM,SSE))

    out$MS <- out$SS/out$df
    out$F <- c(f,NA)
    out$`Pr(>F)`  <- c(stats::pf(q=f,df1=dfModel,df2=dfError,lower.tail=FALSE),NA)
    row.names(out) <- c('Model','Residuals')
    class(out) <- c('anova','data.frame')

    return(out)
    
}#end anovaTable function





