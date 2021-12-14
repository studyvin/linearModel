
#' @name contrastTest
#'
#' @title Test Contrasts
#'
#' @description 
#' 
#' @param object lm or aov model object
#' @param trt character vector indicating the categorical predictor variables
#' @param C numeric matrix, each row is a contrast that should sum to zero
#' @param test character, indicating which testing method should be used
#' @param ... currently ignored
#' 
#'
#' @return 
#'
#' @export
#'
#'
#'
#' @examples
#' data(depression)

load('../data/depression.RData')

head(depression)

depression$gender <- ifelse(depression$gender,'F','M')
depression$history <- ifelse(depression$history,'Y','N')
f <- 'depress~history+gender+trauma+gender*trauma + welbeing + welbeing:gender'
object <- lm(formula(f),data=depression)
summary(object)
anova(object)

coef(object)[c(3)] ## female trauma
sum(coef(object)[c(3,5)]) ## male trauma
coef(object)[c(4)] ## female wellbeing
sum(coef(object)[c(4,6)]) ## male wellbeing


agricolae::HSD.test(object,trt=trt,console=TRUE)

trt <- c('gender','history')
covar <- c('trauma','welbeing')
test='hsd'

TukeyHSD(object,ordered=FALSE)
contrastTest(object,trt=trt,test='hsd')

model.tables(object)


contrastTest <- function(object,trt,C=NULL,test=c('scheffe','bonferroni','hsd','lsd'),covar=NULL,...){


    if(is.null(trt) & is.null(C)){
        stop('Both trt and C cannot be NULL, one must be specified.')
    }#end if

    if(!any(class(object) %in% c('lm','aov'))){
        stop('The argument object must be an lm or aov class.')
    }#end if

    if(!is.character(trt) | !is.vector(trt)){
        stop('The argument trt must be a character vector indicating the categorical predictor variables in the model object.')
    }#end if


    testMethod <- match.arg(test)
    ##testMethod <- 'hsd'

    
    (allCovars <- unique(unlist(strsplit(attr(object$terms,'term.labels'),'\\:'))))
    
    if(any(grepl('I\\(',allCovars))){
        stop('The argument object cannot have any transformations within the formula statement. \n For example: Y ~ X + I(log(A)) is not allowed.')
    }#end if



    if(!is.null(trt)){
        badVar <- trt[!trt%in%allCovars]
        if(length(badVar)>0){
            stop('The following are not variables in the model object: ',paste0(badVar,collapse=', '),'\n Please change the trt argument.')
        }#end if
    }#end if

    if(!is.null(covar)){
        if(is.null(trt) || !is.character(covar) || !is.vector(covar)){
            stop('The argument covar can only be specified when trt is also specified and covar must be a character vector.')
        }#end if

        badVar <- covar[!covar%in%allCovars]
        if(length(badVar)>0){
            stop('The following are not variables in the model object: ',paste0(badVar,collapse=', '),'\n Please change the covar argument.')
        }#end if

    }#end if

    ## ===================================================
    ## internal functions
    makeVec <- function(d,r,t,FUN){
        ##x <- meanDF
        ##r <- responseVar
        x <- aggregate(formula(paste0(r,'~',paste0(t,collapse='+'))),FUN=FUN,data=d)
        out <- x[,r]
        names(out) <- apply(x[,t,drop=FALSE],1,paste0,collapse=':')
        return(out)
    }#end makeVec



    ## ===================================================
    ## extract info from model object
    df <- object$model
    dfError <- object$df.residual
    mse <- sum(object$residuals^2)/dfError
    dfModel <- object$rank-attr(object$terms, "intercept")


    responseVar <- names(df)[!names(df)%in%allCovars]

    if(length(responseVar) !=1){
        stop('Problem with object argument. Only one response variable is allowed.')
    }#end if


    (nVec <- makeVec(df,r=responseVar,trt,FUN=length))


    if(!is.null(covar)){

        newData <- NULL
        for(i in 1:length(covar)){

            newData <- cbind(newData,makeVec(df,covar[i],trt,FUN=mean))
        }#end for i
        newData <- data.frame(newData)
        names(newData) <- covar

        rowSplit <- strsplit(row.names(newData),'\\:')
        for(i in 1:length(trt)){
            newData[,trt[i]] <- unlist(lapply(rowSplit,'[[',i))
        }#end for i

        newData

    }else{

        (betaVec <- makeVec(df,responseVar,trt,FUN=mean))
    
    }#end else if

    
    ## create constract matrix if none provide
    ## does all pairwise comparisons 
    if(is.null(C)){

        ##J$ this does not make sense if testing the beta coefficients
        combResult <- combn(length(betaVec) ,2)

        contrastNames <- NULL
        for(i in 1:ncol(combResult)){
            tmpContrast <- rep(0,length(meanVec))
            tmpContrast[combResult[1,i]] <- 1
            tmpContrast[combResult[2,i]] <- -1
            C <- rbind(C,tmpContrast)

            contrastNames <- c(contrastNames,paste0(names(meanVec)[combResult[1,i]],' vs ',names(meanVec)[combResult[2,i]]))
         
        }#end for i

        row.names(C) <- contrastNames
    }#end if
    

    if(!is.matrix(C) || !is.numeric(C) || !all(rowSums(C)==0) || ncol(C) != length(meanVec)){
        stop('The argument C must either be NULL or a numeric matrix. Each row corresponds to a contrast to be tested. The number of columns must be equal to the number of total number of parameters considered.')
    }#end if


    if(is.null(row.names(C))){
        row.names(C) <- paste0('Contrast ',1:nrow(C))
    }#end if

    ## contrast estimates
    est <- C%*%betaVec

    ## generic constrast variance
    estVar <- (C^2)%*%(1/nVec) *mse 


    if(testMethod=='scheffe'){
        
        out <- data.frame(estimate=est,df1=length(betaVec)-1,
                          df2=dfError,F=est^2/(estVar*(length(betaVec)-1)))
        out[,'Pr(>F)'] <- pf(q=out$F,df1=out$df1,df2=out$df2,lower.tail=TRUE)

    }else if(testMethod=='hsd'){

        out <- data.frame(estimate=est,df1=length(betaVec),
                          df2=dfError,q=abs(est)/sqrt(estVar/2))
        out[,'Pr(>q)'] <- ptukey(q=out$q,nmeans=length(betaVec),df=dfError,
                                 lower.tail=FALSE)
    
    }else{

        out <- data.frame(estimate=est,df=dfError,t=est/sqrt(estVar))
        out[,'Pr(>|t|)'] <- pt(q=abs(out$t),df=dfError,lower.tail=FALSE)

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
