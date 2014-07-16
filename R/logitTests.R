#' Association testing under an additive model
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @return The p-value of a logistic regression of case/control status onto genotype.
logit.additive=function(x,status)
{
    if(length(x) != length(status))
        {
            stop("logit.additive: length(x) != length(status)")
            return;
        }
    GLM=glm( status ~ x, family=binomial("logit") )
    return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
}

#' Association testing under a recessive model
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @note The test works by recoding 1 genotypes as 0 (e.g., one copy of minor allele is the same as none)
#' @return The p-value of a logistic regression of case/control status onto genotype.
logit.recessive=function(x,status)
    {
        if(length(x) != length(status))
            {
                stop("logit.additive: length(x) != length(status)")
                return;
            }
        x[which(x==1)]=0
        GLM=glm( status ~ x, family=binomial("logit") )
        return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
    }

#' Association testing under a dominant model
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @note The test works by recoding 2 genotypes as 1 (e.g., one copy of minor allele is the same as two)
#' @return The p-value of a logistic regression of case/control status onto genotype.
logit.dominant=function(x,status)
    {
        if(length(x) != length(status))
            {
                stop("logit.dominant: length(x) != length(status)")
                return;
            }
        x[which(x==2)]=1
        GLM=glm( status ~ x, family=binomial("logit") )
        return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
    }

#' Obtain single-marker p-values for case/control data
#' @param genos A matrix of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal nrow(x)
#' @param model One of "additive","recessive", or "dominant"
#' @return An array of p-values of logistic regressions of case/control status onto genotype.  The order of the p-values corresponds to the column order in x.
#' @examples
#' data(rec.ccdata)
#' #The function works on the genotype matrix
#' #We assign discrete phenotypes based on number of controls & cases in ccblock.dsdata$genos:
#' pvals.additive = ccpvals(rec.ccdata$genos,c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases)))
#' pvals.recessive = ccpvals(rec.ccdata$genos,c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases)),"recessive")
#' pvals.dominant = ccpvals(rec.ccdata$genos,c(rep(0,rec.ccdata$ncontrols),rep(1,rec.ccdata$ncases)),"dominant")
ccpvals = function(genos,status,model="additive")
    {
        if( length(status) != nrow(genos) )
            {
                stop("ccpvals.additive: nrow(x) != length(status)")
                return;
            }
        if( model == "additive" )
            {
                return(apply(genos,2,logit.additive,status))
            }
        else if ( model == "recessive" )
            {
                return(apply(genos,2,logit.recessive,status))
            }
        else if (model == "dominant")
            {
                return(apply(genos,2,logit.dominant,status))
            }
        else
            {
                stop(paste("ccpvals: model",model,"is not valid. Model must be one of additive, dominant, or recessive"))
            }
        return();
    }
