###############################################################################
#### Calculation of EC50
###############################################################################
#' @title Calculates the EC50 value at steady-state
#' 
#' @description
#' Use the dose-respons curve to calculate the EC50 value 
#' concentration before the addition of oligonucleotide. 
#' 
#' @param KdOT The dissociation constant of the OT complex in nM
#' @param parms list of parameters
#' 
#' @return returns the EC50 value
#' @examples 
#' EC50(1,c(Et = 1,KdOT = 0.3,kOpT = 0.2,KdOTE = 70,
#' kOTpE = 5,vprod = 0.2,kdegrad = 0.04,alpha=0.1,kcleav = 8))
#' @export
EC50 <- function(KdOT,param=parms){ #
  param['KdOT'] <- KdOT
  Otseq <- 10^seq(-3,3,length=40)
  Trelseq <- Trel(Otseq,param=param)
  
  fit <- drm(Trelseq~Otseq,fct=LL.5())
  
  parms  <- coefficients(fit) # drm(Trelseq~Otseq,fct=LL.5())
#   out <- exp((parms[1]*log(parms[4])+
#                 log(exp(log((2*(-parms[3]+parms[2]))/(parms[2]-1))/parms[5])
#                     -1))/parms[1]) 
#   
  out <- summary(fit)[[3]]['e:(Intercept)',1]
  names(out) <- 'EC50'
  out
}

###############################################################################
#### Calculation of EC50 for a stochastic simulation of the ASO model
###############################################################################
#' @title Calculates the EC50 value for a stochastic simulation of the ASO model
#' 
#' @description
#' Use the stochastic dose-respons curve to calculate the EC50 value 
#' concentration before the addition of oligonucleotide. 
#' 
#' @param Trel Sequence of stochastic simulated total target levels as output from \code{Trelstoc}
#' @param Ot Sequence of total oligonucleotide concentrations used to calculate Trel
#' 
#' @return returns the EC50 value
#' @examples 
#' parms1 <- c(kOpT = 2E-5,kOTpE =50E-5 ,vprod = 150,  kdegrad = 0.04,  	  
#' kcleav = 2, kOT =0.06, kOTE=2, kC = 0.1)
#' #Initital state vector
#' x0 <- c(Tt=parms1["vprod"]/parms1["kdegrad"],
#'              OT=0,OTE=0,E=1e3,O=1e5,OCE=0,OC=0)
#' names(x0) <- c('Tt','OT','OTE','E','O','OCE','OC')
#' #Propensity vector
#' a <-  c("vprod","kOpT*O*Tt","kdegrad*Tt","kOT*OT","kOTE*OTE","kdegrad*OT",
#'          "kOTpE*OT*E","kdegrad*OTE","kcleav*OTE","kC*OC","kOTE*OCE" )
#' #State-change matrix
#' nu <- matrix(0,7,length(a))
#' dimnames(nu) <- list(names(x0),a)
#' nu['Tt',c('vprod','kOT*OT')] <- 1
#' nu['Tt',c('kOpT*O*Tt','kdegrad*Tt')] <- -1 
#' nu['OT',c('kOpT*O*Tt','kOTE*OTE')] <- 1
#' nu['OT',c('kOT*OT','kOTpE*OT*E','kdegrad*OT')] <- -1
#' nu['OTE',c('kOTpE*OT*E')] <- 1
#' nu['OTE',c('kOTE*OTE','kdegrad*OTE','kcleav*OTE')] <- -1
#' nu['E',c('kOTE*OTE','kdegrad*OTE','kOTE*OCE')] <- 1
#' nu['E',c('kOTpE*OT*E')] <- -1
#' nu['O',c('kOT*OT','kdegrad*OTE','kdegrad*OT','kC*OC')] <- 1
#' nu['O',c('kOpT*O*Tt')] <- -1
#' nu['OCE',c('kcleav*OTE')] <- 1
#' nu['OCE',c('kOTE*OCE')] <- -1
#' nu['OC',c('kOTE*OCE')] <- 1
#' nu['OC',c('kC*OC')] <- -1
#' 
#' Otseq <- 10^seq(2.5,6,by=0.2)
#' Trelseq <-sapply(Otseq,function(i) Trelstoc(i)$Tstat)})
#' EC50STOC <- EC50stoc(Trelseq, Otseq)
#' @export
EC50stoc <- function(Trel,Ot){ #
parms  <- coefficients(drm(Trel~Ot,fct=LL.5()))
  out <- exp((parms[1]*log(parms[4])+
                log(exp(log((2*(-parms[3]+parms[2]))/(parms[2]-1))/parms[5])-
                      1))/parms[1]) 
  names(out) <- 'EC50'
  out
}

###############################################################################
#### Calculation of EC50 for no coupling of off-rates
###############################################################################
#' @export
EC50NO <- function(KdOT,param=parmsNO){
  param['KdOT'] <- KdOT
  Otseq <- 10^seq(-3,3,length=50)
  Trelseq <- TrelNO(Otseq,param=param)
  parms  <- coefficients(drm(Trelseq~Otseq,fct=LL.5()))
  #b <- parms[1]; c <- parms[2]; d <- parms[3]; e <- parms[4]; f <- parms[5]
  out <- exp((parms[1]*log(parms[4])+log(exp(log((2*(-parms[3]+parms[2]))/(parms[2]-1))/parms[5])-1))/parms[1]) 
  names(out) <- 'EC50'
  out
}

