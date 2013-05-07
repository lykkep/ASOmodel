library(deSolve)
library(rootSolve)
library(nleqslv)
library(TeachingDemos)
library(sfsmisc)
library(drc)

last <- function(x) { tail(x, n = 1) }

#' Calculates the relative Target total target concentration at steady state
#' 
#' Use the algebraic solution to the system to the total target concentration at steady state relative to the target 
#' concentration before the addition of oligonucleotide. 
#' @param otot Total concentration of oligonucleotide added to the system in nM
#' @param parms list of parameters
#' 
#' @return returns the relative total target concentration
#' @example #TODO Please give working example code that does something useful
Trel <- function(otot,param=parms){  #
  
  #### Parameters
  k1 =param['kOpT']; D1 = param['KdOT']; k2 = param['kOTpE']; Et <- param['Et']
  D2 = param['KdOTE']; vt = param['vprod']; k4 = param['vdegrad']
  alpha=param['alpha']; kE=param['kcleav'];
  
  tmp <- function(Ot){

  #### Polynomial oefficients
  n4 <- kE^2*alpha*D2^3*k2^3+2*kE^3*alpha*D2^2*k2^2+kE^4*alpha*D2*k2+k4*k1*D1*kE^3+k4*kE*alpha*D2^3*k2^3+2*k4*kE^2*alpha*D2^2*k2^2+k4*kE^3*alpha*D2*k2+k4*k1*D1*D2^3*k2^3+k1*D1*D2^3*k2^3*kE+3*k1*D1*D2^2*k2^2*kE^2+3*k1*D1*D2*k2*kE^3+k1*D1*kE^4+3*D2^2*k2^2*k1*D1*k4*kE+3*k4*k1*D1*D2*k2*kE^2

  n0 <- vt*k1*D1*Ot*D2^3*k2^3*Et^2
  
  n1 <- -2*vt*k1*D1*Ot*D2^3*k2^3*Et-D2^4*k4^2*k2^3*D1*Et-D2^3*k4^3*k2^2*D1*Et-vt*k1*D1*D2^3*k2^3*Et^2-D2^3*D1^2*Et*k1*k2^2*k4*kE-D2^3*k1*D1*k2^2*vt*Et*k4-D2^3*D1^2*Et*k1*k2^2*k4^2-vt*kE*alpha*D2^3*k2^3*Et^2-D2^3*k4^2*k2^2*k1*D1*Ot*Et-D2^4*k4*k2^3*k1*D1*Ot*Et-D2^3*k1*D1*k2^2*vt*Et*kE-D2^4*k1*D1*k2^3*vt*Et-k4*D1*kE*D2^3*k2^3*Et^2-vt*k1*D1*kE*k2^2*Et^2*D2^2-k4^2*k2^3*D1*D2^3*Et^2-k4*k1*D1*Ot*D2^3*k2^3*Et^2-k1*D1*Ot*D2^3*k2^3*kE*Et^2-2*vt*k1*D1*Ot*D2^2*k2^2*Et*kE-D2^4*D1^2*Et*k1*k2^3*k4-D2^3*k4*k2^2*k1*D1*Ot*Et*kE-D2^3*k4^2*k2^2*D1*Et*kE
  
  n2 <- D1^2*D2^2*k1*k2*k4*kE^2+2*D1^2*D2^3*k1*k2^2*k4*kE+kE^2*alpha*D2^3*k2^3*Et^2+D2^4*k1*D1*k2^3*vt+2*D2^3*k4^2*k2^2*D1*kE+D2^2*k4^3*k2*D1*kE+D2^2*k4^2*k2*D1*kE^2+D2^5*k4*k1*D1*k2^3+2*D2^4*k4^2*k1*D1*k2^2+D2^3*k4^3*k1*D1*k2+D2^4*k4*k2^3*D1^2*k1+2*k1*D1*D2^3*k2^3*vt*Et+2*k4^2*D2^2*k2^2*D1*Et*kE+4*k1*D1*D2^2*k2^2*vt*Et*kE+2*k4*D2^3*k2^3*k1*D1*Ot*Et+D1^2*D2^3*k1*k2^2*k4^2+2*k4*D2^2*k2^2*k1*D1*Ot*Et*kE+2*k4^2*D2^3*k2^3*D1*Et+D2^4*k4^2*k2^3*D1+D2^3*k4^3*k2^2*D1+2*vt*k1*D1*kE^2*k2*Et*D2+2*vt*k1*D1*Ot*D2^2*k2^2*kE+2*k1*D1*Ot*D2^3*k2^3*kE*Et+2*k1*D1*Ot*D2^2*k2^2*kE^2*Et+vt*k1*D1*Ot*D2*k2*kE^2+2*D2^3*k4*k2^2*k1*D1*Ot*kE+k4*k1*D1*kE*k2^2*Et^2*D2^2+4*D2^3*k4*k1*D1*k2^2*Et*kE+D2^2*k1*D1*k2*vt*kE*k4+D2^2*k4^2*k1*D1*kE*k2*Et+D2^2*k4*k1*D1*kE^2*k2*Et+D2^2*k4^2*k2*k1*D1*Ot*kE+D2^2*k4*k2*k1*D1*Ot*kE^2+k1*D1*kE^2*k2^2*Et^2*D2^2+2*k4*D1*kE*D2^3*k2^3*Et+2*k4*D1*kE^2*D2^2*k2^2*Et+D2^3*k4*k1*D1*k2*kE^2+vt*k1*D1*Ot*D2^3*k2^3+k4*kE*alpha*D2^3*k2^3*Et^2+D2^3*k1*D1*k2^2*vt*k4+D2^4*k4*k2^3*k1*D1*Ot+D2^3*k4^2*k2^2*k1*D1*Ot+2*D2^4*k4*k1*D1*k2^3*Et+2*D2^3*k4^2*k1*D1*k2^2*Et+D2^4*k1*D1*k2^3*kE*Et+D2^3*k1*D1*k2^2*kE^2*Et+2*D2^4*k4*k1*D1*k2^2*kE+2*D2^3*k4^2*k1*D1*k2*kE+2*vt*kE*alpha*D2^3*k2^3*Et+2*vt*kE^2*alpha*D2^2*k2^2*Et+k4*k1*D1*D2^3*k2^3*Et^2+k1*D1*D2^3*k2^3*kE*Et^2+D2^2*k4^2*k2*D1^2*k1*kE+2*D2^3*k1*D1*k2^2*vt*kE+D2^2*k1*D1*k2*vt*kE^2+D2^4*k4*kE*alpha*k2^3*Et+D2^3*k4^2*kE*alpha*k2^2*Et+D2^3*k4*kE^2*alpha*k2^2*Et
    
  n3 <- -D2^4*k4*kE*alpha*k2^3-D2^3*k4^2*kE*alpha*k2^2-2*D2^3*k4*kE^2*alpha*k2^2-D2^2*k4^2*kE^2*alpha*k2-D2^2*k4*kE^3*alpha*k2-2*D2^4*k4*k1*D1*k2^3-2*D2^3*k4^2*k1*D1*k2^2-D2^4*k1*D1*k2^3*kE-k4^2*k2*D1*D2*kE^2-2*D2^3*k1*D1*k2^2*kE^2-D2^2*k1*D1*k2*kE^3-D2*k4*k1*D1*kE^3-D2*k4^2*k1*D1*kE^2-vt*kE*alpha*D2^3*k2^3-2*vt*kE^2*alpha*D2^2*k2^2-vt*kE^3*alpha*D2*k2-k4*D1*kE*D2^3*k2^3-2*k4*D1*kE^2*D2^2*k2^2-k4*D1*kE^3*D2*k2-vt*k1*D1*D2^3*k2^3-2*kE^2*alpha*D2^3*k2^3*Et-2*kE^3*alpha*D2^2*k2^2*Et-2*D2^2*k2^2*D1*k4^2*kE-k4^2*k2^3*D1*D2^3-vt*k1*D1*kE^3-4*k4*k1*D1*D2^2*k2^2*Et*kE-2*k4*k1*D1*Ot*D2^2*k2^2*kE-k4*k1*D1*Ot*D2*k2*kE^2-2*k4*k1*D1*kE^2*k2*Et*D2-4*k1*D1*D2^2*k2^2*kE^2*Et-2*k1*D1*kE^3*k2*Et*D2-k4*k1*D1*Ot*D2^3*k2^3-k1*D1*Ot*D2^3*k2^3*kE-2*k1*D1*Ot*D2^2*k2^2*kE^2-k1*D1*Ot*D2*k2*kE^3-2*k4*kE*alpha*D2^3*k2^3*Et-3*D2^2*k4^2*k1*D1*kE*k2-5*D2^2*k4*k1*D1*kE^2*k2-2*k4*kE^2*alpha*D2^2*k2^2*Et-3*vt*k1*D1*kE*D2^2*k2^2-3*vt*k1*D1*kE^2*D2*k2-2*k4*k1*D1*D2^3*k2^3*Et-2*k1*D1*D2^3*k2^3*kE*Et-6*D2^3*k4*k1*D1*kE*k2^2
  
  ROOT <- polyroot(c(n0,n1,n2,n3,n4))
  #### [OTE] at steady-state
  TrelSS <- -k4*(-2*ROOT^2*k4*D1*kE^2*D2^2*k2^2*Et-2*ROOT^2*D2^3*k4*k1*D1*k2*kE^2+ROOT^2*vt*k1*D1*Ot*D2^3*k2^3+ROOT^2*D2^3*k1*D1*k2^2*vt*k4-ROOT^2*D2^4*k1*D1*k2^3*kE*Et-2*ROOT^2*D2^3*k1*D1*k2^2*kE^2*Et-2*ROOT^2*D2^4*k4*k1*D1*k2^2*kE-ROOT^2*D2^3*k4^2*k1*D1*k2*kE+2*ROOT^2*vt*kE*alpha*D2^3*k2^3*Et+2*ROOT^2*vt*kE^2*alpha*D2^2*k2^2*Et+2*ROOT^2*D2^3*k1*D1*k2^2*vt*kE+ROOT^2*D2^2*k1*D1*k2*vt*kE^2-ROOT^2*D2^3*k4*kE^2*alpha*k2^2*Et-ROOT*vt*k1*D1*D2^3*k2^3*Et^2-ROOT*vt*kE*alpha*D2^3*k2^3*Et^2-ROOT*D2^4*k1*D1*k2^3*vt*Et+ROOT*k4*D1*kE*D2^3*k2^3*Et^2+ROOT*D2^3*k4^2*k2^2*D1*Et*kE+2*ROOT^3*D2^2*k4*k1*D1*kE^2*k2-3*ROOT^3*vt*k1*D1*kE*D2^2*k2^2-3*ROOT^3*vt*k1*D1*kE^2*D2*k2+ROOT^3*D2^3*k4*k1*D1*kE*k2^2-ROOT^2*D1^2*D2^2*k1*k2*k4*kE^2-ROOT^2*D1^2*D2^3*k1*k2^2*k4*kE-D2^2*k2*ROOT^2*D1*Et*k1*kE^3-D2^2*k2*ROOT^2*D1*Ot*k1*kE^3-D2^4*k2^3*ROOT^2*D1*Ot*k1*kE-2*D2^3*k2^2*ROOT^2*D1*Ot*k1*kE^2+D2^3*k2^2*Et*ROOT*D1*k4*kE^2+D2^4*k2^3*Et*ROOT*D1*k4*kE+D2^3*k2^2*Et*ROOT*D1^2*k1*kE^2+D2^4*k2^3*Et*ROOT*D1^2*k1*kE+2*ROOT^2*k1*D1*D2^3*k2^3*vt*Et-2*ROOT^2*k4*D1*kE*D2^3*k2^3*Et+ROOT^3*D2^4*alpha*k2^3*kE^2+ROOT^3*D1*D2^3*k2^3*kE^2+2*ROOT^3*D2^3*alpha*k2^2*kE^3+2*ROOT^3*D1*D2^2*k2^2*kE^3+ROOT^3*D2^2*alpha*k2*kE^4+ROOT^3*D1*D2*k1*kE^4+ROOT^3*D1*D2*k2*kE^4-ROOT^3*vt*k1*D1*kE^3+vt*k1*D1*Ot*D2^3*k2^3*Et^2+ROOT^2*D2^2*k1*D1*k2*vt*kE*k4-ROOT*D2^3*k1*D1*k2^2*vt*Et*kE+ROOT^2*vt*k1*D1*Ot*D2*k2*kE^2+4*ROOT^2*k1*D1*D2^2*k2^2*vt*Et*kE+2*ROOT^2*vt*k1*D1*kE^2*k2*Et*D2-ROOT^2*D2^2*k4*k2*k1*D1*Ot*kE^2+2*ROOT^2*vt*k1*D1*Ot*D2^2*k2^2*kE-ROOT^2*D2^3*k4*k1*D1*k2^2*Et*kE-ROOT*D2^3*k1*D1*k2^2*vt*Et*k4-ROOT^2*D2^3*k4*k2^2*k1*D1*Ot*kE+D2^4*k2^3*Et*ROOT*D1*Ot*k1*kE+ROOT*D2^3*D1^2*Et*k1*k2^2*k4*kE-ROOT*vt*k1*D1*kE*k2^2*Et^2*D2^2+D2^3*k2^2*Et*ROOT*D1*Ot*k1*kE^2-ROOT^2*D2^2*k4*k1*D1*kE^2*k2*Et-2*ROOT*vt*k1*D1*Ot*D2^3*k2^3*Et-2*ROOT*vt*k1*D1*Ot*D2^2*k2^2*Et*kE+ROOT*D2^3*k4*k2^2*k1*D1*Ot*Et*kE-D2^3*k2*ROOT^2*D1*k1*kE^3-D2^5*k2^3*ROOT^2*k1*D1*kE-2*D2^4*k2^2*ROOT^2*k1*D1*kE^2-D2^4*k2^3*ROOT^2*D1^2*k1*kE-2*D2^3*k2^2*ROOT^2*D1^2*k1*kE^2+D2^3*k2^3*Et^2*ROOT*D1*kE^2+ROOT^2*D2^4*k1*D1*k2^3*vt-ROOT^2*D2^3*k4^2*k2^2*D1*kE-ROOT^2*D2^2*k4^2*k2*D1*kE^2+ROOT^3*D2^3*k4*kE^2*alpha*k2^2+ROOT^3*D2^2*k4*kE^3*alpha*k2+ROOT^3*D2^4*k1*D1*k2^3*kE+3*ROOT^3*D2^3*k1*D1*k2^2*kE^2+3*ROOT^3*D2^2*k1*D1*k2*kE^3+ROOT^3*D2*k4*k1*D1*kE^3-ROOT^3*vt*kE*alpha*D2^3*k2^3-2*ROOT^3*vt*kE^2*alpha*D2^2*k2^2-ROOT^3*vt*kE^3*alpha*D2*k2+ROOT^3*k4*D1*kE*D2^3*k2^3+2*ROOT^3*k4*D1*kE^2*D2^2*k2^2+ROOT^3*k4*D1*kE^3*D2*k2-ROOT^3*vt*k1*D1*D2^3*k2^3-D2^4*k2^3*ROOT^2*Et*alpha*kE^2-D2^4*k2^3*ROOT^2*D1*k4*kE-2*D2^3*k2^3*ROOT^2*D1*Et*kE^2-D2^3*k2^2*ROOT^2*Et*alpha*kE^3-2*D2^3*k2^2*ROOT^2*D1*k4*kE^2-2*D2^2*k2^2*ROOT^2*D1*Et*kE^3-D2^2*k2*ROOT^2*k4*D1*kE^3-D2^2*k2*ROOT^2*D1^2*k1*kE^3)/((kE+k4)*(kE*ROOT*alpha*D2^2*k2^2*Et-kE*ROOT^2*alpha*D2^2*k2^2-kE^2*ROOT^2*alpha*D2*k2-k1*D1*Ot*D2^2*k2^2*Et+k1*D1*Ot*D2^2*k2^2*ROOT+k1*D1*Ot*D2*k2*kE*ROOT+k1*D1*ROOT*D2^3*k2^2+k4*k1*D1*D2^2*k2*ROOT+k1*D1*ROOT*D2^2*k2*kE+k1*D1*ROOT*D2^2*k2^2*Et-k1*D1*ROOT^2*D2^2*k2^2-2*k1*D1*ROOT^2*D2*k2*kE+k1*D1*kE*ROOT*k2*Et*D2-k1*D1*kE^2*ROOT^2)*vt*(k2*Et*D2-D2*k2*ROOT-kE*ROOT))
  OTSS <- D2*ROOT*(D2*k2+k4+kE)/(k2*Et*D2-D2*k2*ROOT-kE*ROOT)
  ESS <- (k2*Et*D2-D2*k2*ROOT-kE*ROOT)/(D2*k2)
  OTESS <- ROOT
  OCESS <- kE*ROOT/(D2*k2)
  OSS <- (-k1*D1*ROOT*D2^3*k2^2+(((-k1*D1-kE*alpha)*ROOT+k1*D1*Ot)*(Et-ROOT)*k2-k1*D1*ROOT*(kE+k4))*k2*D2^2-kE*((-kE*alpha-2*k1*D1)*ROOT+D1*k1*(Ot+Et))*ROOT*k2*D2+k1*D1*kE^2*ROOT^2)/(D1*(k2*(Et-ROOT)*D2-kE*ROOT)*D2*k2*k1)
  OCSS <- kE*ROOT*alpha/(k1*D1)
  
  Re(TrelSS[ Re(TrelSS) <= 1 & Re(TrelSS) >= 0 & Im(TrelSS) < 1e-5 & 
               Re(OSS) >=0 & Re(OTSS) >= 0 & 
               Re(OTESS) >=0 & Re(OCESS) >=0 & Re(OCSS) >= 0 & Re(ESS) >=0] )
  }

  sapply(otot,tmp)
}


###############################################################################
#### IC50 calculation
###############################################################################
#TODO: document
IC50 <- function(KdOT,param=parms){
  param['KdOT'] <- KdOT
  Otseq <- 10^seq(-3,3,length=50)
  Trelseq <- Trel(Otseq,param=param)
  parms  <- coefficients(drm(Trelseq~Otseq,fct=LL.5()))
  #b <- parms[1]; c <- parms[2]; d <- parms[3]; e <- parms[4]; f <- parms[5]
  out <- exp((parms[1]*log(parms[4])+log(exp(log((2*(-parms[3]+parms[2]))/(parms[2]-1))/parms[5])-1))/parms[1]) 
  names(out) <- 'IC50'
  out
}

#TODO: document
IC50stoc <- function(Trel,Ot){ # TODO: more descriptive function name?
  parms  <- coefficients(drm(Trel~Ot,fct=LL.5()))
  #b <- parms[1]; c <- parms[2]; d <- parms[3]; e <- parms[4]; f <- parms[5]
  out <- exp((parms[1]*log(parms[4])+log(exp(log((2*(-parms[3]+parms[2]))/(parms[2]-1))/parms[5])-1))/parms[1]) 
  names(out) <- 'IC50'
  out
}

IC50NO <- function(KdOT,param=parmsNO){
  param['KdOT'] <- KdOT
  Otseq <- 10^seq(-3,3,length=50)
  Trelseq <- TrelNO(Otseq,param=param)
  parms  <- coefficients(drm(Trelseq~Otseq,fct=LL.5()))
  #b <- parms[1]; c <- parms[2]; d <- parms[3]; e <- parms[4]; f <- parms[5]
  out <- exp((parms[1]*log(parms[4])+log(exp(log((2*(-parms[3]+parms[2]))/(parms[2]-1))/parms[5])-1))/parms[1]) 
  names(out) <- 'IC50'
  out
}

###############################################################################
#### The ODEs
###############################################################################
#TODO: document this function
diffASO <- function(t,y, param){
  k1 =param['kOpT']; D1 = param['KdOT']; k2 = param['kOTpE']; Et <- param['Et']
  D2 = param['KdOTE']; vt = param['vprod']; k4 = param['vdegrad']
  alpha=param['alpha']; kE=param['kcleav'];
  km1 = k1*D1; km2=k2*D2; k3 = km1/alpha
  
  Tt = y[1]; OT = y[2]; OTE= y[3]; E=y[4]; O= y[5]
  OCE= y[6]; OC = y[7]
  
  F <- {}   
  #1= T
  F[1]= vt-k1*O*Tt-k4*Tt+km1*OT
  #2= OT
  F[2]= k1*Tt*O-km1*OT-k2*OT*E+km2*OTE-k4*OT
  #3= OTE
  F[3]= k2*OT*E-km2*OTE-k4*OTE-kE*OTE
  #4= E
  F[4]= -k2*(OT)*E+km2*(OTE+OCE)+k4*(OTE)
  #5= O
  F[5]=km1*OT-k1*Tt*O+k4*(OT+OTE)+k3*OC
  #OCE
  F[6] = kE*OTE-km2*OCE-k3*OCE
  #OC
  F[7]= +km2*OCE-k3*OC
  
  return(list(F))
}

###############################################################################
#### Plot dose-respons for various parameters
###############################################################################
#TODO: document
plot.doseresponse <- function(par,list.par,title,unit,plot=T){
  if(plot){
  curve(Trel,1E-2,1E4,log='x', 
        xlab=expression(Total~oligo~concentration~'(nM)'),
        ylab=expression(T[rel]), 
        las=1,ylim=c(0,1),xaxt='n',mgp=c(2.2,0.7,0))
  axis(1,at=10^c(-2,0,2,4),labels=pretty10expLP(10^c(-2,0,2,4),drop.1=T))
  for(i in 1:length(list.par)){
    parmsM <- parms ; parmsM[par] <- list.par[i]
    runM <- function(x)Trel(x,param=parmsM)
    curve(runM,1E-2,1E4,col=rainbow(length(list.par))[i],add=TRUE)
  }
  Legend=as.expression(sapply(c(parms[par],list.par),
                              function(x) substitute(x*y,list(x=x,y=unit))))
  legend("bottomleft", lwd=2, col=c('black',rainbow(length(list.par))),
         legend=Legend, bty="n",cex=0.6)
  legend(2E-3,0.6,title,bty='n',yjust=0.55,cex=0.7)
  }
}

###############################################################################
#### Pretty logarithmic axis-labels
###############################################################################
#TODO: document
pretty10expLP <-  function (x, drop.1 = FALSE, digits.fuzz = 7) 
{
  eT <- floor(log10(abs(x)) + 10^-digits.fuzz)
  mT <- signif(x/10^eT, digits.fuzz)
  ss <- vector("list", length(x))
  for (i in seq(along = x)) ss[[i]] <- if (x[i] == 0) 
    quote(0)
  else if (drop.1 && mT[i] == 1) 
    substitute(10^E, list(E = eT[i]))
  else if (drop.1 && mT[i] == -1) 
    substitute(-10^E, list(E = eT[i]))
  else substitute(A %.% 10^E, list(A = mT[i], E = eT[i]))
  do.call("expression", ss)
}




###############################################################################
#### Calculation of Trel for no coupling of off-rates
###############################################################################

TrelNO <- function(otot,param=parmsNO){
  
  #### Parameters
  k1 =param['kOpT']; D1 = param['KdOT']; k2 = param['kOTpE']; Et <- param['Et']
  D2 = param['KdOTE']; vt = param['vprod']; k4 = param['vdegrad']
  alpha=param['alpha']; kE=param['kcleav']; k3=param['kC']
  
  
  tmp <- function(Ot){
    
    #### Polynomial oefficients
    n4 <- 3*k1*k3*D2*k2*kE^3+3*k1*k4*k3*kE*D2^2*k2^2+3*k1*k4*k3*kE^2*D2*k2+k1*k4*k3*D2^3*k2^3+3*k1*k3*D2^2*k2^2*kE^2+k1*kE^2*D2^3*k2^3+2*k1*kE^3*D2^2*k2^2+k1*kE^4*D2*k2+k1*k3*kE^4+k1*k4*kE*D2^3*k2^3+2*k1*k4*kE^2*D2^2*k2^2+k1*k4*kE^3*D2*k2+k1*k4*k3*kE^3+k1*k3*D2^3*k2^3*kE
    
    n3 <- -4*k1*k4*k3*kE*D2^2*k2^2*Et-2*k1*k4*k3*kE^2*k2*Et*D2-2*k1*k4*k3*Ot*D2^2*k2^2*kE-k1*k4*k3*Ot*D2*k2*kE^2-k4^2*k3*D2^3*k2^3-k1*vt*k3*kE^3-3*D2^2*k4^2*k3*k2*k1*kE-3*k1*vt*k3*kE*D2^2*k2^2-3*k1*vt*k3*kE^2*D2*k2-2*k1*k4*kE*D2^3*k2^3*Et-2*k1*k4*kE^2*D2^2*k2^2*Et-k1*k3*Ot*D2^3*k2^3*kE-2*k1*k3*Ot*D2^2*k2^2*kE^2-k1*k3*Ot*D2*k2*kE^3-k1*k4*k3*Ot*D2^3*k2^3-2*k1*k3*D2^3*k2^3*kE*Et-4*k1*k3*D2^2*k2^2*kE^2*Et-2*k1*k4*k3*D2^3*k2^3*Et-2*k1*k3*kE^3*k2*Et*D2-6*D2^3*k1*kE*k2^2*k3*k4-5*D2^2*k1*k3*k2*kE^2*k4-D2*k1*k4*k3*kE^3-D2^4*k1*k3*k2^3*kE-D2^2*k1*k3*k2*kE^3-2*D2^4*k4*k3*k2^3*k1-2*D2^3*k4^2*k3*k2^2*k1-k3*k4*kE*D2^3*k2^3-2*k3*k4*kE^2*D2^2*k2^2-k3*k4*kE^3*D2*k2-2*k4^2*k3*D2^2*k2^2*kE-k4^2*k3*D2*k2*kE^2-D2^4*k1*k4*kE*k2^3-D2^3*k1*k4^2*kE*k2^2-2*D2^3*k1*k4*kE^2*k2^2-D2^2*k1*k4^2*kE^2*k2-D2^2*k1*k4*kE^3*k2-k1*vt*kE*D2^3*k2^3-2*k1*vt*kE^2*D2^2*k2^2-k1*vt*kE^3*D2*k2-k1*vt*k3*D2^3*k2^3-2*k1*kE^2*D2^3*k2^3*Et-2*k1*kE^3*D2^2*k2^2*Et-2*D2^3*k1*kE^2*k2^2*k3-D2*k1*k4^2*k3*kE^2
    
    n2 <- 2*k4^2*k3*D2^3*k2^3*Et+2*k1*k4*k3*Ot*D2^2*k2^2*Et*kE+2*k1*vt*k3*Ot*D2^2*k2^2*kE+k1*vt*k3*Ot*D2*k2*kE^2+D2^2*k1*k4^2*k3*Ot*k2*kE+D2^2*k1*k4*k3*Ot*k2*kE^2+D2^2*vt*k3*k2*k1*kE*k4+2*D2^3*k1*D1*k4*k3*k2^2*kE+D2^2*k1*D1*k4^2*k3*k2*kE+D2^2*k1*D1*k4*k3*k2*kE^2+4*D2^3*k1*k4*k3*kE*k2^2*Et+D2^2*k1*k4^2*k3*kE*k2*Et+D2^2*k1*k4*k3*kE^2*k2*Et+2*D2^3*k1*k4*k3*Ot*k2^2*kE+2*k1*vt*k3*kE^2*k2*Et*D2+2*k1*k3*Ot*D2^3*k2^3*kE*Et+2*k1*k3*Ot*D2^2*k2^2*kE^2*Et+k1*k4*k3*kE*k2^2*Et^2*D2^2+D2^4*k4^2*k3*k2^3+D2^3*k4^3*k3*k2^2+2*k3*k4*kE*D2^3*k2^3*Et+D2^3*k1*k4*k3*k2*kE^2+D2^4*k1*k4*k3*Ot*k2^3+D2^3*k1*k4^2*k3*Ot*k2^2+D2^4*k1*k3*k2^3*kE*Et+D2^3*k1*k3*k2^2*kE^2*Et+2*D2^4*k4*k3*k2^3*k1*Et+2*D2^3*k4^2*k3*k2^2*k1*Et+D2^4*k1*D1*k4*k3*k2^3+D2^3*k1*D1*k4^2*k3*k2^2+D2^2*vt*k3*k2*k1*kE^2+D2^4*k1*k4*kE*k2^3*Et+D2^3*k1*k4^2*kE*k2^2*Et+D2^3*k1*k4*kE^2*k2^2*Et+2*D2^4*k1*k4*k3*k2^2*kE+2*k1*vt*kE*D2^3*k2^3*Et+2*k1*vt*kE^2*D2^2*k2^2*Et+k1*k4*kE*D2^3*k2^3*Et^2+k1*vt*k3*Ot*D2^3*k2^3+k1*k3*D2^3*k2^3*kE*Et^2+k1*k4*k3*D2^3*k2^3*Et^2+k1*k3*kE^2*k2^2*Et^2*D2^2+2*D2^3*vt*k3*k2^2*k1*kE+2*k3*k4*kE^2*D2^2*k2^2*Et+2*D2^3*k1*k4^2*k3*k2*kE+D2^3*vt*k3*k2^2*k1*k4+D2^5*k1*k4*k3*k2^3+2*D2^4*k1*k4^2*k3*k2^2+D2^3*k1*k4^3*k3*k2+D2^4*vt*k3*k2^3*k1+2*D2^3*k4^2*k3*k2^2*kE+D2^2*k4^3*k3*k2*kE+D2^2*k4^2*k3*k2*kE^2+k1*kE^2*D2^3*k2^3*Et^2+2*vt*k3*D2^3*k2^3*k1*Et+2*k4^2*k3*D2^2*k2^2*Et*kE+4*vt*k3*D2^2*k2^2*k1*Et*kE+2*k1*k4*k3*Ot*D2^3*k2^3*Et
    
    n1 <- -D2^3*k1*k4*k3*Ot*k2^2*Et*kE-D2^3*k1*D1*k4*k3*k2^2*Et*kE-D2^4*k1*k4*k3*Ot*k2^3*Et-D2^3*vt*k3*k2^2*k1*Et*k4-D2^3*vt*k3*k2^2*k1*Et*kE-D2^4*k4^2*k3*k2^3*Et-k1*vt*k3*D2^3*k2^3*Et^2-D2^3*k1*k4^2*k3*Ot*k2^2*Et-k1*k3*Ot*D2^3*k2^3*kE*Et^2-k1*vt*kE*D2^3*k2^3*Et^2-D2^3*k4^3*k3*k2^2*Et-k1*k4*k3*Ot*D2^3*k2^3*Et^2-k4^2*k3*D2^3*k2^3*Et^2-D2^3*k1*D1*k4^2*k3*k2^2*Et-k1*vt*k3*kE*k2^2*Et^2*D2^2-2*k1*vt*k3*Ot*D2^2*k2^2*kE*Et-D2^4*k1*D1*k4*k3*k2^3*Et-k3*k4*kE*D2^3*k2^3*Et^2-D2^4*vt*k3*k2^3*k1*Et-D2^3*k4^2*k3*k2^2*Et*kE-2*k1*vt*k3*Ot*D2^3*k2^3*Et
    
    n0 <- k1*vt*k3*Ot*D2^3*k2^3*Et^2
    
    ROOT <- polyroot(c(n0,n1,n2,n3,n4))
    #### steady-state concentrations
    TrelSS <- -k4*(k1*vt*k3*Ot*D2^3*k2^3*Et^2+ROOT*D2^3*k1*D1*k4*k3*k2^2*Et*kE-2*ROOT*k1*vt*k3*Ot*D2^2*k2^2*kE*Et+ROOT^3*D2^4*k1*k2^3*kE^2+2*ROOT^3*D2^3*k1*k2^2*kE^3+ROOT^3*D2^3*k2^3*k3*kE^2+ROOT^3*D2^2*k1*k2*kE^4+2*ROOT^3*D2^2*k2^2*k3*kE^3+ROOT^3*D2*k1*k3*kE^4+ROOT^3*D2*k2*k3*kE^4-ROOT^3*k1*vt*k3*kE^3+4*ROOT^2*vt*k3*D2^2*k2^2*k1*Et*kE-ROOT*D2^3*vt*k3*k2^2*k1*Et*k4-ROOT*D2^3*vt*k3*k2^2*k1*Et*kE+D2^3*k2^2*Et*ROOT*Ot*k1*k3*kE^2+D2^3*k2^2*Et*ROOT*D1*k1*k3*kE^2+D2^4*k2^3*Et*ROOT*Ot*k1*k3*kE+D2^4*k2^3*Et*ROOT*D1*k1*k3*kE-ROOT*k1*vt*k3*kE*k2^2*Et^2*D2^2+ROOT^2*k1*vt*k3*Ot*D2*k2*kE^2-2*ROOT*k1*vt*k3*Ot*D2^3*k2^3*Et-ROOT^2*D2^2*k1*k4*k3*Ot*k2*kE^2+ROOT^2*D2^2*vt*k3*k2*k1*kE*k4+2*ROOT^2*k1*vt*k3*Ot*D2^2*k2^2*kE-ROOT^2*D2^3*k1*D1*k4*k3*k2^2*kE-ROOT^2*D2^2*k1*D1*k4*k3*k2*kE^2-ROOT^2*D2^3*k1*k4*k3*kE*k2^2*Et-ROOT^2*D2^2*k1*k4*k3*kE^2*k2*Et-ROOT^2*D2^3*k1*k4*k3*Ot*k2^2*kE+2*ROOT^2*k1*vt*k3*kE^2*k2*Et*D2-2*D2^2*k2^2*ROOT^2*Et*k3*kE^3-2*D2^4*k2^2*ROOT^2*k1*k3*kE^2-D2^5*k2^3*ROOT^2*k1*k3*kE+D2^3*k2^3*Et^2*ROOT*k3*kE^2+ROOT^2*D2^4*vt*k3*k2^3*k1-ROOT^2*D2^3*k4^2*k3*k2^2*kE-ROOT^2*D2^2*k4^2*k3*k2*kE^2+ROOT^3*D2*k1*k4*k3*kE^3+ROOT^3*D2^4*k1*k3*k2^3*kE+3*ROOT^3*D2^2*k1*k3*k2*kE^3+ROOT^3*k3*k4*kE*D2^3*k2^3+2*ROOT^3*k3*k4*kE^2*D2^2*k2^2+ROOT^3*k3*k4*kE^3*D2*k2+ROOT^3*D2^3*k1*k4*kE^2*k2^2+ROOT^3*D2^2*k1*k4*kE^3*k2-ROOT^3*k1*vt*kE*D2^3*k2^3-2*ROOT^3*k1*vt*kE^2*D2^2*k2^2-ROOT^3*k1*vt*kE^3*D2*k2-ROOT^3*k1*vt*k3*D2^3*k2^3+3*ROOT^3*D2^3*k1*kE^2*k2^2*k3-2*D2^3*k2^3*ROOT^2*Et*kE^2*k3-D2^4*k2^3*ROOT^2*k3*k4*kE-2*D2^3*k2^2*ROOT^2*k3*kE^2*k4-D2^3*k2*ROOT^2*k1*k3*kE^3-D2^2*k2*ROOT^2*k3*k4*kE^3-D2^4*k2^3*ROOT^2*Et*k1*kE^2-D2^3*k2^2*ROOT^2*Et*k1*kE^3-3*ROOT^3*k1*vt*k3*kE*D2^2*k2^2-3*ROOT^3*k1*vt*k3*kE^2*D2*k2-2*D2^3*k2^2*ROOT^2*Ot*k1*k3*kE^2-D2^2*k2*ROOT^2*D1*k1*k3*kE^3-D2^2*k2*ROOT^2*Et*k1*k3*kE^3-D2^2*k2*ROOT^2*Ot*k1*k3*kE^3+D2^3*k2^2*Et*ROOT*k4*k3*kE^2+D2^4*k2^3*Et*ROOT*k3*k4*kE-D2^4*k2^3*ROOT^2*k3*k1*D1*kE-2*D2^3*k2^2*ROOT^2*k3*kE^2*k1*D1-D2^4*k2^3*ROOT^2*Ot*k1*k3*kE+2*ROOT^2*D2^3*vt*k3*k2^2*k1*kE-2*ROOT^2*k3*k4*kE^2*D2^2*k2^2*Et-ROOT^2*D2^3*k1*k4^2*k3*k2*kE+ROOT^2*D2^3*vt*k3*k2^2*k1*k4+2*ROOT^2*vt*k3*D2^3*k2^3*k1*Et-ROOT*k1*vt*k3*D2^3*k2^3*Et^2-ROOT*k1*vt*kE*D2^3*k2^3*Et^2+ROOT*k3*k4*kE*D2^3*k2^3*Et^2-ROOT*D2^4*vt*k3*k2^3*k1*Et+ROOT*D2^3*k4^2*k3*k2^2*Et*kE+ROOT^3*D2^3*k1*kE*k2^2*k3*k4+2*ROOT^3*D2^2*k1*k3*k2*kE^2*k4-2*ROOT^2*k3*k4*kE*D2^3*k2^3*Et-2*ROOT^2*D2^3*k1*k4*k3*k2*kE^2-ROOT^2*D2^4*k1*k3*k2^3*kE*Et-2*ROOT^2*D2^3*k1*k3*k2^2*kE^2*Et+ROOT^2*D2^2*vt*k3*k2*k1*kE^2-ROOT^2*D2^3*k1*k4*kE^2*k2^2*Et-2*ROOT^2*D2^4*k1*k4*k3*k2^2*kE+2*ROOT^2*k1*vt*kE*D2^3*k2^3*Et+2*ROOT^2*k1*vt*kE^2*D2^2*k2^2*Et+ROOT^2*k1*vt*k3*Ot*D2^3*k2^3+ROOT*D2^3*k1*k4*k3*Ot*k2^2*Et*kE)/((kE+k4)*(kE*ROOT*D2^2*k2^2*Et-kE*ROOT^2*D2^2*k2^2-kE^2*ROOT^2*D2*k2-k3*Ot*D2^2*k2^2*Et+k3*Ot*D2^2*k2^2*ROOT+k3*Ot*D2*k2*kE*ROOT+k3*ROOT*D2^3*k2^2+k3*ROOT*D2^2*k2*k4+k3*ROOT*D2^2*k2*kE+k3*ROOT*D2^2*k2^2*Et-k3*ROOT^2*D2^2*k2^2-2*k3*ROOT^2*D2*k2*kE+k3*kE*ROOT*k2*Et*D2-k3*kE^2*ROOT^2)*vt*k1*(k2*Et*D2-D2*k2*ROOT-kE*ROOT))
    OTSS <- D2*ROOT*(D2*k2+k4+kE)/(k2*Et*D2-D2*k2*ROOT-kE*ROOT)
    ESS <- (k2*Et*D2-D2*k2*ROOT-kE*ROOT)/(D2*k2)
    OTESS <- ROOT
    OCESS <- kE*ROOT/(D2*k2)
    OSS <- (-k3*ROOT*D2^3*k2^2+k2*((Et-ROOT)*((-kE-k3)*ROOT+k3*Ot)*k2-k3*ROOT*(kE+k4))*D2^2-k2*((-kE-2*k3)*ROOT+k3*(Ot+Et))*kE*ROOT*D2+k3*kE^2*ROOT^2)/(k2*(k2*(Et-ROOT)*D2-kE*ROOT)*k3*D2)
    OCSS <- kE*ROOT/k3
    
    Re(TrelSS[ Re(TrelSS) <= 1 & Re(TrelSS) >= 0 & Im(TrelSS) < 1e-5 & 
                 Re(OSS) >=0 & Re(OTSS) >= 0 & 
                 Re(OTESS) >=0 & Re(OCESS) >=0 & Re(OCSS) >= 0 & Re(ESS) >=0] )
  }
  
  sapply(otot,tmp)
}

###############################################################################
#### Calculation of Trel for stochastic simulation of the ASOmodel
###############################################################################

Trelstoc <- function(Ot,km1=0.06){
  parms1['km1'] <- km1; parms1['k3'] <- km1/0.6
  x0['O'] <- Ot
  Gillespie <- ssa(x0=x0,a=a,nu=nu,parms = parms1,tf=2E2,method = "ETL")
  data <- rowSums(Gillespie$data[,2:4])/x0['Tt']
  Tmean <- mean(data[200:nrow(Gillespie$data)])
  Tsd <- sd(data[200:nrow(Gillespie$data)])
  return(list(Trel=Gillespie$data,Tstat=c('O'=x0['O'],TrelM=Tmean,TrelSD=Tsd)))
}
