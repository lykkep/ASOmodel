###############################################################################
#### Calculation of Trel in steday-state
###############################################################################
#' @title
#' Calculates the relative Target total target concentration at steady-state
#' 
#' @description
#' Use the algebraic solution to the system to the total target concentration at steady state relative to the target 
#' concentration before the addition of oligonucleotide. 
#' @param Ot Total concentration of oligonucleotide added to the system in nM
#' @param parms list of parameters
#' 
#' @return returns the relative total target concentration
#' @examples 
#' Trel(0.1,c(Et = 1,KdOT = 0.3,kOpT = 0.2,KdOTE = 70,
#' kOTpE = 5,vprod = 0.2,vdegrad = 0.04,alpha=0.1,kcleav = 8))
#' @export
Trel <- function(Ot,param=parms){  #
  
  #### Parameters
  k1 =param['kOpT']; D1 = param['KdOT']; k2 = param['kOTpE']; Et <- param['Et']
  D2 = param['KdOTE']; vt = param['vprod']; k4 = param['vdegrad']
  alpha=param['alpha']; kE=param['kcleav'];
  
  tmp <- function(Otot){

  #### Polynomial oefficients
  n4 <- kE^2*alpha*D2^3*k2^3+2*kE^3*alpha*D2^2*k2^2+kE^4*alpha*D2*k2+k4*k1*D1*kE^3+k4*kE*alpha*D2^3*k2^3+2*k4*kE^2*alpha*D2^2*k2^2+k4*kE^3*alpha*D2*k2+k4*k1*D1*D2^3*k2^3+k1*D1*D2^3*k2^3*kE+3*k1*D1*D2^2*k2^2*kE^2+3*k1*D1*D2*k2*kE^3+k1*D1*kE^4+3*D2^2*k2^2*k1*D1*k4*kE+3*k4*k1*D1*D2*k2*kE^2

  n0 <- vt*k1*D1*Otot*D2^3*k2^3*Et^2
  
  n1 <- -2*vt*k1*D1*Otot*D2^3*k2^3*Et-D2^4*k4^2*k2^3*D1*Et-D2^3*k4^3*k2^2*D1*Et-vt*k1*D1*D2^3*k2^3*Et^2-D2^3*D1^2*Et*k1*k2^2*k4*kE-D2^3*k1*D1*k2^2*vt*Et*k4-D2^3*D1^2*Et*k1*k2^2*k4^2-vt*kE*alpha*D2^3*k2^3*Et^2-D2^3*k4^2*k2^2*k1*D1*Otot*Et-D2^4*k4*k2^3*k1*D1*Otot*Et-D2^3*k1*D1*k2^2*vt*Et*kE-D2^4*k1*D1*k2^3*vt*Et-k4*D1*kE*D2^3*k2^3*Et^2-vt*k1*D1*kE*k2^2*Et^2*D2^2-k4^2*k2^3*D1*D2^3*Et^2-k4*k1*D1*Otot*D2^3*k2^3*Et^2-k1*D1*Otot*D2^3*k2^3*kE*Et^2-2*vt*k1*D1*Otot*D2^2*k2^2*Et*kE-D2^4*D1^2*Et*k1*k2^3*k4-D2^3*k4*k2^2*k1*D1*Otot*Et*kE-D2^3*k4^2*k2^2*D1*Et*kE
  
  n2 <- D1^2*D2^2*k1*k2*k4*kE^2+2*D1^2*D2^3*k1*k2^2*k4*kE+kE^2*alpha*D2^3*k2^3*Et^2+D2^4*k1*D1*k2^3*vt+2*D2^3*k4^2*k2^2*D1*kE+D2^2*k4^3*k2*D1*kE+D2^2*k4^2*k2*D1*kE^2+D2^5*k4*k1*D1*k2^3+2*D2^4*k4^2*k1*D1*k2^2+D2^3*k4^3*k1*D1*k2+D2^4*k4*k2^3*D1^2*k1+2*k1*D1*D2^3*k2^3*vt*Et+2*k4^2*D2^2*k2^2*D1*Et*kE+4*k1*D1*D2^2*k2^2*vt*Et*kE+2*k4*D2^3*k2^3*k1*D1*Otot*Et+D1^2*D2^3*k1*k2^2*k4^2+2*k4*D2^2*k2^2*k1*D1*Otot*Et*kE+2*k4^2*D2^3*k2^3*D1*Et+D2^4*k4^2*k2^3*D1+D2^3*k4^3*k2^2*D1+2*vt*k1*D1*kE^2*k2*Et*D2+2*vt*k1*D1*Otot*D2^2*k2^2*kE+2*k1*D1*Otot*D2^3*k2^3*kE*Et+2*k1*D1*Otot*D2^2*k2^2*kE^2*Et+vt*k1*D1*Otot*D2*k2*kE^2+2*D2^3*k4*k2^2*k1*D1*Otot*kE+k4*k1*D1*kE*k2^2*Et^2*D2^2+4*D2^3*k4*k1*D1*k2^2*Et*kE+D2^2*k1*D1*k2*vt*kE*k4+D2^2*k4^2*k1*D1*kE*k2*Et+D2^2*k4*k1*D1*kE^2*k2*Et+D2^2*k4^2*k2*k1*D1*Otot*kE+D2^2*k4*k2*k1*D1*Otot*kE^2+k1*D1*kE^2*k2^2*Et^2*D2^2+2*k4*D1*kE*D2^3*k2^3*Et+2*k4*D1*kE^2*D2^2*k2^2*Et+D2^3*k4*k1*D1*k2*kE^2+vt*k1*D1*Otot*D2^3*k2^3+k4*kE*alpha*D2^3*k2^3*Et^2+D2^3*k1*D1*k2^2*vt*k4+D2^4*k4*k2^3*k1*D1*Otot+D2^3*k4^2*k2^2*k1*D1*Otot+2*D2^4*k4*k1*D1*k2^3*Et+2*D2^3*k4^2*k1*D1*k2^2*Et+D2^4*k1*D1*k2^3*kE*Et+D2^3*k1*D1*k2^2*kE^2*Et+2*D2^4*k4*k1*D1*k2^2*kE+2*D2^3*k4^2*k1*D1*k2*kE+2*vt*kE*alpha*D2^3*k2^3*Et+2*vt*kE^2*alpha*D2^2*k2^2*Et+k4*k1*D1*D2^3*k2^3*Et^2+k1*D1*D2^3*k2^3*kE*Et^2+D2^2*k4^2*k2*D1^2*k1*kE+2*D2^3*k1*D1*k2^2*vt*kE+D2^2*k1*D1*k2*vt*kE^2+D2^4*k4*kE*alpha*k2^3*Et+D2^3*k4^2*kE*alpha*k2^2*Et+D2^3*k4*kE^2*alpha*k2^2*Et
    
  n3 <- -D2^4*k4*kE*alpha*k2^3-D2^3*k4^2*kE*alpha*k2^2-2*D2^3*k4*kE^2*alpha*k2^2-D2^2*k4^2*kE^2*alpha*k2-D2^2*k4*kE^3*alpha*k2-2*D2^4*k4*k1*D1*k2^3-2*D2^3*k4^2*k1*D1*k2^2-D2^4*k1*D1*k2^3*kE-k4^2*k2*D1*D2*kE^2-2*D2^3*k1*D1*k2^2*kE^2-D2^2*k1*D1*k2*kE^3-D2*k4*k1*D1*kE^3-D2*k4^2*k1*D1*kE^2-vt*kE*alpha*D2^3*k2^3-2*vt*kE^2*alpha*D2^2*k2^2-vt*kE^3*alpha*D2*k2-k4*D1*kE*D2^3*k2^3-2*k4*D1*kE^2*D2^2*k2^2-k4*D1*kE^3*D2*k2-vt*k1*D1*D2^3*k2^3-2*kE^2*alpha*D2^3*k2^3*Et-2*kE^3*alpha*D2^2*k2^2*Et-2*D2^2*k2^2*D1*k4^2*kE-k4^2*k2^3*D1*D2^3-vt*k1*D1*kE^3-4*k4*k1*D1*D2^2*k2^2*Et*kE-2*k4*k1*D1*Otot*D2^2*k2^2*kE-k4*k1*D1*Otot*D2*k2*kE^2-2*k4*k1*D1*kE^2*k2*Et*D2-4*k1*D1*D2^2*k2^2*kE^2*Et-2*k1*D1*kE^3*k2*Et*D2-k4*k1*D1*Otot*D2^3*k2^3-k1*D1*Otot*D2^3*k2^3*kE-2*k1*D1*Otot*D2^2*k2^2*kE^2-k1*D1*Otot*D2*k2*kE^3-2*k4*kE*alpha*D2^3*k2^3*Et-3*D2^2*k4^2*k1*D1*kE*k2-5*D2^2*k4*k1*D1*kE^2*k2-2*k4*kE^2*alpha*D2^2*k2^2*Et-3*vt*k1*D1*kE*D2^2*k2^2-3*vt*k1*D1*kE^2*D2*k2-2*k4*k1*D1*D2^3*k2^3*Et-2*k1*D1*D2^3*k2^3*kE*Et-6*D2^3*k4*k1*D1*kE*k2^2
  
  ROOT <- polyroot(c(n0,n1,n2,n3,n4))
  #### [OTE] at steady-state
  TrelSS <- -k4*(-2*ROOT^2*k4*D1*kE^2*D2^2*k2^2*Et-2*ROOT^2*D2^3*k4*k1*D1*k2*kE^2+ROOT^2*vt*k1*D1*Otot*D2^3*k2^3+ROOT^2*D2^3*k1*D1*k2^2*vt*k4-ROOT^2*D2^4*k1*D1*k2^3*kE*Et-2*ROOT^2*D2^3*k1*D1*k2^2*kE^2*Et-2*ROOT^2*D2^4*k4*k1*D1*k2^2*kE-ROOT^2*D2^3*k4^2*k1*D1*k2*kE+2*ROOT^2*vt*kE*alpha*D2^3*k2^3*Et+2*ROOT^2*vt*kE^2*alpha*D2^2*k2^2*Et+2*ROOT^2*D2^3*k1*D1*k2^2*vt*kE+ROOT^2*D2^2*k1*D1*k2*vt*kE^2-ROOT^2*D2^3*k4*kE^2*alpha*k2^2*Et-ROOT*vt*k1*D1*D2^3*k2^3*Et^2-ROOT*vt*kE*alpha*D2^3*k2^3*Et^2-ROOT*D2^4*k1*D1*k2^3*vt*Et+ROOT*k4*D1*kE*D2^3*k2^3*Et^2+ROOT*D2^3*k4^2*k2^2*D1*Et*kE+2*ROOT^3*D2^2*k4*k1*D1*kE^2*k2-3*ROOT^3*vt*k1*D1*kE*D2^2*k2^2-3*ROOT^3*vt*k1*D1*kE^2*D2*k2+ROOT^3*D2^3*k4*k1*D1*kE*k2^2-ROOT^2*D1^2*D2^2*k1*k2*k4*kE^2-ROOT^2*D1^2*D2^3*k1*k2^2*k4*kE-D2^2*k2*ROOT^2*D1*Et*k1*kE^3-D2^2*k2*ROOT^2*D1*Otot*k1*kE^3-D2^4*k2^3*ROOT^2*D1*Otot*k1*kE-2*D2^3*k2^2*ROOT^2*D1*Otot*k1*kE^2+D2^3*k2^2*Et*ROOT*D1*k4*kE^2+D2^4*k2^3*Et*ROOT*D1*k4*kE+D2^3*k2^2*Et*ROOT*D1^2*k1*kE^2+D2^4*k2^3*Et*ROOT*D1^2*k1*kE+2*ROOT^2*k1*D1*D2^3*k2^3*vt*Et-2*ROOT^2*k4*D1*kE*D2^3*k2^3*Et+ROOT^3*D2^4*alpha*k2^3*kE^2+ROOT^3*D1*D2^3*k2^3*kE^2+2*ROOT^3*D2^3*alpha*k2^2*kE^3+2*ROOT^3*D1*D2^2*k2^2*kE^3+ROOT^3*D2^2*alpha*k2*kE^4+ROOT^3*D1*D2*k1*kE^4+ROOT^3*D1*D2*k2*kE^4-ROOT^3*vt*k1*D1*kE^3+vt*k1*D1*Otot*D2^3*k2^3*Et^2+ROOT^2*D2^2*k1*D1*k2*vt*kE*k4-ROOT*D2^3*k1*D1*k2^2*vt*Et*kE+ROOT^2*vt*k1*D1*Otot*D2*k2*kE^2+4*ROOT^2*k1*D1*D2^2*k2^2*vt*Et*kE+2*ROOT^2*vt*k1*D1*kE^2*k2*Et*D2-ROOT^2*D2^2*k4*k2*k1*D1*Otot*kE^2+2*ROOT^2*vt*k1*D1*Otot*D2^2*k2^2*kE-ROOT^2*D2^3*k4*k1*D1*k2^2*Et*kE-ROOT*D2^3*k1*D1*k2^2*vt*Et*k4-ROOT^2*D2^3*k4*k2^2*k1*D1*Otot*kE+D2^4*k2^3*Et*ROOT*D1*Otot*k1*kE+ROOT*D2^3*D1^2*Et*k1*k2^2*k4*kE-ROOT*vt*k1*D1*kE*k2^2*Et^2*D2^2+D2^3*k2^2*Et*ROOT*D1*Otot*k1*kE^2-ROOT^2*D2^2*k4*k1*D1*kE^2*k2*Et-2*ROOT*vt*k1*D1*Otot*D2^3*k2^3*Et-2*ROOT*vt*k1*D1*Otot*D2^2*k2^2*Et*kE+ROOT*D2^3*k4*k2^2*k1*D1*Otot*Et*kE-D2^3*k2*ROOT^2*D1*k1*kE^3-D2^5*k2^3*ROOT^2*k1*D1*kE-2*D2^4*k2^2*ROOT^2*k1*D1*kE^2-D2^4*k2^3*ROOT^2*D1^2*k1*kE-2*D2^3*k2^2*ROOT^2*D1^2*k1*kE^2+D2^3*k2^3*Et^2*ROOT*D1*kE^2+ROOT^2*D2^4*k1*D1*k2^3*vt-ROOT^2*D2^3*k4^2*k2^2*D1*kE-ROOT^2*D2^2*k4^2*k2*D1*kE^2+ROOT^3*D2^3*k4*kE^2*alpha*k2^2+ROOT^3*D2^2*k4*kE^3*alpha*k2+ROOT^3*D2^4*k1*D1*k2^3*kE+3*ROOT^3*D2^3*k1*D1*k2^2*kE^2+3*ROOT^3*D2^2*k1*D1*k2*kE^3+ROOT^3*D2*k4*k1*D1*kE^3-ROOT^3*vt*kE*alpha*D2^3*k2^3-2*ROOT^3*vt*kE^2*alpha*D2^2*k2^2-ROOT^3*vt*kE^3*alpha*D2*k2+ROOT^3*k4*D1*kE*D2^3*k2^3+2*ROOT^3*k4*D1*kE^2*D2^2*k2^2+ROOT^3*k4*D1*kE^3*D2*k2-ROOT^3*vt*k1*D1*D2^3*k2^3-D2^4*k2^3*ROOT^2*Et*alpha*kE^2-D2^4*k2^3*ROOT^2*D1*k4*kE-2*D2^3*k2^3*ROOT^2*D1*Et*kE^2-D2^3*k2^2*ROOT^2*Et*alpha*kE^3-2*D2^3*k2^2*ROOT^2*D1*k4*kE^2-2*D2^2*k2^2*ROOT^2*D1*Et*kE^3-D2^2*k2*ROOT^2*k4*D1*kE^3-D2^2*k2*ROOT^2*D1^2*k1*kE^3)/((kE+k4)*(kE*ROOT*alpha*D2^2*k2^2*Et-kE*ROOT^2*alpha*D2^2*k2^2-kE^2*ROOT^2*alpha*D2*k2-k1*D1*Otot*D2^2*k2^2*Et+k1*D1*Otot*D2^2*k2^2*ROOT+k1*D1*Otot*D2*k2*kE*ROOT+k1*D1*ROOT*D2^3*k2^2+k4*k1*D1*D2^2*k2*ROOT+k1*D1*ROOT*D2^2*k2*kE+k1*D1*ROOT*D2^2*k2^2*Et-k1*D1*ROOT^2*D2^2*k2^2-2*k1*D1*ROOT^2*D2*k2*kE+k1*D1*kE*ROOT*k2*Et*D2-k1*D1*kE^2*ROOT^2)*vt*(k2*Et*D2-D2*k2*ROOT-kE*ROOT))
  OTSS <- D2*ROOT*(D2*k2+k4+kE)/(k2*Et*D2-D2*k2*ROOT-kE*ROOT)
  ESS <- (k2*Et*D2-D2*k2*ROOT-kE*ROOT)/(D2*k2)
  OTESS <- ROOT
  OCESS <- kE*ROOT/(D2*k2)
  OSS <- (-k1*D1*ROOT*D2^3*k2^2+(((-k1*D1-kE*alpha)*ROOT+k1*D1*Otot)*(Et-ROOT)*k2-k1*D1*ROOT*(kE+k4))*k2*D2^2-kE*((-kE*alpha-2*k1*D1)*ROOT+D1*k1*(Otot+Et))*ROOT*k2*D2+k1*D1*kE^2*ROOT^2)/(D1*(k2*(Et-ROOT)*D2-kE*ROOT)*D2*k2*k1)
  OCSS <- kE*ROOT*alpha/(k1*D1)
  
  Re(TrelSS[ Re(TrelSS) <= 1 & Re(TrelSS) >= 0 & Im(TrelSS) < 1e-5 & 
               Re(OSS) >=0 & Re(OTSS) >= 0 & 
               Re(OTESS) >=0 & Re(OCESS) >=0 & Re(OCSS) >= 0 & Re(ESS) >=0] )
  }

  sapply(Ot,tmp)
}


###############################################################################
#### Calculation of Trel for no coupling of off-rates
###############################################################################
#' @title
#' Calculates the relative Target total target concentration at steady-state
#' 
#' @description
#' Use the algebraic solution to the system to the total target concentration at steady state relative to the target 
#' concentration before the addition of oligonucleotide. 
#' @param Ot Total concentration of oligonucleotide added to the system in nM
#' @param parms list of parameters
#' 
#' @return returns the relative total target concentration
#' @examples 
#' Trel(0.1,c(Et = 1,KdOT = 0.3,kOpT = 0.2,KdOTE = 70,
#' kOTpE = 5,vprod = 0.2,vdegrad = 0.04,alpha=0.1,kcleav = 8))
#' @export
TrelNO <- function(Ot,param=parmsNO){
  
  #### Parameters
  k1 =param['kOpT']; D1 = param['KdOT']; k2 = param['kOTpE']; Et <- param['Et']
  D2 = param['KdOTE']; vt = param['vprod']; k4 = param['vdegrad']
  alpha=param['alpha']; kE=param['kcleav']; k3=param['kC']
  
  
  tmp <- function(Otot){
    
    #### Polynomial oefficients
    n4 <- 3*k1*k3*D2*k2*kE^3+3*k1*k4*k3*kE*D2^2*k2^2+3*k1*k4*k3*kE^2*D2*k2+k1*k4*k3*D2^3*k2^3+3*k1*k3*D2^2*k2^2*kE^2+k1*kE^2*D2^3*k2^3+2*k1*kE^3*D2^2*k2^2+k1*kE^4*D2*k2+k1*k3*kE^4+k1*k4*kE*D2^3*k2^3+2*k1*k4*kE^2*D2^2*k2^2+k1*k4*kE^3*D2*k2+k1*k4*k3*kE^3+k1*k3*D2^3*k2^3*kE
    
    n3 <- -4*k1*k4*k3*kE*D2^2*k2^2*Et-2*k1*k4*k3*kE^2*k2*Et*D2-2*k1*k4*k3*Otot*D2^2*k2^2*kE-k1*k4*k3*Otot*D2*k2*kE^2-k4^2*k3*D2^3*k2^3-k1*vt*k3*kE^3-3*D2^2*k4^2*k3*k2*k1*kE-3*k1*vt*k3*kE*D2^2*k2^2-3*k1*vt*k3*kE^2*D2*k2-2*k1*k4*kE*D2^3*k2^3*Et-2*k1*k4*kE^2*D2^2*k2^2*Et-k1*k3*Otot*D2^3*k2^3*kE-2*k1*k3*Otot*D2^2*k2^2*kE^2-k1*k3*Otot*D2*k2*kE^3-k1*k4*k3*Otot*D2^3*k2^3-2*k1*k3*D2^3*k2^3*kE*Et-4*k1*k3*D2^2*k2^2*kE^2*Et-2*k1*k4*k3*D2^3*k2^3*Et-2*k1*k3*kE^3*k2*Et*D2-6*D2^3*k1*kE*k2^2*k3*k4-5*D2^2*k1*k3*k2*kE^2*k4-D2*k1*k4*k3*kE^3-D2^4*k1*k3*k2^3*kE-D2^2*k1*k3*k2*kE^3-2*D2^4*k4*k3*k2^3*k1-2*D2^3*k4^2*k3*k2^2*k1-k3*k4*kE*D2^3*k2^3-2*k3*k4*kE^2*D2^2*k2^2-k3*k4*kE^3*D2*k2-2*k4^2*k3*D2^2*k2^2*kE-k4^2*k3*D2*k2*kE^2-D2^4*k1*k4*kE*k2^3-D2^3*k1*k4^2*kE*k2^2-2*D2^3*k1*k4*kE^2*k2^2-D2^2*k1*k4^2*kE^2*k2-D2^2*k1*k4*kE^3*k2-k1*vt*kE*D2^3*k2^3-2*k1*vt*kE^2*D2^2*k2^2-k1*vt*kE^3*D2*k2-k1*vt*k3*D2^3*k2^3-2*k1*kE^2*D2^3*k2^3*Et-2*k1*kE^3*D2^2*k2^2*Et-2*D2^3*k1*kE^2*k2^2*k3-D2*k1*k4^2*k3*kE^2
    
    n2 <- 2*k4^2*k3*D2^3*k2^3*Et+2*k1*k4*k3*Otot*D2^2*k2^2*Et*kE+2*k1*vt*k3*Otot*D2^2*k2^2*kE+k1*vt*k3*Otot*D2*k2*kE^2+D2^2*k1*k4^2*k3*Otot*k2*kE+D2^2*k1*k4*k3*Otot*k2*kE^2+D2^2*vt*k3*k2*k1*kE*k4+2*D2^3*k1*D1*k4*k3*k2^2*kE+D2^2*k1*D1*k4^2*k3*k2*kE+D2^2*k1*D1*k4*k3*k2*kE^2+4*D2^3*k1*k4*k3*kE*k2^2*Et+D2^2*k1*k4^2*k3*kE*k2*Et+D2^2*k1*k4*k3*kE^2*k2*Et+2*D2^3*k1*k4*k3*Otot*k2^2*kE+2*k1*vt*k3*kE^2*k2*Et*D2+2*k1*k3*Otot*D2^3*k2^3*kE*Et+2*k1*k3*Otot*D2^2*k2^2*kE^2*Et+k1*k4*k3*kE*k2^2*Et^2*D2^2+D2^4*k4^2*k3*k2^3+D2^3*k4^3*k3*k2^2+2*k3*k4*kE*D2^3*k2^3*Et+D2^3*k1*k4*k3*k2*kE^2+D2^4*k1*k4*k3*Otot*k2^3+D2^3*k1*k4^2*k3*Otot*k2^2+D2^4*k1*k3*k2^3*kE*Et+D2^3*k1*k3*k2^2*kE^2*Et+2*D2^4*k4*k3*k2^3*k1*Et+2*D2^3*k4^2*k3*k2^2*k1*Et+D2^4*k1*D1*k4*k3*k2^3+D2^3*k1*D1*k4^2*k3*k2^2+D2^2*vt*k3*k2*k1*kE^2+D2^4*k1*k4*kE*k2^3*Et+D2^3*k1*k4^2*kE*k2^2*Et+D2^3*k1*k4*kE^2*k2^2*Et+2*D2^4*k1*k4*k3*k2^2*kE+2*k1*vt*kE*D2^3*k2^3*Et+2*k1*vt*kE^2*D2^2*k2^2*Et+k1*k4*kE*D2^3*k2^3*Et^2+k1*vt*k3*Otot*D2^3*k2^3+k1*k3*D2^3*k2^3*kE*Et^2+k1*k4*k3*D2^3*k2^3*Et^2+k1*k3*kE^2*k2^2*Et^2*D2^2+2*D2^3*vt*k3*k2^2*k1*kE+2*k3*k4*kE^2*D2^2*k2^2*Et+2*D2^3*k1*k4^2*k3*k2*kE+D2^3*vt*k3*k2^2*k1*k4+D2^5*k1*k4*k3*k2^3+2*D2^4*k1*k4^2*k3*k2^2+D2^3*k1*k4^3*k3*k2+D2^4*vt*k3*k2^3*k1+2*D2^3*k4^2*k3*k2^2*kE+D2^2*k4^3*k3*k2*kE+D2^2*k4^2*k3*k2*kE^2+k1*kE^2*D2^3*k2^3*Et^2+2*vt*k3*D2^3*k2^3*k1*Et+2*k4^2*k3*D2^2*k2^2*Et*kE+4*vt*k3*D2^2*k2^2*k1*Et*kE+2*k1*k4*k3*Otot*D2^3*k2^3*Et
    
    n1 <- -D2^3*k1*k4*k3*Otot*k2^2*Et*kE-D2^3*k1*D1*k4*k3*k2^2*Et*kE-D2^4*k1*k4*k3*Otot*k2^3*Et-D2^3*vt*k3*k2^2*k1*Et*k4-D2^3*vt*k3*k2^2*k1*Et*kE-D2^4*k4^2*k3*k2^3*Et-k1*vt*k3*D2^3*k2^3*Et^2-D2^3*k1*k4^2*k3*Otot*k2^2*Et-k1*k3*Otot*D2^3*k2^3*kE*Et^2-k1*vt*kE*D2^3*k2^3*Et^2-D2^3*k4^3*k3*k2^2*Et-k1*k4*k3*Otot*D2^3*k2^3*Et^2-k4^2*k3*D2^3*k2^3*Et^2-D2^3*k1*D1*k4^2*k3*k2^2*Et-k1*vt*k3*kE*k2^2*Et^2*D2^2-2*k1*vt*k3*Otot*D2^2*k2^2*kE*Et-D2^4*k1*D1*k4*k3*k2^3*Et-k3*k4*kE*D2^3*k2^3*Et^2-D2^4*vt*k3*k2^3*k1*Et-D2^3*k4^2*k3*k2^2*Et*kE-2*k1*vt*k3*Otot*D2^3*k2^3*Et
    
    n0 <- k1*vt*k3*Otot*D2^3*k2^3*Et^2
    
    ROOT <- polyroot(c(n0,n1,n2,n3,n4))
    #### steady-state concentrations
    TrelSS <- -k4*(k1*vt*k3*Otot*D2^3*k2^3*Et^2+ROOT*D2^3*k1*D1*k4*k3*k2^2*Et*kE-2*ROOT*k1*vt*k3*Otot*D2^2*k2^2*kE*Et+ROOT^3*D2^4*k1*k2^3*kE^2+2*ROOT^3*D2^3*k1*k2^2*kE^3+ROOT^3*D2^3*k2^3*k3*kE^2+ROOT^3*D2^2*k1*k2*kE^4+2*ROOT^3*D2^2*k2^2*k3*kE^3+ROOT^3*D2*k1*k3*kE^4+ROOT^3*D2*k2*k3*kE^4-ROOT^3*k1*vt*k3*kE^3+4*ROOT^2*vt*k3*D2^2*k2^2*k1*Et*kE-ROOT*D2^3*vt*k3*k2^2*k1*Et*k4-ROOT*D2^3*vt*k3*k2^2*k1*Et*kE+D2^3*k2^2*Et*ROOT*Otot*k1*k3*kE^2+D2^3*k2^2*Et*ROOT*D1*k1*k3*kE^2+D2^4*k2^3*Et*ROOT*Otot*k1*k3*kE+D2^4*k2^3*Et*ROOT*D1*k1*k3*kE-ROOT*k1*vt*k3*kE*k2^2*Et^2*D2^2+ROOT^2*k1*vt*k3*Otot*D2*k2*kE^2-2*ROOT*k1*vt*k3*Otot*D2^3*k2^3*Et-ROOT^2*D2^2*k1*k4*k3*Otot*k2*kE^2+ROOT^2*D2^2*vt*k3*k2*k1*kE*k4+2*ROOT^2*k1*vt*k3*Otot*D2^2*k2^2*kE-ROOT^2*D2^3*k1*D1*k4*k3*k2^2*kE-ROOT^2*D2^2*k1*D1*k4*k3*k2*kE^2-ROOT^2*D2^3*k1*k4*k3*kE*k2^2*Et-ROOT^2*D2^2*k1*k4*k3*kE^2*k2*Et-ROOT^2*D2^3*k1*k4*k3*Otot*k2^2*kE+2*ROOT^2*k1*vt*k3*kE^2*k2*Et*D2-2*D2^2*k2^2*ROOT^2*Et*k3*kE^3-2*D2^4*k2^2*ROOT^2*k1*k3*kE^2-D2^5*k2^3*ROOT^2*k1*k3*kE+D2^3*k2^3*Et^2*ROOT*k3*kE^2+ROOT^2*D2^4*vt*k3*k2^3*k1-ROOT^2*D2^3*k4^2*k3*k2^2*kE-ROOT^2*D2^2*k4^2*k3*k2*kE^2+ROOT^3*D2*k1*k4*k3*kE^3+ROOT^3*D2^4*k1*k3*k2^3*kE+3*ROOT^3*D2^2*k1*k3*k2*kE^3+ROOT^3*k3*k4*kE*D2^3*k2^3+2*ROOT^3*k3*k4*kE^2*D2^2*k2^2+ROOT^3*k3*k4*kE^3*D2*k2+ROOT^3*D2^3*k1*k4*kE^2*k2^2+ROOT^3*D2^2*k1*k4*kE^3*k2-ROOT^3*k1*vt*kE*D2^3*k2^3-2*ROOT^3*k1*vt*kE^2*D2^2*k2^2-ROOT^3*k1*vt*kE^3*D2*k2-ROOT^3*k1*vt*k3*D2^3*k2^3+3*ROOT^3*D2^3*k1*kE^2*k2^2*k3-2*D2^3*k2^3*ROOT^2*Et*kE^2*k3-D2^4*k2^3*ROOT^2*k3*k4*kE-2*D2^3*k2^2*ROOT^2*k3*kE^2*k4-D2^3*k2*ROOT^2*k1*k3*kE^3-D2^2*k2*ROOT^2*k3*k4*kE^3-D2^4*k2^3*ROOT^2*Et*k1*kE^2-D2^3*k2^2*ROOT^2*Et*k1*kE^3-3*ROOT^3*k1*vt*k3*kE*D2^2*k2^2-3*ROOT^3*k1*vt*k3*kE^2*D2*k2-2*D2^3*k2^2*ROOT^2*Otot*k1*k3*kE^2-D2^2*k2*ROOT^2*D1*k1*k3*kE^3-D2^2*k2*ROOT^2*Et*k1*k3*kE^3-D2^2*k2*ROOT^2*Otot*k1*k3*kE^3+D2^3*k2^2*Et*ROOT*k4*k3*kE^2+D2^4*k2^3*Et*ROOT*k3*k4*kE-D2^4*k2^3*ROOT^2*k3*k1*D1*kE-2*D2^3*k2^2*ROOT^2*k3*kE^2*k1*D1-D2^4*k2^3*ROOT^2*Otot*k1*k3*kE+2*ROOT^2*D2^3*vt*k3*k2^2*k1*kE-2*ROOT^2*k3*k4*kE^2*D2^2*k2^2*Et-ROOT^2*D2^3*k1*k4^2*k3*k2*kE+ROOT^2*D2^3*vt*k3*k2^2*k1*k4+2*ROOT^2*vt*k3*D2^3*k2^3*k1*Et-ROOT*k1*vt*k3*D2^3*k2^3*Et^2-ROOT*k1*vt*kE*D2^3*k2^3*Et^2+ROOT*k3*k4*kE*D2^3*k2^3*Et^2-ROOT*D2^4*vt*k3*k2^3*k1*Et+ROOT*D2^3*k4^2*k3*k2^2*Et*kE+ROOT^3*D2^3*k1*kE*k2^2*k3*k4+2*ROOT^3*D2^2*k1*k3*k2*kE^2*k4-2*ROOT^2*k3*k4*kE*D2^3*k2^3*Et-2*ROOT^2*D2^3*k1*k4*k3*k2*kE^2-ROOT^2*D2^4*k1*k3*k2^3*kE*Et-2*ROOT^2*D2^3*k1*k3*k2^2*kE^2*Et+ROOT^2*D2^2*vt*k3*k2*k1*kE^2-ROOT^2*D2^3*k1*k4*kE^2*k2^2*Et-2*ROOT^2*D2^4*k1*k4*k3*k2^2*kE+2*ROOT^2*k1*vt*kE*D2^3*k2^3*Et+2*ROOT^2*k1*vt*kE^2*D2^2*k2^2*Et+ROOT^2*k1*vt*k3*Otot*D2^3*k2^3+ROOT*D2^3*k1*k4*k3*Otot*k2^2*Et*kE)/((kE+k4)*(kE*ROOT*D2^2*k2^2*Et-kE*ROOT^2*D2^2*k2^2-kE^2*ROOT^2*D2*k2-k3*Otot*D2^2*k2^2*Et+k3*Otot*D2^2*k2^2*ROOT+k3*Otot*D2*k2*kE*ROOT+k3*ROOT*D2^3*k2^2+k3*ROOT*D2^2*k2*k4+k3*ROOT*D2^2*k2*kE+k3*ROOT*D2^2*k2^2*Et-k3*ROOT^2*D2^2*k2^2-2*k3*ROOT^2*D2*k2*kE+k3*kE*ROOT*k2*Et*D2-k3*kE^2*ROOT^2)*vt*k1*(k2*Et*D2-D2*k2*ROOT-kE*ROOT))
    OTSS <- D2*ROOT*(D2*k2+k4+kE)/(k2*Et*D2-D2*k2*ROOT-kE*ROOT)
    ESS <- (k2*Et*D2-D2*k2*ROOT-kE*ROOT)/(D2*k2)
    OTESS <- ROOT
    OCESS <- kE*ROOT/(D2*k2)
    OSS <- (-k3*ROOT*D2^3*k2^2+k2*((Et-ROOT)*((-kE-k3)*ROOT+k3*Otot)*k2-k3*ROOT*(kE+k4))*D2^2-k2*((-kE-2*k3)*ROOT+k3*(Otot+Et))*kE*ROOT*D2+k3*kE^2*ROOT^2)/(k2*(k2*(Et-ROOT)*D2-kE*ROOT)*k3*D2)
    OCSS <- kE*ROOT/k3
    
    Re(TrelSS[ Re(TrelSS) <= 1 & Re(TrelSS) >= 0 & Im(TrelSS) < 1e-5 & 
                 Re(OSS) >=0 & Re(OTSS) >= 0 & 
                 Re(OTESS) >=0 & Re(OCESS) >=0 & Re(OCSS) >= 0 & Re(ESS) >=0] )
  }
  
  sapply(Ot,tmp)
}

###############################################################################
#### Calculation of Trel for stochastic simulation of the ASOmodel
###############################################################################
#' @title Calculates the relative total target level for a stochastic simulation of the ASO model
#' 
#' @description
#' Uses a Gillespie algorithm (from the GillespieSSA package) to stochastically simulate the relative total target level.
#' 
#' @param Ot The total oligonucleotide concentration
#' @param kOT The dissociation rate between the oligonucleotide and the target
#' 
#' @return returns timeseries for the stochatic simluation ($data) and a vector ($Tstat) with the inital Ot, the mean of Trel and the standard deviation of Trel after the system has reached steady-state.
#' @examples 
#' parms1 <- c(kOpT = 2E-5,kOTpE =50E-5 ,vprod = 150,  vdegrad = 0.04,      
#' kcleav = 2, kOT =0.06, kOTE=2, kC = 0.1)
#' #Initital state vector
#' x0 <- c(Tt=parms1["vprod"]/parms1["vdegrad"],
#'              OT=0,OTE=0,E=1e3,O=1e5,OCE=0,OC=0)
#' names(x0) <- c('Tt','OT','OTE','E','O','OCE','OC')
#' #Propensity vector
#' a <-  c("vprod","kOpT*O*Tt","vdegrad*Tt","kOT*OT","kOTE*OTE","vdegrad*OT",
#'          "kOTpE*OT*E","vdegrad*OTE","kcleav*OTE","kC*OC","kOTE*OCE" )
#' #State-change matrix
#' nu <- matrix(0,7,length(a))
#' dimnames(nu) <- list(names(x0),a)
#' nu['Tt',c('vprod','kOT*OT')] <- 1
#' nu['Tt',c('kOpT*O*Tt','vdegrad*Tt')] <- -1 
#' nu['OT',c('kOpT*O*Tt','kOTE*OTE')] <- 1
#' nu['OT',c('kOT*OT','kOTpE*OT*E','vdegrad*OT')] <- -1
#' nu['OTE',c('kOTpE*OT*E')] <- 1
#' nu['OTE',c('kOTE*OTE','vdegrad*OTE','kcleav*OTE')] <- -1
#' nu['E',c('kOTE*OTE','vdegrad*OTE','kOTE*OCE')] <- 1
#' nu['E',c('kOTpE*OT*E')] <- -1
#' nu['O',c('kOT*OT','vdegrad*OTE','vdegrad*OT','kC*OC')] <- 1
#' nu['O',c('kOpT*O*Tt')] <- -1
#' nu['OCE',c('kcleav*OTE')] <- 1
#' nu['OCE',c('kOTE*OCE')] <- -1
#' nu['OC',c('kOTE*OCE')] <- 1
#' nu['OC',c('kC*OC')] <- -1
#' 
#' Trelstoc(0.1)$Tstat
#' @export
Trelstoc <- function(Ot,kOT=0.06){
  parms1['kOT'] <- kOT; parms1['kC'] <- kOT/0.6
  x0['O'] <- Ot
  Gillespie <- ssa(x0=x0,a=a,nu=nu,parms = parms1,tf=2E2,method = "ETL")
  data <- rowSums(Gillespie$data[,2:4])/x0['Tt']
  Tmean <- mean(data[200:nrow(Gillespie$data)])
  Tsd <- sd(data[200:nrow(Gillespie$data)])
  return(list(Trel=Gillespie$data,Tstat=c('O'=x0['O'],TrelM=Tmean,TrelSD=Tsd)))
}


