last <- function(x) { tail(x, n = 1) }
parms <- c(Et = 1,KdOT = 0.3,kOpT = 0.2,KdOTE = 70,
           kOTpE = 5,vprod = 0.2,kdegrad = 0.04,alpha=0.1,kcleav = 8)

###############################################################################
#### The ODEs of the ASO model
###############################################################################
#TODO: document this function
#' @export
diffASO <- function(t,y, param){
  k1 =param['kOpT']; D1 = param['KdOT']; k2 = param['kOTpE']; Et <- param['Et']
  D2 = param['KdOTE']; vt = param['vprod']; k4 = param['kdegrad']
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
#### Plot dose-response curves for various parameters
###############################################################################
#TODO: document
#' @export
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
#' @export
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

