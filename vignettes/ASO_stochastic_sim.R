### INTEGRATE THIS CODE INTO THE SUPPLEMENTARY VIGNETTE. Or make it Rnw and source it into the Vignette

library(GillespieSSA)

###############################################################################
#### Stochastic
###############################################################################
a <-  c("vt","k1*O*Tt","kd*Tt","km1*OT","km2*OTE","kd*OT",
        "k2*OT*E","kd*OTE","kE*OTE","k3*OD","km2*ODE","k3*ODE" )

nu <- matrix(NA,8,length(a))
dimnames(nu) <- list(c("Tt","OT","OTE","E","O","ODE","OD","OE"),a)
#T
nu[1,] <- c(+1,-1,-1,+1,0,-1,rep(0,7))
#OT
nu[2,] <- c(0,+1,0,-1,+1,-1,-1,rep(0,6))
#OTE
nu[3,] <- c(0,0,0,0,-1,0,+1,-1,-1,rep(0,4))
#E
nu[4,] <- c(0,0,0,0,1,0,-1,1,0,1,0,1,0)
#O
nu[5,] <- c(0,-1,0,1,0,1,0,1,0,1,1,0,0)
#ODE
nu[6,] <- c(rep(0,8),1,0,0,-1,-1)
#OD
nu[7,] <- c(rep(0,10),-1,1,0)
#OE
nu[8,] <- c(rep(0,9),-1,0,0,1)

Ag = 6.0221415E23
Kj <- 2
#5E6/(2.25E-12*Ag)

parms1 <- c(k1 = 0.06/135000,k2 =5/135000 ,vt = 0.1*13500,	kd = 0.2,		  
              kE = 2 km1 =0.03, km2=0.01, k3 = 0.09)
x0 <- c(parms1["vt"]/parms1["kd"],OT=0,OTE=0,E=10*1350,O=50*1350,ODE=0,OD=0,OE=0)
names(x0)[1] <- 'Tt'

Gillespie <- ssa( x0=x0,# initial state vector
      a=a, # propensity vector
      nu=nu, # state-change matrix
      parms = parms1, # model parameters
      tf=5E2, # final time
      method = "ETL" # SSA method
)
range(rowSums(Gillespie$data[,c(3,4,6:9)])-x0['O'])
range(rowSums(Gillespie$data[,c(4,5,7,9)])-x0['E'])
plot(Gillespie$data[,1],rowSums(Gillespie$data[,2:4])/x0['Tt'],
		xlab='time',ylab='Relative target level',type='l')

Trel <- function(Ot,km1=0.31){  # TODO: should this be in the package functions instead?  Perhaps made as similar possible to the ODE version?
	# TODO: Can we avoid hard-coded parameters in functions?
  parms1 <- c(k1 = 0.06/135000,k2 =5/135000 ,vt = 0.5*13500,  kd = 0.2,		  
	             kE = 0.1,kB = 0.1, km1 =km1, km2=0.01, k3 = km1/0.5)
	x0 <- c(parms1["vt"]/parms1["kd"],OT=0,OTE=0,E=10*1350,O=Ot*1350,ODE=0,OD=0,OE=0)
	names(x0) <- c('Tt','OT','OTE','E','O','ODE','OD','OE')
	Gillespie <- ssa(x0=x0,a=a,nu=nu,parms = parms1,tf=2E2,method = "ETL")
	data <- rowSums(Gillespie$data[,2:4])/x0['Tt']
	Tmean <- mean(data[250:nrow(Gillespie$data)])
	Tsd <- sd(data[250:nrow(Gillespie$data)])
	return(list(Trel=Gillespie$data,Tstat=c('O'=x0['O'],TrelM=Tmean,TrelSD=Tsd)))
}

tmp <- Trel(100)
plot(tmp$Trel[,1],rowSums(tmp$Trel[,2:4])/(parms1["vt"]/parms1["kd"]),ylab='Trel',xlab='time')

#time T    OT   OTE     E      O    ODE     OD    OE
# 105 0.04 1.97 0.00239 0.0973 97.4 0.00034 0.664 1.23e-05

###############################################################################
#### IC50
###############################################################################

#rm(DRcurve)
KM <- 10^seq(-7,4,length=10)
DRcurve <- lapply(KM,function(kmi){ 
				sapply(10^seq(1,3.8,by=0.1),function(i) Trel(i,km1=kmi)$Tstat)})

IC50_KM <- sapply(1:length(DRcurve),function(x){
  IC50stoc(DRcurve[[x]][2,!is.na(DRcurve[[x]][2,])],DRcurve[[x]][1,!is.na(DRcurve[[x]][2,])])})

#quartz('IC50_vs_Dot')
plot(KM,IC50_KM,log='x',ylab=expression(IC[50]),
		xlab=expression(k[OT->O+T]))

#pdf('Stochastic.pdf',width=8,height=5)
par(lwd=2)
N <- length(DRcurve)
plot(DRcurve[[1]][1,],DRcurve[[1]][2,],ylab='Relative target level',
		xlab='Total oligo (molecules)',log='x',type='b',col=rainbow(N)[1],xlim=c(1E4,1.4E7))
legend('bottomleft',paste('k(OT->O+T)=',signif(KM,2)),
		col=rainbow(N),lwd=2,bty='n')
for(i in 2:N){
	points(DRcurve[[i]][1,],DRcurve[[i]][2,],col=rainbow(N)[i],type='b')
}
#dev.off()

###############################################################################
#### Stochastic for 2 target sites  # TODO: not part of first paper. Move out.
###############################################################################
a <-  c("(1-OT1)*(1-OT2)*(1-T1)*(1-T2)*(1-OT1E)*(1-OT2E)*vt",
		"kd*(T1+T2+OT1+OT2+OT1E+OT2E)/2",
		"k1*O*T1*(1-OT1)","km1*OT1","k2*OT1*E*(1-OT1E)","km2*OT1E",
		"k1*O*T2*(1-OT2)","km1*OT2","k2*OT2*E*(1-OT2E)","km2*OT2E",
        "kE*OT1E*OT2E","kE*OT1E*OT2","kE*OT1E*T2",
        "kE*OT2E*OT1","kE*OT2E*T1")

aa <- c(	gsub('T2','T21',gsub('T1','T11',a)),
		gsub('T2','T22',gsub('T1','T12',a)),"kB*OE","k3*OD","km2*ODE","k3*ODE" )

nu <- matrix(NA,6,length(a))
dimnames(nu) <- list(c("T1","T2","OT1","OT2","OT1E","OT2E"),a)

#T1
nu[1,] <- c(+1,0,-1,1,rep(0,10),-1)
#T2
nu[2,] <- c(+1,0,rep(0,4),-1,1,rep(0,4),-1,rep(0,2))
#OT1
nu[3,] <- c(rep(0,2),+1,-1,-1,1,rep(0,7),-1,0)
#OT2
nu[4,] <- c(rep(0,6),+1,-1,-1,+1,0,-1,rep(0,3))
#OT1E
nu[5,] <- c(rep(0,4),1,-1,rep(0,4),-1,-1,-1,rep(0,2))
#OT2E
nu[6,] <- c(rep(0,8),1,-1,-1,0,0,-1,-1)

nunu <- matrix(0,2*6+5,length(aa))
colnames(nunu) <- aa

nunu[1:6,1:length(a)] <- nu
nunu[7:12,(length(a)+1):(2*length(a))] <- nu

#ODE
nunu[13,] <- c(rep(c(rep(0,10),1,1,1,1,1),2),0,0,-1,-1)
#OD
nunu[14,32:34] <- c(-1,1,0)
#OE
nunu[15,31:34] <- c(-1,0,0,1)
#E
nunu[16,] <- c(rep(c(rep(0,4),-1,1,0,0,-1,1,1,0,0,0,0),2),1,0,1,0)
#O
nunu[17,] <- c(rep(c(rep(0,2),-1,1,0,0,-1,1,0,0, 1,1,0,1,0  ),2),1,1,0,0)

parms1 <- c(k1 = 0.6,k2 =0.1 ,vt = 0.1,	kd = 0.02,		  
              kE = 0.1,kB = 0.1, km1 =0.03, km2=0.01, k3 = 0.09)
i0 <- c(rep(0,18),5,100)
navn <- c('T1','T2','OT1','OT2','OT1E','OT2E','ODE','OD','OE')
names(i0) <- c(	gsub('T2','T21',gsub('T1','T11',navn)),
			gsub('T2','T22',gsub('T1','T12',navn)),'E','O')

Gillespie <- ssa( # initial state vector
      x0=i0,
      a=aa, # propensity vector
      nu=nunu, # state-change matrix
      # model parameters
      parms = parms1, 
     tf=50, # final time
     method = "D" # SSA method
    )

Gillespie$data[nrow(Gillespie$data)-10:0,]

range(rowSums(Gillespie$data[,c(3,4,6:9)])-x0['O'])
range(rowSums(Gillespie$data[,c(4,5,7,9)])-x0['E'])
plot(Gillespie$data[,1],rowSums(Gillespie$data[,2:4])/x0['Tt'],
		xlab='time',ylab='Relative target level',type='l')


###############################################################################
#### Deterministic
###############################################################################
# diffOTEex <- function(t,y, parms1){
  # k1 =parms1["k1"]; D1 = parms1["D1"]; k2 = parms1["k2"]
  # D2 = parms1["D2"]; vt = parms1["vt"]; k4 = parms1["kd"]; 
  # alpha=parms1["alpha"]; kE=parms1["kE"]; kB=parms1["kB"]
  # km1 = k1*D1; km2=k2*D2; k3 = km1/alpha
  
  # Tt = y[1]; OT = y[2]; OTE= y[3]; E=y[4]; O= y[5]
  # ODE= y[6]; OD = y[7]; OE= y[8]
  
  # F <- {}   
  # 1= T
  # F[1]= vt-k1*O*Tt-k4*Tt+km1*OT
  # 2= OT
  # F[2]= k1*Tt*O-km1*OT-k2*OT*E+km2*OTE-k4*OT
  # 3= OTE
  # F[3]= k2*OT*E-km2*OTE-k4*OTE-kE*OTE
  # 4= E
  # F[4]= -k2*(OT)*E+km2*(OTE+ODE)+k4*(OTE)+kB*OE
  # 5= O
  # F[5]=km1*OT-k1*Tt*O+k4*(OT+OTE)+k3*OD+kB*OE
  # ODE
  # F[6] = kE*OTE-km2*ODE-k3*ODE
  # OD
  # F[7]= +km2*ODE-k3*OD
  # OE
  # F[8] = k3*ODE-kB*OE
  
  # return(list(F))
# }

xEX <- c(T=parms["vt"]/parms["kd"],OT=0,OTE=0,E=parms["etot"],O=100,ODE=0,OD=0,OE=0)
names(xEX) <- c("T","OT","OTE","E","O","ODE","OD","OE")
seqList <- list(seq(0,3.5,by=0.05)/60,seq(0.1+5E-2,4.3,by=5E-2),seq(5,105,by=1))
timeseqLenght <- sapply(seqList,length)
solEXt <- vode(xEX,unlist(seqList),diffOTEex,parms)
