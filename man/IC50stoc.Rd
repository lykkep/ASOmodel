\name{IC50stoc}
\alias{IC50stoc}
\title{Calculates the IC50 value for a stochastic simulation of the ASO model}
\usage{
  IC50stoc(Trel, Ot)
}
\arguments{
  \item{Trel}{The dissociation constant of the OT complex
  in nM}

  \item{Ot}{The total concentration of oligonucleotide}
}
\value{
  returns the IC50 value
}
\description{
  Use the stochastic dose-respons curve to calculate the
  IC50 value concentration before the addition of
  oligonucleotide.
}
\examples{
parms1 <- c(kOpT = 2E-5,kOTpE =50E-5 ,vprod = 150,  vdegrad = 0.04,
kcleav = 2, kOT =0.06, kOTE=2, kC = 0.1)
#Initital state vector
x0 <- c(Tt=parms1["vprod"]/parms1["vdegrad"],
             OT=0,OTE=0,E=1e3,O=1e5,OCE=0,OC=0)
names(x0) <- c('Tt','OT','OTE','E','O','OCE','OC')
#Propensity vector
a <-  c("vprod","kOpT*O*Tt","vdegrad*Tt","kOT*OT","kOTE*OTE","vdegrad*OT",
         "kOTpE*OT*E","vdegrad*OTE","kcleav*OTE","kC*OC","kOTE*OCE" )
#State-change matrix
nu <- matrix(0,7,length(a))
dimnames(nu) <- list(names(x0),a)
nu['Tt',c('vprod','kOT*OT')] <- 1
nu['Tt',c('kOpT*O*Tt','vdegrad*Tt')] <- -1
nu['OT',c('kOpT*O*Tt','kOTE*OTE')] <- 1
nu['OT',c('kOT*OT','kOTpE*OT*E','vdegrad*OT')] <- -1
nu['OTE',c('kOTpE*OT*E')] <- 1
nu['OTE',c('kOTE*OTE','vdegrad*OTE','kcleav*OTE')] <- -1
nu['E',c('kOTE*OTE','vdegrad*OTE','kOTE*OCE')] <- 1
nu['E',c('kOTpE*OT*E')] <- -1
nu['O',c('kOT*OT','vdegrad*OTE','vdegrad*OT','kC*OC')] <- 1
nu['O',c('kOpT*O*Tt')] <- -1
nu['OCE',c('kcleav*OTE')] <- 1
nu['OCE',c('kOTE*OCE')] <- -1
nu['OC',c('kOTE*OCE')] <- 1
nu['OC',c('kC*OC')] <- -1

Otseq <- 10^seq(2.5,6,by=0.2)
Trelseq <-sapply(Otseq,function(i) Trelstoc(i)$Tstat)})
IC50STOC <- IC50stoc(Trelseq, Otseq)
}

