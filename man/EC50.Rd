\name{EC50}
\alias{EC50}
\title{Calculates the EC50 value at steady-state}
\usage{
  EC50(KdOT, param = parms)
}
\arguments{
  \item{KdOT}{The dissociation constant of the OT complex
  in nM}

  \item{parms}{list of parameters}
}
\value{
  returns the EC50 value
}
\description{
  Use the dose-respons curve to calculate the EC50 value
  concentration before the addition of oligonucleotide.
}
\examples{
EC50(1,c(Et = 1,KdOT = 0.3,kOpT = 0.2,KdOTE = 70,
kOTpE = 5,vprod = 0.2,vdegrad = 0.04,alpha=0.1,kcleav = 8))
}

