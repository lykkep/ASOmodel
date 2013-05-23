\name{IC50}
\alias{IC50}
\title{Calculates the IC50 value at steady-state}
\usage{
  IC50(KdOT, param = parms)
}
\arguments{
  \item{KdOT}{The dissociation constant of the OT complex
  in nM}

  \item{parms}{list of parameters}
}
\value{
  returns the IC50 value
}
\description{
  Use the dose-respons curve to calculate the IC50 value
  concentration before the addition of oligonucleotide.
}
\examples{
IC50(1,c(Et = 1,KdOT = 0.3,kOpT = 0.2,KdOTE = 70,
kOTpE = 5,vprod = 0.2,vdegrad = 0.04,alpha=0.1,kcleav = 8))
}

