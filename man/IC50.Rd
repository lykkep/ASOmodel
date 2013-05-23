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
  Use the algebraic solution to the system to the total
  target concentration at steady state relative to the
  target concentration before the addition of
  oligonucleotide.
}
\examples{
IC50(1,parms)
}

