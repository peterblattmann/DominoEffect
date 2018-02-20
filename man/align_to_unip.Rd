\name{align_to_unip}
\alias{align_to_unip}
\title{Align protein segnent around the hotspot to the UniProt/Swiss-Prot KB sequence.}
\usage{
  align_to_unip(ens.seq, uni.seq, ensembl_mut_position)
}
\arguments{
  \item{ens.seq}{AAString object with the Ensembl protein sequence corresponding to the representative transcript.}

  \item{uni.seq}{AAString with the UniProt sequence for the identifier matching the Ensembl gene name.}

  \item{ensembl_mut_position}{Numeric vector defining the hotspot position in the Ensembl sequence, i.e. in the ens.seq}

}
\value{Returns a list where the first element is a character vector defining whether there was a significant alignment and the second element provides the hotspot position in the UniProt sequence.}

\description{This function alignes the Ensembl protein region with a hotspot to the UniProt sequence. The Ensembl region encompasses 15 amino acids where the hotspot is in the middle. If the hotspot was at the protein start or end the region is still 15 amino acids long, but the hotspot position is shifted.}

\examples{
library(Biostrings)

ens.seq <- AAString("MDLSALREEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLK")
uni.seq <- AAString("MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLA")
ensembl_mut_position <- 25

align_to_unip(ens.seq, uni.seq, ensembl_mut_position)
}
\author{
  Marija Buljan <buljan@imsb.biol.ethz.ch>
  Peter Blattmann <blattmann@imsb.biol.ethz.ch>
}
