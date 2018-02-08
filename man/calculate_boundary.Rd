\name{calculate_boundary}
\alias{calculate_boundary}
\title{calculate_boundary}
\usage{
calculate_boundary(mut_pos_numeric, length_aa, flanking_region)
}
\arguments{
  \item{mut_pos_numeric}{Amino acid position of mutation}

  \item{length_aa}{Length of transcript in amino acids}

  \item{flanking_region}{Vector containing two flanking regions}

}
\value{returns a list with the boundaries for the two regions}
\description{The function calculates boundaries of sequence windows around the mutation. It is possible to define up to two window lengths. If the mutation is close to the start or end of the protein, the region is extended into the other direction to keep the indicated size}
\examples{
calculate_boundary(250, 500, c(200, 300))
calculate_boundary(250, 500, 300)
}
\author{
  Marija Buljan <buljan@imsb.biol.ethz.ch>
  Peter Blattmann <blattmann@imsb.biol.ethz.ch>
}
