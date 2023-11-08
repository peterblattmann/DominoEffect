\name{map_to_func_elem}
\alias{map_to_func_elem}
\title{Functional annotation of significant hotspot residues.}
\usage{
  map_to_func_elem(hotspot_results, write_to_file = "NO", ens_release = "109")
}
\arguments{
  \item{hotspot_results}{Object containining information about the hotspot residues identified using the function identify_hotspots().}

  \item{write_to_file}{A character vector defining whether the resulting annotated hotspots should be saved in a file (YES or NO). }

  \item{ens_release}{A character vector defining whether the default gene annotations are used, i.e. Ensembl release 109, or if the gene_data correspond to a different Ensembl release. For the current Ensembl version this should be set to: ens_release="www.ensembl.org". For the archive versions to the corresponding archive website.}

}

\value{Updated results file containing an additional columns with the information on the annotated functional and structural region within which the mutation is mapped.}
\description{This function retrieves Uniprot annotations for th efunctional elements in the proteins with significant hotspots and overlaps and maps the hotspot residues to these.}
\examples{
data("TestData", package = "DominoEffect")
data("DominoData", package = "DominoEffect")
data("SnpData", package = "DominoEffect")

hotspot_mutations <- identify_hotspots(TestData, DominoData, SnpData)
hotspot_mutations <- map_to_func_elem(hotspot_mutations, 
write_to_file = "NO", ens_release = "109")

head(hotspot_mutations)
}
\author{
  Marija Buljan <buljan@imsb.biol.ethz.ch>
  Peter Blattmann <blattmann@imsb.biol.ethz.ch>
}
