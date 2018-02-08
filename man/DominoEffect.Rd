\name{DominoEffect}
\alias{DominoEffect}
\title{Identification of significant mutation hotspot residues.}
\usage{
  DominoEffect(mutation_dataset = TestData, gene_data = DominoData, snp_data = SnpData,
min_n_muts = 5, MAF_thresh = 0.01,
flanking_region = c(200, 300),
poisson.thr = 0.01, percentage.thr = 0.15,
ratio.thr = 45, approach = "percentage", write_to_file = "WRITE")
}
\arguments{
  \item{mutation_dataset}{Object containing a table with the mutation data (e.g. TCGA point mutations mapped to protein level).}

  \item{gene_data}{DominoData object containing information about Ensembl gene annotations: gene identifiers and regresentative transcript cDNA length.}

  \item{snp_data}{Object containing a table with information on population SNPs.}

  \item{min_n_muts}{Numeric vector defining a minimum number of mutations that need to occur at the same residue. Default: 5}

  \item{MAF_thresh}{Numeric vector that defines Minor allele frequency threshold for considering reported mutations as population SNPs.}

 \item{flanking_region}{Numeric vector that defines size of a window around the mutation that will be considered. Up to two window sizes are allowed.}

 \item{poisson.thr}{Numeric vector that defines a treshold for the adjusted p-value.
  Residues with an associated p-value that is lower than the defined
  value are reported. Default: 0.01}

  \item{percentage.thr}{Number defining the fraction of mutations within the window that need to fall on a single residue in order for it to be classified as a hotspot. Default: 0.15}

  \item{ratio.thr}{Number defining a requirement that a number of mutations on a single residue should exceed what would be expected by chance given a background mutation rate in the window (i.e. the surrounding region). Default: 45}

  \item{approach}{Option to define selection criteria to use precentage.thr or ratio.thr as criterion for finding single residue mutation clusters. Options: "both", "percentage" or "ratio". Default = "percentage"}

  \item{write_to_file}{Option if the identified and annotated hotspots should be written to a file. Default: WRITE}

}
\value{Results table}
\description{
The function identifies individual amino acid residues, which accumulate a high fraction of the overall mutation load within a protein. After detecting mutation hotspots, the function obtains UniProt/Swiss-Prot KB functional and structural annotations for the affected protein regions and checks for the sequence agreement.
}
\examples{
data("SnpData", package = "DominoEffect")
data("TestData", package = "DominoEffect")
data("DominoData", package = "DominoEffect")

hotspot_mutations <- DominoEffect(mutation_dataset = TestData, 
gene_data = DominoData, snp_data = SnpData)
hotspot_mutations <- DominoEffect(mutation_dataset = TestData)
}
\author{
  Marija Buljan <buljan@imsb.biol.ethz.ch>
  Peter Blattmann <blattmann@imsb.biol.ethz.ch>
}
