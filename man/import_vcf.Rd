\name{import_vcf}
\alias{import_vcf}
\title{Import vcf data and convert into format for DominoEffect package}
\usage{
  import_vcf(vcf, MAF_thresh = 0.01)
}
\arguments{
  \item{vcf}{Object of class CollapsedVCF}

  \item{MAF_thresh}{Threshold for frequency of observation. Default = 0.01. }

}

\value{Data frame that is used in DominoEffect package.}
\description{This function imports vcf data and converts into a data frame used in the DominoEffect package only extracting the required information from the vcf object.}
\examples{
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")

SnpData <- import_vcf(vcf)
head(SnpData)
}
\author{
  Marija Buljan <buljan@imsb.biol.ethz.ch>
  Peter Blattmann <blattmann@imsb.biol.ethz.ch>
}
