\name{import_txdb}
\alias{import_txdb}
\title{Imports txdb data and converts it into format required for DominoEffect package}
\usage{
  import_txdb(txdb_object)
}
\arguments{
  \item{txdb_object}{TxDB Object, e.g. from makeTxDbFromEnsembl}

}

\value{Data frame that can be used in DominoEffect package.}

\description{This function imports txdb data and converts into a data frame used in the DominoEffect package only extracting the required information from the txdb object.}
\examples{

#EnsTxDB <- makeTxDbFromEnsembl(organism="Homo sapiens", release=73, 
#                               server="ensembldb.ensembl.org")
#DominoData <- import_txdb(EnsTxDB)
#head(DominoData)

}
\author{
  Marija Buljan <buljan@imsb.biol.ethz.ch>
  Peter Blattmann <blattmann@imsb.biol.ethz.ch>
}
