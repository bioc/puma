\name{puma-package}
\alias{puma-package}
\alias{puma}
\docType{package}
\title{
	puma - Propagating Uncertainty in Microarray Analysis
}
\description{
Most analyses of Affymetrix GeneChip data are based on point estimates of expression levels and ignore the uncertainty of such estimates. By propagating uncertainty to downstream analyses we can improve results from microarray analyses. For the first time, the puma package makes a suite of uncertainty propagation methods available to a general audience. puma also offers improvements in terms of scope and speed of execution over previously available uncertainty propagation methods. Included are summarisation, differential expression detection, clustering and PCA methods, together with useful plotting functions.
}
\details{
	\tabular{ll}{
		Package: \tab puma\cr
		Type: \tab Package\cr
		Version: \tab 1.2.1\cr
		Date: \tab 2007-06-26\cr
		License: \tab LGPL excluding donlp2\cr
	}
	For details of using the package please refer to the Vignette
}
\author{
	Richard Pearson, Xuejun Liu, Guido Sanguinetti, Marta Milo, Neil D. Lawrence, Magnus Rattray

	Maintainer: Richard Pearson <richard.pearson@postgrad.manchester.ac.uk>
}
\references{
	Milo, M., Niranjan, M., Holley, M. C., Rattray, M. and Lawrence, N. D. (2004) A probabilistic approach for summarising oligonucleotide gene expression data, technical report available upon request.
	
	Liu, X., Milo, M., Lawrence, N. D. and Rattray, M. (2005) A tractable probabilistic model for Affymetrix probe-level analysis across multiple chips, Bioinformatics, 21(18):3637-3644.
	
	Sanguinetti, G., Milo, M., Rattray, M. and Lawrence, N. D. (2005) Accounting for probe-level noise in principal component analysis of microarray data, Bioinformatics, 21(19):3748-3754.

	Rattray, M., Liu, X., Sanguinetti, G., Milo, M. and Lawrence, N. D. (2006) Propagating uncertainty in Microarray data analysis, Briefings in Bioinformatics, 7(1):37-47. 
	
	Liu, X., Milo, M., Lawrence, N. D. and Rattray, M. (2006) Probe-level measurement error improves accuracy in detecting differential gene expression, Bioinformatics, 22(17):2107-2113.
	
	Liu, X. Lin, K., Andersen, B. Rattray, M. (2007) Including probe-level uncertainty in model-based gene expression clustering, BMC Bioinformatics, 8(98).
	
	Pearson, R. D., Liu, X., Sanguinetti, G., Milo, M., Lawrence, N. D., Rattray, M. (2007) In preparation.
}
\examples{
	data(Dilution)
	varLabels(Dilution) <- "level"
	pData(Dilution) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(Dilution)))
	eset_mmgmos <- mmgmos(Dilution)
	pumapca_mmgmos <- pumaPCA(eset_mmgmos)
	plot(pumapca_mmgmos)
	eset_mmgmos_normd <- pumaNormalize(eset_mmgmos)
	eset_comb <- pumaComb(eset_mmgmos_normd)
	esetDE <- pumaDE(eset_comb)
}
\keyword{ package }
