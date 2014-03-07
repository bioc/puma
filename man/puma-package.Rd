\name{puma-package}
\alias{puma-package}
\alias{puma}
\docType{package}
\title{
	puma - Propagating Uncertainty in Microarray Analysis
}
\description{
Most analyses of Affymetrix GeneChip data (including tranditional 3' arrays and exon arrays) are based on point estimates of expression levels and ignore the uncertainty of such estimates. By propagating uncertainty to downstream analyses we can improve results from microarray analyses. For the first time, the puma package makes a suite of uncertainty propagation methods available to a general audience. In additon to calculte gene expression from Affymetrix 3' arrays, puma also provides methods to process exon arrays and produces gene and isoform expression for alternative splicing study. puma also offers improvements in terms of scope and speed of execution over previously available uncertainty propagation methods. Included are summarisation, differential expression detection, clustering and PCA methods, together with useful plotting functions.
}
\details{
	\tabular{ll}{
		Package: \tab puma\cr
		Type: \tab Package\cr
		Version: \tab 3.4.3\cr
		Date: \tab 2013-11-04\cr
		License: \tab LGPL excluding donlp2\cr
	}
	For details of using the package please refer to the Vignette
}
\author{
	Richard Pearson, Xuejun Liu, Guido Sanguinetti, Marta Milo, Neil D. Lawrence, Magnus Rattray, Li Zhang

	Maintainer: Richard Pearson <richard.pearson@postgrad.manchester.ac.uk>, Li Zhang <leo.zhang@nuaa.edu.cn>
}
\references{
	Milo, M., Niranjan, M., Holley, M. C., Rattray, M. and Lawrence, N. D. (2004) A probabilistic approach for summarising oligonucleotide gene expression data, technical report available upon request.
	
	Liu, X., Milo, M., Lawrence, N. D. and Rattray, M. (2005) A tractable probabilistic model for Affymetrix probe-level analysis across multiple chips, Bioinformatics, 21(18):3637-3644.
	
	Sanguinetti, G., Milo, M., Rattray, M. and Lawrence, N. D. (2005) Accounting for probe-level noise in principal component analysis of microarray data, Bioinformatics, 21(19):3748-3754.

	Rattray, M., Liu, X., Sanguinetti, G., Milo, M. and Lawrence, N. D. (2006) Propagating uncertainty in Microarray data analysis, Briefings in Bioinformatics, 7(1):37-47. 
	
	Liu, X., Milo, M., Lawrence, N. D. and Rattray, M. (2006) Probe-level measurement error improves accuracy in detecting differential gene expression, Bioinformatics, 22(17):2107-2113.
	
	Liu, X. Lin, K., Andersen, B. Rattray, M. (2007) Including probe-level uncertainty in model-based gene expression clustering, BMC Bioinformatics, 8(98).

	Pearson, R. D., Liu, X., Sanguinetti, G., Milo, M., Lawrence, N. D., Rattray, M. (2008)  puma: a Bioconductor package for Propagating Uncertainty in Microarray Analysis, BMC Bioinformatics, 2009, 10:211.

          Zhang,L. and Liu,X. (2009) An improved probabilistic model for finding differential gene expression, the 2nd BMEI 17-19 oct. 2009. Tianjin. China.	

          Liu,X. and Rattray,M. (2009) Including probe-level measurement error in robust mixture clustering of replicated microarray gene expression, Statistical Application in Genetics and Molecular Biology, 9(1), Article 42.

          puma 3.0: improved uncertainty propagation methods for gene and transcript expression analysis, Liu et al. BMC Bioinformatics, 2013, 14:39.


}
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)
	pumapca_mmgmos <- pumaPCA(eset_mmgmos)
	plot(pumapca_mmgmos)
	eset_mmgmos_100 <- eset_mmgmos[1:100,]
	eset_comb <- pumaComb(eset_mmgmos_100)
        eset_combImproved <- pumaCombImproved(eset_mmgmos_100)
	esetDE <- pumaDE(eset_comb)
        esetDEImproved <- pumaDE(eset_combImproved)
}
\keyword{ package }
