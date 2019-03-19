options(Ncpus=8L, repos="https://cloud.r-project.org/")
update.packages(ask=FALSE)
# default packages: might fail with permissions error
install.packages(c(
	"foreign",
	"MASS",
	"survival"))
install.packages(c(
	"tidyverse",
	"devtools",
	"assertthat",
	"testthat",
	"NMF",
	"stringdist",
	"stringr",
	"argparser",
	"R.cache"))
if (!requireNamespace("BiocManager", quietly=TRUE)) {
	install.packages("BiocManager")
}
BiocManager::install(ask=FALSE,
	#version = "devel",
	pkgs=c(
	"copynumber",
	"VariantAnnotation",
	"rtracklayer",
	"BSgenome",
	"Rsamtools",
	"biomaRt",
	"org.Hs.eg.db",
	"TxDb.Hsapiens.UCSC.hg19.knownGene",
	"BSgenome.Hsapiens.UCSC.hg19"))
devtools::install_github("PapenfussLab/StructuralVariantAnnotation")
