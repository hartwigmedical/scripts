options(Ncpus=8L, repos="https://cloud.r-project.org/")
# default packages: might fail with permissions error
install.packages(c("foreign","MASS", "survival"))
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
	pkgs=c(
	"copynumber",
	"StructuralVariantAnnotation",
	"VariantAnnotation",
	"rtracklayer",
	"BSgenome",
	"Rsamtools",
	"biomaRt",
	"org.Hs.eg.db",
	"TxDb.Hsapiens.UCSC.hg19.knownGene",
	"TxDb.Hsapiens.UCSC.hg38.knownGene",
	"BSgenome.Hsapiens.UCSC.hg19",
	"BSgenome.Hsapiens.UCSC.hg38"))
