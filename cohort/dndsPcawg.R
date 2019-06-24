dndsFiles = list.files(path = "~/hmf/dnds", pattern = "dNdScv_output_")

PcawgRefCDSCv = data.frame()
for (filename in dndsFiles) {
  cat("Processing ", filename, "\n")
  file = paste("~/hmf/dnds/", filename, sep = "")
  df = read.delim(file, stringsAsFactors = FALSE)
  PcawgRefCDSCv = rbind(PcawgRefCDSCv, df)
}

save(PcawgRefCDSCv, file = "~/hmf/dnds/PcawgRefCDSCv.RData")
