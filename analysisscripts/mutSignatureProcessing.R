library(devtools)
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library(grid)
library(gridExtra)
library("NMF")
library(ggplot2)

myCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
             "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
             "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
             "#dea185","#a0729d","#8a392f")


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

}


plot_fitted_signatures<-function(fit_contribution,cancer_signatures,cancerType,fileId,chart_mode="absolute",writePDF=FALSE){
  if (writePDF){
    pdf(file=paste(cancerType,fileId,chart_mode,".pdf",sep=""),width=10)
  }
  fit_contribution[prop.table(fit_contribution, margin=2)<0 | fit_contribution<0]<-0
  orderVector<-colSums(fit_contribution)
  fit_contribution=fit_contribution[, order(orderVector,decreasing=F),drop=FALSE]
  colnames(fit_contribution)<-paste(colnames(fit_contribution),round(colSums(fit_contribution)))
  numSamples = ncol(fit_contribution)
  maxSamples = 90
  print(numSamples)
  if (numSamples > 1) {
    for (n in 1:ceiling(numSamples/maxSamples)){
      minCol = ((n-1)*maxSamples+1)
      maxCol = min((n*maxSamples),numSamples)
      p1<-plot_contribution(fit_contribution[,minCol:maxCol,drop=F], cancer_signatures,coord_flip = F, mode = chart_mode)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7),legend.text=element_text(size=5),axis.title.y = element_text(size=6))+
           scale_fill_manual( values= myCOLORS)+labs(fill="")+ggtitle(paste(cancerType," sorted by total mutational load asc"))
      print(p1)
    }
  }

  if (writePDF){
    dev.off()
  }
}

###################################
dataFile = "~/hmf/mutSignature3.RData"
load(dataFile)

# distinct cancer types, excluding NA
cancer_types = na.omit(unique(cohort$cancerType))
cancer_types = sort(cancer_types)

# filter by cancerType
sigMatrix=sapply(signatures, function(x) x[,"subclonal"])
#sigMatrix[prop.table(sigMatrix, margin=2)<0 | sigMatrix<0]<-0
#orderVector<-colSums(sigMatrix)
#sigMatrix=sigMatrix[, order(orderVector,decreasing=F),drop=FALSE]

for (cancer_type in cancer_types) {
  print(cancer_type)
  samples=cohort[cancerType == cancer_type]$sampleId
  plot_fitted_signatures(sigMatrix[,c(samples),drop=FALSE],cancer_signatures,cancer_type,'subclonal','relative',TRUE)
}

write.csv(signatures,file = '~/hmf/Signatures3.csv')
