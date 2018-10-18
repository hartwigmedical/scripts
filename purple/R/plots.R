

plot_absolute_contribution <- function(contribution) {
  require(tidyr)
  require(dplyr)

  m_contribution = contribution %>% gather(variable, value, -1) %>% mutate(variable = factor(variable, levels = names(contribution)[-1]))
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  plot = ggplot(m_contribution, aes(x = factor(Sample),  y = Contribution, fill = factor(Signature), order = Sample)) +
    geom_bar(stat = "identity", colour = "black") +
    xlim(rev(levels(factor(m_contribution$Sample)))) +
    labs(x = "", y = "Absolute contribution") +
    scale_fill_discrete(name = "") + theme_bw() +
    theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  require(gridExtra)

  plots <- c(list(...), plotlist)
  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
