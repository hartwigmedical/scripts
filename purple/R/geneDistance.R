
gene_distance<-function(startGeneName, endGeneName, allGenes) {

  startIndex = match(startGeneName, allGenes$gene)
  endIndex = match(endGeneName, allGenes$gene)
  return (gene_distance_index(startIndex, endIndex, allGenes))
}

gene_distance_index<-function(startIndex, endIndex, allGenes) {
  if (startIndex == endIndex) {
    return (0)
  }

  if (startIndex > endIndex) {
    return (-1 * gene_distance_index(endIndex, startIndex, allGenes))
  }

  endGeneName = allGenes[endIndex, ]$gene

  distance = 0
  previous = allGenes[startIndex, ]
  for (i in ((startIndex+1):endIndex)) {
    current = allGenes[i,]
    if (current$end > previous$end || current$gene == endGeneName) {
      distance = distance + 1
      previous = current
    }
  }
  return (distance)
}




