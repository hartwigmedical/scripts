
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

  distance = 0
  previous = allGenes[startIndex, ]
  end = allGenes[endIndex, ]
  for (i in ((startIndex+1):endIndex)) {
    current = allGenes[i,]
    if (current$end > previous$end || current$gene == end$gene) {
      distance = distance + 1
      previous = current
    }
  }
  return (distance)
}




