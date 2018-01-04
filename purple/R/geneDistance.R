
gene_distance<-function(startGeneName, endGeneName, allGenes) {
  if (startGeneName == endGeneName) {
    return (0)
  }

  startIndex = match(startGeneName, allGenes$gene)
  endIndex = match(endGeneName, allGenes$gene)

  if (startIndex > endIndex) {
    return (-1 * gene_distance(endGeneName, startGeneName, allGenes))
  }

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




