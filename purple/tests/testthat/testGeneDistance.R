library(purple)
context("GeneDistance")

test_that("Non overlapping genes are easy", {

  genes = data.frame(chromosome=character(), start=integer(), end=integer(), gene=character())
  genes = data.frame(chromosome="X", start=197859, end=220023, gene="PLCXD1")
  genes <- rbind(genes, data.frame(chromosome="X", start=220025, end=230886, gene="GTPBP6"))
  genes <- rbind(genes, data.frame(chromosome="X", start=294698, end=347445, gene="PPP2R3B"))

  expect_equal(gene_distance("PLCXD1", "PLCXD1", genes), 0)
  expect_equal(gene_distance("PLCXD1", "GTPBP6", genes), 1)
  expect_equal(gene_distance("PLCXD1", "PPP2R3B", genes), 2)

  expect_equal(gene_distance("GTPBP6", "PLCXD1", genes), -1)
  expect_equal(gene_distance("GTPBP6", "GTPBP6", genes), 0)
  expect_equal(gene_distance("GTPBP6", "PPP2R3B", genes), 1)

  expect_equal(gene_distance("PPP2R3B", "PLCXD1", genes), -2)
  expect_equal(gene_distance("PPP2R3B", "GTPBP6", genes), -1)
  expect_equal(gene_distance("PPP2R3B", "PPP2R3B", genes), 0)
})

create_overlapping_genes<-function() {
  genes = data.frame(chromosome=character(), start=integer(), end=integer(), gene=character())
  genes = data.frame(chromosome="X", start=31089360, end=31090170, gene="FTHL17")
  genes <- rbind(genes, data.frame(chromosome="X", start=31137345, end=33229636, gene="DMD"))
  genes <- rbind(genes, data.frame(chromosome="X", start=32601773, end=32601869, gene="MIR3915"))
  genes <- rbind(genes, data.frame(chromosome="X", start=32659591, end=32659676, gene="MIR548F5"))
  genes <- rbind(genes, data.frame(chromosome="X", start=34147869, end=34150447, gene="FAM47A"))
  return (genes[order(genes$start), ])
}

test_that("Overlapping genes are harder", {
  allGenes = create_overlapping_genes()
  expect_equal(gene_distance("FTHL17", "FTHL17", allGenes), 0)
  expect_equal(gene_distance("FTHL17", "DMD", allGenes), 1)
  expect_equal(gene_distance("FTHL17", "MIR3915", allGenes), 2)
  expect_equal(gene_distance("FTHL17", "MIR548F5", allGenes), 2)
  expect_equal(gene_distance("FTHL17", "FAM47A", allGenes), 2)

  expect_equal(gene_distance("DMD", "FTHL17", allGenes), -1)
  expect_equal(gene_distance("DMD", "DMD", allGenes), 0)
  expect_equal(gene_distance("DMD", "MIR3915", allGenes), 1)
  expect_equal(gene_distance("DMD", "MIR548F5", allGenes), 1)
  expect_equal(gene_distance("DMD", "FAM47A", allGenes), 1)

  expect_equal(gene_distance("MIR3915", "FTHL17", allGenes), -2)
  expect_equal(gene_distance("MIR3915", "DMD", allGenes), -1)
  expect_equal(gene_distance("MIR3915", "MIR3915", allGenes), 0)
  expect_equal(gene_distance("MIR3915", "MIR548F5", allGenes), 1)
  expect_equal(gene_distance("MIR3915", "FAM47A", allGenes), 2)

  expect_equal(gene_distance("MIR548F5", "FTHL17", allGenes), -2)
  expect_equal(gene_distance("MIR548F5", "DMD", allGenes), -1)
  expect_equal(gene_distance("MIR548F5", "MIR3915", allGenes), -1)
  expect_equal(gene_distance("MIR548F5", "MIR548F5", allGenes), 0)
  expect_equal(gene_distance("MIR548F5", "FAM47A", allGenes), 1)

  expect_equal(gene_distance("FAM47A", "FTHL17", allGenes), -2)
  expect_equal(gene_distance("FAM47A", "DMD", allGenes), -1)
  expect_equal(gene_distance("FAM47A", "MIR3915", allGenes), -2)
  expect_equal(gene_distance("FAM47A", "MIR548F5", allGenes), -1)
  expect_equal(gene_distance("FAM47A", "FAM47A", allGenes), 0)
})
