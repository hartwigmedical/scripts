library(purple)
library(dplyr)
library(GenomicRanges)
context("NearHotspot")

test_that("Del before hotspot", {

  chromosome = c(1, 1)
  position = c(100, 107)
  ref = c("GAT", "A")

  mutations = data.frame(chromosome = chromosome, position = position, hotspot = c(F, T), ref = ref, stringsAsFactors = F)
  expect_equal(nearHotspot(mutations, distance = 3), c(F, F))
  expect_equal(nearHotspot(mutations, distance = 4), c(F, F))
  expect_equal(nearHotspot(mutations, distance = 5), c(T, F))
})

test_that("Del after hotspot", {

  chromosome = c(1, 1)
  position = c(112, 107)
  ref = c("GAT", "A")

  mutations = data.frame(chromosome = chromosome, position = position, hotspot = c(F, T), ref = ref, stringsAsFactors = F)
  expect_equal(nearHotspot(mutations, distance = 3), c(F, F))
  expect_equal(nearHotspot(mutations, distance = 4), c(F, F))
  expect_equal(nearHotspot(mutations, distance = 5), c(T, F))
})

test_that("Really long del", {

  chromosome = c(1, 1)
  position = c(100, 107)
  ref = c("GATACAGATACAGATACAGATACA", "A")

  mutations = data.frame(chromosome = chromosome, position = position, hotspot = c(F, T), ref = ref, stringsAsFactors = F)
  expect_equal(nearHotspot(mutations, distance = 3), c(T, F))
  expect_equal(nearHotspot(mutations, distance = 4), c(T, F))
  expect_equal(nearHotspot(mutations, distance = 5), c(T, F))
})
