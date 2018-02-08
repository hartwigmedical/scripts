library(purple)
context("Common")

test_that("PatientId in manual overrides", {
  lookup = data.frame(sample=character(), patient=character())

  expect_equal(sample_to_patient_id("CPCT02020192", lookup), "CPCT02020438")
  expect_equal(sample_to_patient_id("CPCT02020192T", lookup), "CPCT02020438")
  expect_equal(sample_to_patient_id("CPCT02020192TII", lookup), "CPCT02020438")
  expect_equal(sample_to_patient_id("CPCT02020192TIII", lookup), "CPCT02020438")
})

test_that("PatientId not in lookup ", {
  lookup = data.frame(sample=character(), patient=character())

  expect_equal(sample_to_patient_id("CPCT02020193", lookup), "CPCT02020193")
  expect_equal(sample_to_patient_id("CPCT02020193T", lookup), "CPCT02020193")
  expect_equal(sample_to_patient_id("CPCT02020193TII", lookup), "CPCT02020193")
  expect_equal(sample_to_patient_id("CPCT02020193TIII", lookup), "CPCT02020193")
})

test_that("PatientId in lookup ", {
  sample  = c("CPCT02020193")
  patient = c("CPCT02020194")
  lookup = data.frame(sample, patient, stringsAsFactors = FALSE)

  expect_equal(sample_to_patient_id("CPCT02020193", lookup), "CPCT02020194")
  expect_equal(sample_to_patient_id("CPCT02020193T", lookup), "CPCT02020194")
  expect_equal(sample_to_patient_id("CPCT02020193TII", lookup), "CPCT02020194")
  expect_equal(sample_to_patient_id("CPCT02020193TIII", lookup), "CPCT02020194")
})

