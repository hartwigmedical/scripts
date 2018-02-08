library(purple)
context("Common")

test_that("Multiple Biopsy Special Case", {
  lookup = data.frame(sample=character(), patient=character())

  expect_equal(sample_to_patient_id("CPCT02020192T", lookup), "CPCT02020438")
  expect_equal(sample_to_patient_id("CPCT02030224T", lookup), "CPCT02030292")
  expect_equal(sample_to_patient_id("DRUP01010007T", lookup), "DRUP01010044")
  expect_equal(sample_to_patient_id("DRUP01070024T", lookup), "CPCT02070110")
  expect_equal(sample_to_patient_id("DRUP01050008T", lookup), "CPCT02050116")

  expect_equal(sample_to_patient_id("DRUP01010065T", lookup), "CPCT02010639")
  expect_equal(sample_to_patient_id("DRUP01330002T", lookup), "CPCT02330049")
  expect_equal(sample_to_patient_id("DRUP01340004T", lookup), "CPCT02340029")
  expect_equal(sample_to_patient_id("DRUP01340003T", lookup), "CPCT02340014")
  expect_equal(sample_to_patient_id("DRUP01340002T", lookup), "CPCT02340026")
  expect_equal(sample_to_patient_id("DRUP01070008T", lookup), "CPCT02070023")

})

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

