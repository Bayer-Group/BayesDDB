test_that("elastic_fun returns numeric vector of same length as x", {
  ppp <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  out <- elastic_fun(ppp, map_a = 1, map_b = 1)
  expect_type(out, "double")
  expect_length(out, length(ppp))
})

test_that("elastic_fun gives 1 at x = 0.5 for any map_a and map_b", {
  for (map_a in c(0.5, 1, 5, 20)) {
    for (map_b in c(0.5, 1, 3)) {
      expect_equal(elastic_fun(0.5, map_a = map_a, map_b = map_b), 1)
    }
  }
})

test_that("elastic_fun output is in [0, 1] and symmetric around 0.5", {
  x <- seq(0, 1, by = 0.05)
  for (map_b in c(0.5, 1, 3)) {
    out <- elastic_fun(x, map_a = 5, map_b = map_b)
    expect_true(all(out >= 0 & out <= 1))
    expect_equal(
      out,
      elastic_fun(1 - x, map_a = 5, map_b = map_b),
      tolerance = 1e-6
    )
  }
})

test_that("elastic_fun implements the two-parameter arctangent form", {
  ppp <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  expect_equal(
    elastic_fun(ppp, map_a = 2, map_b = 3),
    atan(2 * sin(ppp * pi))^3 / atan(2)^3
  )
  expect_equal(
    elastic_fun(ppp, map_a = 1, map_b = 1),
    atan(sin(ppp * pi)) / atan(1)
  )
})

test_that("elastic_fun reproduces reference values for map_a = 1, map_b = 3", {
  out <- elastic_fun(c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), map_a = 1, map_b = 3)
  expect_equal(round(out, 3), c(0, 0.056, 0.481, 1, 0.481, 0.056, 0))
})

test_that("larger map_a is more lenient towards outcome differences", {
  ppp <- c(0.05, 0.1, 0.2, 0.8, 0.9, 0.95)
  strict <- elastic_fun(ppp, map_a = 0.5, map_b = 1)
  mild <- elastic_fun(ppp, map_a = 5, map_b = 1)
  lenient <- elastic_fun(ppp, map_a = 20, map_b = 1)
  expect_true(all(strict < mild))
  expect_true(all(mild < lenient))
})

test_that("larger map_b discounts more away from x = 0.5", {
  ppp <- c(0.1, 0.2, 0.8, 0.9)
  expect_true(all(
    elastic_fun(ppp, map_a = 1, map_b = 3) <
      elastic_fun(ppp, map_a = 1, map_b = 1)
  ))
})

test_that("elastic_fun errors on invalid input", {
  expect_error(elastic_fun(NA, 1, 1), "NAs are not allowed")
  expect_error(elastic_fun(c(0.5, NA), 1, 1), "NAs are not allowed")
  expect_error(elastic_fun(-0.1, 1, 1), "out of the interval")
  expect_error(elastic_fun(1.5, 1, 1), "out of the interval")
  expect_error(elastic_fun(0.5, map_a = 0, map_b = 1), "map_a must be positive")
  expect_error(
    elastic_fun(0.5, map_a = -1, map_b = 1),
    "map_a must be positive"
  )
  expect_error(elastic_fun(0.5, map_a = 1, map_b = 0), "map_b must be positive")
  expect_error(
    elastic_fun(0.5, map_a = 1, map_b = -1),
    "map_b must be positive"
  )
})

test_that("elastic_fun rejects non-numeric x and malformed tuning parameters", {
  expect_error(elastic_fun("0.5", 1, 1), "must be numeric")
  expect_error(elastic_fun(NULL, 1, 1), "must be numeric")
  for (bad in list(NA, NA_real_, NULL, "1", c(1, 2), Inf, NaN)) {
    expect_error(
      elastic_fun(0.5, map_a = bad, map_b = 1),
      "map_a must be positive"
    )
    expect_error(
      elastic_fun(0.5, map_a = 1, map_b = bad),
      "map_b must be positive"
    )
  }
  expect_length(elastic_fun(numeric(0), 1, 1), 0)
})
