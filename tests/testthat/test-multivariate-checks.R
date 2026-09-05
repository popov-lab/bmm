test_that("formula_re_frame() extracts grouping variables with their |ID| labels", {
  re <- formula_re_frame(
    bmf(kappa ~ 1 + (1 | p | id), thetat ~ cond + (cond | q | gr(session, by = grp)))
  )
  expect_equal(re$group, c("id", "session"))
  expect_equal(re$id, c("p", "q"))

  unlabeled <- formula_re_frame(bmf(kappa ~ 1 + (1 | id), thetat ~ 1 + (1 || id)))
  expect_equal(unlabeled$group, "id")
  expect_true(is.na(unlabeled$id))

  expect_equal(nrow(formula_re_frame(bmf(kappa ~ 1, mu1 = 0))), 0)
})

test_that("a missing shared random-effects ID triggers a warning", {
  d1 <- data.frame(id = factor(1:10), y1 = rnorm(10))
  d2 <- data.frame(id = factor(1:10), y2 = rnorm(10))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | id)), family = gaussian(), data = d2)
  expect_warning(configure_fit(joint), "No random-effects ID is shared")

  # the same label on different grouping variables does not tie the components
  d2$session <- factor(1:10)
  joint_groups <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | session)), family = gaussian(), data = d2)
  expect_warning(configure_fit(joint_groups), "No random-effects ID is shared")

  shared <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_silent(configure_fit(shared))
})

test_that("components with no shared grouping values trigger a warning", {
  d1 <- data.frame(id = factor(1:10), y1 = rnorm(10))
  d2 <- data.frame(id = factor(paste0("P", 1:10)), y2 = rnorm(10))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_warning(
    configure_fit(joint),
    "values of the grouping variable 'id' in common"
  )
})

test_that("partially overlapping grouping values are reported", {
  d1 <- data.frame(id = factor(1:6), y1 = rnorm(6))
  d2 <- data.frame(id = factor(4:9), y2 = rnorm(6))
  joint <- bmm_component(bmf(y1 ~ 1 + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ 1 + (1 | p | id)), family = gaussian(), data = d2)
  expect_message(configure_fit(joint), "3 of 9 unique values of 'id'")
})

test_that("factor predictors with different levels across components warn", {
  d1 <- data.frame(
    id = factor(rep(1:5, 2)), y1 = rnorm(10),
    cond = factor(rep(c("a", "b"), 5))
  )
  d2 <- data.frame(
    id = factor(rep(1:5, 2)), y2 = rnorm(10),
    cond = factor(rep(c("a", "c"), 5))
  )
  joint <- bmm_component(bmf(y1 ~ cond + (1 | p | id)), family = gaussian(), data = d1) +
    bmm_component(bmf(y2 ~ cond + (1 | p | id)), family = gaussian(), data = d2)
  expect_warning(configure_fit(joint), "different levels across")
})
