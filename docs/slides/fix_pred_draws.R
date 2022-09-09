library(tidyverse)

fix_draws = function(object, newdata, ..., func = tidybayes::predicted_draws) {
  draws = func(object, newdata, ...)

  n = names(draws)

  full_join(
    draws %>%
      select(-.chain, -.iteration),
    tidybayes::tidy_draws(b) %>%
      select(.chain, .iteration, .draw),
    by = ".draw"
  ) %>%
    select(all_of(n)) %>%
    ungroup()
}

predicted_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::predicted_draws)
}

epred_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::epred_draws)
}
