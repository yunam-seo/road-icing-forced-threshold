# Defines the analysis winter window shared by every script that reads the 1-min RWIS data.
# Only winters with data from the start of the freezing season (11 November) are kept, so 2022/23 is excluded.
# Provides filter_winter() for tables with a datetime column or separate year and month columns.

WINTER_MIN <- 2023L

winter_label <- function(y, m) {
  return(as.integer(ifelse(m >= 9L, y, y - 1L)))
}

winter_of_datetime <- function(x) {
  x <- as.numeric(x)
  return(winter_label(as.integer(x %/% 1e8), as.integer((x %/% 1e6) %% 100)))
}

filter_winter <- function(d, col = "datetime", ym = NULL) {
  w <- if (is.null(ym)) winter_of_datetime(d[[col]]) else winter_label(ym$y, ym$m)
  return(d[w >= WINTER_MIN])
}
