
as_longDF_and_rm_NA <- function(obj) {
  # Convert to long data frame format:
  df <- hyperSpec::as.long.df(obj, rownames = TRUE, na.rm = FALSE)
  # Remove NA values:
  df <- df[!is.na(df$spc), , drop = FALSE]
  return(df)
}

extract_wl_ranges <- function(obj, wl.range, format = "%g") {

  wl.range <- c(wl.range)

  if (length(wl.range) != 1) {
    stop("`wl.range` must contain only one range.",
      "Now the ranges are:\n ",
      paste(wl.range, collapse = ", ")
    )
  }
  # extract range
  x <- (obj[, , wl.range])

  # Make label for range
  range_label <- wl(x) %>% range %>% sprintf(format, .) %>%
    paste(collapse = "-")

  # Create variable `.wl.range` which wil be used for facetting
  x$.wl.range <- as.factor(range_label)
  return(x)

}


ggplot.hyperSpec <- function(data, mapping = aes(), ..., wl.range = NULL,
  format = "%g", environment) {

  chk.hy(data)
  validObject(data)

  force(format)

  mapping <- modifyList(
    aes_string(
      x = ".wavelength",
      y = "spc",
      group = ".rownames"),
    mapping)


  # If `wl.range` is provided a subset of object is taken
  if (!is.null(wl.range))   {
    wl.range <- c(wl.range) # Convert to list
    obj_list <- lapply(wl.range, extract_wl_ranges, obj = data, format = format)

    # Transform to data frame
    DF <- lapply(obj_list, as_longDF_and_rm_NA) %>% Reduce(rbind, .)

  } else {
    # if `wl.range` is missing whole dataset is used

    # Transform to data frame
    DF <- as_longDF_and_rm_NA(data)
  }

  # Make a ggplot
  p <- ggplot2::ggplot(DF, mapping = mapping)

  # Modify x and y axis labels in ggplot
  p <- p +
    xlab(hyperSpec::labels(data, ".wavelength")) +
    ylab(hyperSpec::labels(data, "spc"))

  # Use facetting if two or more wl.ranges are selected
  if (!is.null(wl.range) & length(wl.range) > 1) {
    p <- p + facet_wl()
  }

  return(p)
}
