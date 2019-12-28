# Functions ------------------------------------------------------------------

# set standardized color palettes ============================================
seq.palette <- colorRampPalette(c("white", "dark green"), space = "Lab")

YG.palette <- function(n = 20) {
  rgb(colorRamp(c(
    "#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476",
    "#41AB5D", "#238B45", "#006D2C", "#00441B"
  ), space = "Lab")
  # was: brewer.pal (9, "Greens")
  (seq(1 / 3, 1, length.out = n)^2), maxColorValue = 255)
}


div.palette <- colorRampPalette(c(
  "#00008B", "#351C96", "#5235A2", "#6A4CAE", "#8164BA", "#967CC5",
  "#AC95D1", "#C1AFDC", "#D5C9E8", "#E0E3E3", "#F8F8B0", "#F7E6C2",
  "#EFCFC6", "#E6B7AB", "#DCA091", "#D08977", "#C4725E", "#B75B46",
  "#A9432F", "#9A2919", "#8B0000"
), space = "Lab")


# Packages ===================================================================
pkg_suggests <- function(...) {
  strsplit(packageDescription(..., fields = "Suggests"), ",\\s*")[[1]]
}

pkg_exists <- function(pkg = stop("package name needed"), lib.loc = NULL) {
  dir <- sapply(pkg, function(p) system.file(package = p, lib.loc = lib.loc))
  nzchar(dir) > 0L
}

is_basepkg <- function(pkg) {
  pkg_exists(pkg) && grepl("^base$", packageDescription(pkg, fields = "Priority"))
}

pkg_or_base <- function(pkg) {
  pkg[sapply(pkg, is_basepkg)] <- "base"

  pkg
}

citation_or_file <- function(pkg, svd.cit = sprintf("%s.CITATION", pkg)) {
  if (pkg_exists(pkg)) {
    citation(pkg)
  } else if (file.exists(svd.cit)) {
    readCitationFile(file = svd.cit)
  } else {
    NULL
  }
}

make_cite_keys <- function(pkg, entries) {
  pkg <- pkg_or_base(pkg)

  if (!pkg_exists(pkg)) {
    return(pkg)
  }

  if (missing(entries)) {
    entries <- citation_or_file(pkg)
  }

  keys <- sapply(unclass(entries), attr, "key")
  noname <- which(sapply(keys, is.null))
  if (length(keys) == 1L && noname == 1L) {
    keys <- pkg
  } else {
    for (i in noname) {
      keys[[i]] <- paste(pkg, i, sep = ".")
    }
  }
  keys <- make.unique(unlist(keys))
  keys
}

citation_with_key <- function(pkg = "base") {
  pkg  <- pkg_or_base(pkg)
  tmp  <- citation_or_file(pkg)
  keys <- make_cite_keys(pkg, tmp)
  for (entry in seq_along(tmp)) {
    tmp[entry]$"key" <- keys[[entry]]
  }
  tmp
}

cite_pkg <- function(p, entries) {
  paste0("[@", paste(make_cite_keys(p, entries), collapse = ";"), "]")
}
# cite.pkg <- function(p, entries, citefun = "cite") {
#   paste("\\\\", citefun, "{", paste(make_cite_keys(p, entries), collapse = ", "), "}", sep = "")
# }

make_bib <- function(..., file = NULL) {
  pkg <- c(...)

  if (length(pkg) == 0L) {
    pkg <- loadedNamespaces()

    pkg <- unique(pkg_or_base(pkg))
  }

  l <- lapply(pkg, citation_with_key)
  l <- do.call("c", l[!sapply(l, is.null)])

  if (!is.null(file)) {
    if (is.null(l)) {
      cat(NULL, file = file)
    } # touches file
    else {
      cat(toBibtex(l), file = file, sep = "\n")
    }
  }

  invisible(l)
}
