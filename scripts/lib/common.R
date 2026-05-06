# Common utilities shared by figure scripts

get_repo_root <- function() {
  # Prefer running from repo root; fall back to script location if possible.
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  # If user runs from code/ or code/figures/, walk up.
  if (basename(wd) %in% c("figures", "lib")) wd <- dirname(wd)
  if (basename(wd) == "code") wd <- dirname(wd)

  wd
}

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nInstall them (e.g. install.packages(...)) and re-run."
    )
  }
}

save_pdf <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "pdf"
  )
}

add_significance <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  "ns"
}

cohens_d_with_ci <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2 || n2 < 2) {
    return(list(d = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, d_str = "N/A"))
  }
  mean1 <- mean(x, na.rm = TRUE)
  mean2 <- mean(y, na.rm = TRUE)
  var1 <- stats::var(x, na.rm = TRUE)
  var2 <- stats::var(y, na.rm = TRUE)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
  if (pooled_sd == 0 || is.na(pooled_sd)) {
    return(list(d = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, d_str = "N/A"))
  }
  d <- (mean1 - mean2) / pooled_sd
  se_d <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2)))
  ci_lower <- d - 1.96 * se_d
  ci_upper <- d + 1.96 * se_d
  d_str <- paste0(round(d, 2), " [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]")
  list(d = d, ci_lower = ci_lower, ci_upper = ci_upper, d_str = d_str)
}

clr_transform <- function(x) {
  # Add small pseudocount to avoid zeros
  x <- x + 0.001
  geom_mean <- exp(rowMeans(log(x)))
  log(x / geom_mean)
}

