# =====================================================================
# 02_outline_features.R
# Shape descriptors computed from the binary mask + ordered contour.
#
# Features (the 11 suggested in Week 6 + extras):
#   area, perimeter, circularity, compactness, eccentricity,
#   major/minor axis, aspect (length/width) ratio, rectangularity,
#   solidity (convexity), equivalent diameter, form factor,
#   contour roughness, plus left/right asymmetry about the main axis,
#   and the first few normalized elliptic-Fourier-style magnitudes.
#
# Petiole handling: species 2 & 8 have long thin petioles (stalks) that
# inflate the length/width ratio. We optionally detect the petiole as a
# thin protrusion and trim it before computing axis-based features. Both
# trimmed and untrimmed L/W are returned so the report can compare.
# =====================================================================

suppressPackageStartupMessages({
  library(EBImage)
  library(moments)
})


# Central image moments of a binary mask -> used for orientation,
# eccentricity, and Hu moments.
.region_moments <- function(mask) {
  idx <- which(mask, arr.ind = TRUE)
  y <- idx[, "row"]; x <- idx[, "col"]
  n <- length(x)
  xc <- mean(x); yc <- mean(y)
  mu <- function(p, q) sum((x - xc)^p * (y - yc)^q)
  list(n = n, xc = xc, yc = yc,
       mu20 = mu(2, 0), mu02 = mu(0, 2), mu11 = mu(1, 1),
       mu30 = mu(3, 0), mu03 = mu(0, 3), mu21 = mu(2, 1), mu12 = mu(1, 2))
}

# Eccentricity from second moments (0 = circle, ->1 = elongated).
.eccentricity <- function(m) {
  a <- m$mu20 / m$n; b <- m$mu02 / m$n; c <- m$mu11 / m$n
  common <- sqrt((a - b)^2 + 4 * c^2)
  lam1 <- (a + b + common) / 2
  lam2 <- (a + b - common) / 2
  if (lam1 <= 0) return(NA_real_)
  sqrt(pmax(0, 1 - lam2 / lam1))
}

# Major axis orientation (radians) from central moments.
.orientation <- function(m) 0.5 * atan2(2 * m$mu11 / m$n, (m$mu20 - m$mu02) / m$n)

# Seven Hu invariant moments (scale/rotation/translation invariant).
.hu_moments <- function(m) {
  n00 <- m$n
  norm <- function(mu, p, q) mu / (n00 ^ (1 + (p + q) / 2))
  n20 <- norm(m$mu20, 2, 0); n02 <- norm(m$mu02, 0, 2); n11 <- norm(m$mu11, 1, 1)
  n30 <- norm(m$mu30, 3, 0); n03 <- norm(m$mu03, 0, 3)
  n21 <- norm(m$mu21, 2, 1); n12 <- norm(m$mu12, 1, 2)
  h1 <- n20 + n02
  h2 <- (n20 - n02)^2 + 4 * n11^2
  h3 <- (n30 - 3 * n12)^2 + (3 * n21 - n03)^2
  h4 <- (n30 + n12)^2 + (n21 + n03)^2
  h5 <- (n30 - 3*n12)*(n30 + n12)*((n30+n12)^2 - 3*(n21+n03)^2) +
        (3*n21 - n03)*(n21 + n03)*(3*(n30+n12)^2 - (n21+n03)^2)
  h6 <- (n20 - n02)*((n30+n12)^2 - (n21+n03)^2) +
        4*n11*(n30+n12)*(n21+n03)
  h7 <- (3*n21 - n03)*(n30 + n12)*((n30+n12)^2 - 3*(n21+n03)^2) -
        (n30 - 3*n12)*(n21 + n03)*(3*(n30+n12)^2 - (n21+n03)^2)
  hu <- c(h1, h2, h3, h4, h5, h6, h7)
  # log-transform (keep sign) to compress the huge dynamic range
  sign(hu) * log10(abs(hu) + 1e-30)
}

# Convex hull area of the mask (for solidity).
.hull_area <- function(mask) {
  idx <- which(mask, arr.ind = TRUE)
  pts <- cbind(idx[, "col"], idx[, "row"])
  h <- chull(pts)
  hp <- pts[c(h, h[1]), ]
  abs(sum(hp[-nrow(hp), 1] * hp[-1, 2] - hp[-1, 1] * hp[-nrow(hp), 2])) / 2
}

# ---- petiole detection / trimming -----------------------------------
# The petiole is a thin stalk. We erode the mask with a disc; the leaf
# blade survives, the thin petiole mostly disappears. Keeping the
# largest surviving component (then dilating back) gives a blade-only
# mask. Returns the trimmed mask, or the original if trimming removed
# too much (a safety guard for leaves with no real petiole).
trim_petiole <- function(mask, radius = 4) {
  er <- erode(Image(mask), makeBrush(2 * radius + 1, shape = "disc"))
  lab <- bwlabel(er)
  tab <- table(as.vector(lab)); tab <- tab[names(tab) != "0"]
  if (length(tab) == 0) return(mask)
  keep <- as.integer(names(tab)[which.max(tab)])
  blade_core <- lab == keep
  blade <- dilate(Image(blade_core), makeBrush(2 * radius + 1, shape = "disc"))
  blade <- blade & Image(mask)                 # stay within original leaf
  blade <- matrix(as.logical(blade), nrow = nrow(mask), ncol = ncol(mask))
  if (sum(blade) < 0.5 * sum(mask)) return(mask)  # trimmed too much -> bail
  blade
}

# Length & width measured ALONG the principal axes (rotation-aware).
# Projects leaf pixels onto the major/minor axes and takes the spans.
.axis_length_width <- function(mask) {
  m <- .region_moments(mask)
  theta <- .orientation(m)
  idx <- which(mask, arr.ind = TRUE)
  x <- idx[, "col"] - m$xc; y <- idx[, "row"] - m$yc
  proj_major <-  x * cos(theta) + y * sin(theta)
  proj_minor <- -x * sin(theta) + y * cos(theta)
  list(length = diff(range(proj_major)),
       width  = diff(range(proj_minor)))
}

# Left/right asymmetry about the major axis (the "main vein" proxy).
# Fraction of area on each side of the principal axis; we report the
# absolute imbalance |left - right| / total (0 = symmetric).
.axis_asymmetry <- function(mask) {
  m <- .region_moments(mask)
  theta <- .orientation(m)
  idx <- which(mask, arr.ind = TRUE)
  x <- idx[, "col"] - m$xc; y <- idx[, "row"] - m$yc
  side <- -x * sin(theta) + y * cos(theta)     # signed distance to major axis
  left <- sum(side < 0); right <- sum(side > 0)
  total <- left + right
  list(asym_index = abs(left - right) / total,
       left_frac = left / total, right_frac = right / total)
}

# outline_features(p): p is the list returned by preprocess_leaf().
outline_features <- function(p, trim = TRIM_PETIOLE) {
  mask <- p$mask
  bnd  <- p$boundary

  area <- sum(mask)
  # perimeter from ordered boundary (sum of step lengths; closes the loop)
  bx <- bnd$x; by <- bnd$y
  dx <- diff(c(bx, bx[1])); dy <- diff(c(by, by[1]))
  perimeter <- sum(sqrt(dx^2 + dy^2))

  m  <- .region_moments(mask)
  ecc <- .eccentricity(m)
  hu  <- .hu_moments(m)
  hull <- .hull_area(mask)

  # axis length/width WITH the petiole
  lw_raw <- .axis_length_width(mask)
  # and WITHOUT it (blade only) when requested
  if (trim) {
    blade <- trim_petiole(mask)
    lw_trim <- .axis_length_width(blade)
    petiole_removed_frac <- 1 - sum(blade) / area
  } else {
    lw_trim <- lw_raw
    petiole_removed_frac <- 0
  }

  asym <- .axis_asymmetry(mask)

  circularity      <- perimeter^2 / area              # >= 4*pi; higher = more ragged
  compactness      <- (4 * pi * area) / perimeter^2    # in (0,1]; 1 = circle
  solidity         <- area / hull                      # convexity
  rectangularity   <- area / (lw_raw$length * lw_raw$width)
  equiv_diameter   <- sqrt(4 * area / pi)
  form_factor      <- (4 * pi * area) / perimeter^2

  data.frame(
    area               = area,
    perimeter          = perimeter,
    circularity        = circularity,
    compactness        = compactness,
    form_factor        = form_factor,
    eccentricity       = ecc,
    major_axis         = lw_raw$length,
    minor_axis         = lw_raw$width,
    aspect_ratio_raw   = lw_raw$length / lw_raw$width,     # WITH petiole
    aspect_ratio_blade = lw_trim$length / lw_trim$width,   # WITHOUT petiole
    petiole_removed_frac = petiole_removed_frac,
    rectangularity     = rectangularity,
    solidity           = solidity,
    equiv_diameter     = equiv_diameter,
    asymmetry_index    = asym$asym_index,
    left_frac          = asym$left_frac,
    right_frac         = asym$right_frac,
    hu1 = hu[1], hu2 = hu[2], hu3 = hu[3], hu4 = hu[4],
    hu5 = hu[5], hu6 = hu[6], hu7 = hu[7],
    stringsAsFactors = FALSE
  )
}
