# =====================================================================
# 03_color_features.R
# Color descriptors computed ONLY over leaf pixels (mask == TRUE).
# Covers the three systems from the course material:
#   RGB   -- raw channel statistics
#   HSV   -- hue / saturation / value (closer to human description)
#   CIELAB-- perceptually-uniform lightness (L*) and color (a*, b*)
#
# Also: per-region color (center / edge / tip / base) and a uniformity
# measure, to answer "is color uniform within a leaf or does it vary?"
# =====================================================================

suppressPackageStartupMessages({
  library(grDevices)   # convertColor, rgb2hsv
})

.summ <- function(v, prefix) {
  v <- v[is.finite(v)]
  setNames(
    list(mean(v), sd(v), median(v),
         as.numeric(quantile(v, 0.25)), as.numeric(quantile(v, 0.75))),
    paste0(prefix, c("_mean", "_sd", "_median", "_q25", "_q75"))
  )
}

# Convert leaf RGB (0-255) to CIELAB. convertColor expects 0-1 sRGB.
.rgb_to_lab <- function(R, G, B) {
  srgb <- cbind(R, G, B) / 255
  lab <- convertColor(srgb, from = "sRGB", to = "Lab")
  list(L = lab[, 1], a = lab[, 2], b = lab[, 3])
}

# Assign each leaf pixel to a coarse region: center / edge / tip / base.
# tip = top 25% of leaf rows, base = bottom 25%, then within the middle
# band we split edge vs center by distance to the boundary.
.region_labels <- function(p) {
  idx <- which(p$mask, arr.ind = TRUE)
  y <- idx[, "row"]; x <- idx[, "col"]
  yr <- range(y); span <- diff(yr)
  region <- rep("center", length(x))
  region[y <= yr[1] + 0.25 * span] <- "tip"     # image top = leaf tip (convention)
  region[y >= yr[2] - 0.25 * span] <- "base"

  # edge vs center within the non-tip/base pixels, via distance map
  dm <- EBImage::distmap(EBImage::Image(p$mask))
  dm <- matrix(as.numeric(EBImage::imageData(dm)),
               nrow = nrow(p$mask), ncol = ncol(p$mask))
  dvals <- dm[cbind(y, x)]                       # [row, col] indexing
  thr <- as.numeric(quantile(dvals, 0.5))       # outer half = edge
  middle <- region == "center"
  region[middle & dvals <= thr] <- "edge"
  region
}

color_features <- function(p) {
  idx <- which(p$mask)
  R <- p$R[idx]; G <- p$G[idx]; B <- p$B[idx]

  # HSV (rgb2hsv wants a 3 x n matrix of 0-255 values)
  hsv <- grDevices::rgb2hsv(rbind(R, G, B), maxColorValue = 255)
  H <- hsv["h", ]; S <- hsv["s", ]; V <- hsv["v", ]

  # CIELAB
  lab <- .rgb_to_lab(R, G, B)

  feats <- c(
    .summ(R, "R"), .summ(G, "G"), .summ(B, "B"),
    .summ(H, "H"), .summ(S, "S"), .summ(V, "V"),
    .summ(lab$L, "L"), .summ(lab$a, "a"), .summ(lab$b, "b")
  )

  # Greenness proxies that are robust to brightness:
  #   ExG (excess green) and the green channel fraction.
  exg <- (2 * G - R - B)
  green_frac <- G / (R + G + B + 1e-9)
  feats <- c(feats,
             ExG_mean = mean(exg), ExG_sd = sd(exg),
             green_frac_mean = mean(green_frac))

  # Proportion of pixels in coarse color bins (greener / yellower / darker)
  greener  <- mean(G > R & G > B)
  yellower <- mean(R > B & G > B & abs(R - G) < 40)
  darker   <- mean(V < 0.3)
  feats <- c(feats,
             prop_greener = greener, prop_yellower = yellower,
             prop_darker = darker)

  # Within-leaf color uniformity: coefficient of variation of Hue & Value
  feats <- c(feats,
             hue_cv   = sd(H) / (mean(H) + 1e-9),
             value_cv = sd(V) / (mean(V) + 1e-9))

  # Per-region mean greenness (center/edge/tip/base) + tip-vs-base shift
  reg <- .region_labels(p)
  region_green <- tapply(green_frac, reg, mean)
  get_reg <- function(nm) if (nm %in% names(region_green)) region_green[[nm]] else NA_real_
  feats <- c(feats,
             green_center = get_reg("center"),
             green_edge   = get_reg("edge"),
             green_tip    = get_reg("tip"),
             green_base   = get_reg("base"),
             green_tip_minus_base = get_reg("tip") - get_reg("base"))

  as.data.frame(feats, stringsAsFactors = FALSE)
}
