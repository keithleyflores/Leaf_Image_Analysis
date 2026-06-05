# =====================================================================
# 04_vein_features.R
# Vein extraction by gray-scale morphology, following the idea in
# Zheng & Wang (2010), "Leaf Vein Extraction Based on Gray-Scale
# Morphology".
#
# Idea: veins are thin structures slightly different in intensity from
# the surrounding blade. A morphological white/black top-hat at several
# small structuring-element sizes responds strongly to those thin
# structures while suppressing the smooth blade background. We threshold
# the combined top-hat response, restrict to the leaf interior, and
# summarize vein density and total vein length.
#
# These are texture/structure features -- "inner structure" in the
# Soderkvist taxonomy -- and are returned per image for the report.
# =====================================================================

suppressPackageStartupMessages({
  library(EBImage)
})

# Multi-scale black top-hat: closing(img) - img. Veins on a leaf scan
# are typically slightly DARKER than the blade, so the black top-hat
# (a.k.a. bottom-hat) highlights them. We also add the white top-hat to
# be robust to lighter veins, and combine.
.tophat_response <- function(gray_in_leaf, radii = c(1, 2, 3)) {
  img <- Image(gray_in_leaf)
  resp <- matrix(0, nrow = nrow(gray_in_leaf), ncol = ncol(gray_in_leaf))
  for (r in radii) {
    se <- makeBrush(2 * r + 1, shape = "disc")
    black_th <- closing(img, se) - img          # dark thin structures
    white_th <- img - opening(img, se)           # light thin structures
    resp <- resp + as.matrix(black_th) + as.matrix(white_th)
  }
  resp / length(radii)
}

# .skeletonize(): use EBImage::thinImage when the installed version
# provides it; otherwise fall back to a distance-map ridge proxy. The
# fallback is decided at call time (robust across EBImage versions).
.skeletonize <- function(x) {
  if ("thinImage" %in% getNamespaceExports("EBImage")) {
    return(EBImage::thinImage(x))
  }
  # fallback: ridge = pixels whose distance-transform value is a local
  # maximum vs. left/up neighbours (cheap centreline proxy).
  dmat <- as.matrix(distmap(x))
  up   <- rbind(0, dmat[-nrow(dmat), ])
  left <- cbind(0, dmat[, -ncol(dmat)])
  ridge <- dmat > 0 & dmat >= up & dmat >= left
  Image(ridge)
}

# vein_features(p): operate inside the leaf only, ignore the boundary
# ring (erode the mask a little so the leaf edge isn't counted as vein).
vein_features <- function(p, radii = c(1, 2, 3)) {
  inner <- erode(Image(p$mask), makeBrush(5, shape = "disc"))
  inner <- matrix(as.logical(inner), nrow = nrow(p$mask), ncol = ncol(p$mask))

  gray <- p$gray
  gray[!inner] <- NA                              # restrict to interior

  # fill NAs with the interior mean so morphology doesn't see hard edges
  fill_val <- mean(gray[inner], na.rm = TRUE)
  g_filled <- gray; g_filled[!inner] <- fill_val

  resp <- .tophat_response(g_filled, radii)
  resp[!inner] <- 0

  # adaptive threshold: mean + k*sd of the interior response
  rv <- resp[inner]
  thr <- mean(rv) + 1.0 * sd(rv)
  vein <- resp > thr & inner

  # thin to skeleton so "length" is meaningful, then measure
  skel_img <- .skeletonize(Image(vein > 0))
  skel <- matrix(as.logical(skel_img), nrow = nrow(vein), ncol = ncol(vein))

  leaf_area  <- sum(inner)
  vein_area  <- sum(vein)
  vein_len   <- sum(skel)                         # skeleton pixel count ~ length

  data.frame(
    vein_density   = vein_area / leaf_area,       # fraction of blade that is vein
    vein_length    = vein_len,                    # total skeleton length (px)
    vein_len_per_area = vein_len / leaf_area,     # length normalized by leaf size
    vein_response_mean = mean(rv),
    vein_response_sd   = sd(rv),
    stringsAsFactors = FALSE
  )
}

