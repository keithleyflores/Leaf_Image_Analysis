# =====================================================================
# 01_preprocess.R
# Read a leaf image and turn it into the objects every downstream step
# needs: the color channels, a clean binary leaf mask, and the ordered
# boundary contour.
#
# Pipeline (matches the course discussion files):
#   color -> grayscale -> threshold -> binary mask
#         -> keep largest connected component -> fill holes
#         -> contour (8-neighbourhood boundary)
#
# NOTE on image reading: TIFFs are read with the `tiff` package (which
# links libtiff directly) instead of ImageMagick/magick. Some ImageMagick
# builds -- especially the one bundled with conda -- ship without the TIFF
# delegate and fail on these scans with:
#     MustSpecifyImageSize ... ReadRAWImage
# Reading via tiff::readTIFF sidesteps that entirely.
# =====================================================================

# ---- dependencies ---------------------------------------------------
# Install only when missing (don't reinstall on every source()).
if (!requireNamespace("EBImage", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("EBImage", update = FALSE, ask = FALSE)
}
if (!requireNamespace("tiff", quietly = TRUE)) install.packages("tiff")

suppressPackageStartupMessages({
  library(tiff)
  library(EBImage)
})

# ---------------------------------------------------------------------
# read_leaf(): load a TIFF as an H x W x 3 array (values 0-1) via libtiff,
# optionally downsizing the long side to resize_max. No ImageMagick.
#   - grayscale TIFFs (H x W) are promoted to 3 channels
#   - alpha (4th channel) is dropped
#   - resizing uses EBImage::resize
# ---------------------------------------------------------------------
read_leaf <- function(path, resize_max = RESIZE_MAX_PX) {
  arr <- tiff::readTIFF(path, native = FALSE)   # H x W (gray) or H x W x C, 0-1

  # promote grayscale to 3 channels so get_channels() always has R/G/B
  if (length(dim(arr)) == 2) arr <- array(arr, dim = c(dim(arr), 3))
  # drop alpha / extra channels, keep first 3
  if (dim(arr)[3] > 3) arr <- arr[, , 1:3, drop = FALSE]
  # single-channel stored as H x W x 1 -> replicate to 3
  if (dim(arr)[3] == 1) arr <- array(arr[, , 1], dim = c(dim(arr)[1:2], 3))

  if (!is.na(resize_max)) {
    h <- dim(arr)[1]; w <- dim(arr)[2]
    long_side <- max(w, h)
    if (long_side > resize_max) {
      scale <- resize_max / long_side
      # EBImage uses (x=cols, y=rows); aperm to its layout, resize, aperm back
      img <- EBImage::Image(aperm(arr, c(2, 1, 3)), colormode = "Color")
      img <- EBImage::resize(img, w = round(w * scale), h = round(h * scale))
      arr <- aperm(EBImage::imageData(img), c(2, 1, 3))
    }
  }
  arr   # H x W x 3, values in 0-1
}

# ---------------------------------------------------------------------
# get_channels(): pull R, G, B matrices (0-255) and a grayscale matrix
# (0-1) from the H x W x 3 array. Rows = image rows (y), cols = cols (x).
# Grayscale uses the standard Rec. 601 luma weights.
# ---------------------------------------------------------------------
get_channels <- function(arr) {
  h <- dim(arr)[1]; w <- dim(arr)[2]
  R <- round(arr[, , 1] * 255)
  G <- round(arr[, , 2] * 255)
  B <- round(arr[, , 3] * 255)
  gray <- 0.299 * arr[, , 1] + 0.587 * arr[, , 2] + 0.114 * arr[, , 3]  # 0-1
  list(R = R, G = G, B = B, gray = gray, width = w, height = h)
}

# ---------------------------------------------------------------------
# leaf_mask(): binary leaf region from the grayscale matrix.
# Leaf is dark (gray <= threshold); keep the single largest blob and
# fill interior holes so vein gaps / specular spots don't punch holes.
# ---------------------------------------------------------------------
leaf_mask <- function(gray, threshold = GRAY_THRESHOLD) {
  raw <- gray <= threshold                       # logical matrix, TRUE = leaf

  lab <- bwlabel(Image(raw))                     # label connected components
  tab <- table(as.vector(lab))
  tab <- tab[names(tab) != "0"]                  # drop background label 0
  if (length(tab) == 0) stop("No leaf region found -- check threshold.")
  largest <- as.integer(names(tab)[which.max(tab)])

  mask <- lab == largest
  mask <- fillHull(mask)                         # fill holes inside the leaf
  matrix(as.logical(mask), nrow = nrow(gray), ncol = ncol(gray))
}

# ---------------------------------------------------------------------
# mask_contour(): boundary pixels of the mask using an 8-neighbourhood.
# A leaf pixel is on the boundary if eroding the mask removes it.
# Returns BOTH a logical edge matrix and ordered (x, y) boundary coords.
# ---------------------------------------------------------------------
mask_contour <- function(mask) {
  m <- Image(mask)
  edge <- mask & !erode(m, makeBrush(3, shape = "box"))
  edge <- matrix(as.logical(edge), nrow = nrow(mask), ncol = ncol(mask))

  # ordered boundary via EBImage's ocontour (walks the perimeter in order,
  # which the Fourier / curvature features need).
  oc <- ocontour(Image(mask))
  if (length(oc) == 0) {
    coords <- which(edge, arr.ind = TRUE)
    ordered <- data.frame(x = coords[, "col"], y = coords[, "row"])
  } else {
    biggest <- oc[[which.max(vapply(oc, nrow, integer(1)))]]
    # ocontour returns 0-based (col,row) == (x,y); make 1-based
    ordered <- data.frame(x = biggest[, 1] + 1, y = biggest[, 2] + 1)
  }
  list(edge = edge, boundary = ordered)
}

# ---------------------------------------------------------------------
# preprocess_leaf(): one call that returns everything downstream needs.
# `img` now holds the H x W x 3 array (was a magick object previously);
# nothing downstream reads p$img, but it's kept for the visual check.
# ---------------------------------------------------------------------
preprocess_leaf <- function(path, resize_max = RESIZE_MAX_PX,
                            threshold = GRAY_THRESHOLD) {
  arr  <- read_leaf(path, resize_max)
  ch   <- get_channels(arr)
  mask <- leaf_mask(ch$gray, threshold)
  ct   <- mask_contour(mask)
  list(
    path = path, img = arr,
    R = ch$R, G = ch$G, B = ch$B, gray = ch$gray,
    width = ch$width, height = ch$height,
    mask = mask, edge = ct$edge, boundary = ct$boundary
  )
}

# ---- quick visual check (run interactively on one image) ------------
# Wrapped in if (FALSE) so that source()-ing this file only DEFINES the
# functions and does not try to read an image. To run the check, set the
# condition to TRUE (or paste the body at the console) after config.R is
# loaded so DATA_DIR exists.
if (FALSE) {
  p <- preprocess_leaf(file.path(DATA_DIR, "leaf1", "l1nr001.tif"))
  par(mfrow = c(1, 3))
  image(t(p$gray[nrow(p$gray):1, ]), col = gray.colors(256), main = "gray")
  image(t(p$mask[nrow(p$mask):1, ]), col = c("white", "black"), main = "mask")
  image(t(p$edge[nrow(p$edge):1, ]), col = c("white", "black"), main = "edge")
}
