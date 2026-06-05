# =====================================================================
# 05_run_all.R
# Loop over every selected species and image, extract all feature
# groups, and write one tidy table to output/features.csv.
#
# Run AFTER sourcing config.R + scripts 01-04. This is the slow step
# (a few minutes for 5 species at 512px). It only needs to run once;
# the report and EDA read the cached CSV.
# =====================================================================

stopifnot(exists("DATA_DIR"), exists("SPECIES_IDS"))

# build the file path for species id, image index n (1..75)
.leaf_path <- function(id, n) {
  folder <- gsub("\\{id\\}", id, FOLDER_TEMPLATE)
  fname  <- FILE_TEMPLATE
  fname  <- gsub("\\{id\\}", id, fname)
  fname  <- gsub("\\{n\\}", sprintf("%03d", n), fname)
  file.path(DATA_DIR, folder, fname)
}

process_one <- function(id, n) {
  path <- .leaf_path(id, n)
  if (!file.exists(path)) return(NULL)
  p <- tryCatch(preprocess_leaf(path), error = function(e) NULL)
  if (is.null(p)) {
    warning("preprocess failed: ", path); return(NULL)
  }
  out <- tryCatch(
    cbind(
      data.frame(species = id,
                 species_name = unname(SPECIES_NAMES[as.character(id)]),
                 image = n, file = basename(path),
                 stringsAsFactors = FALSE),
      outline_features(p),
      color_features(p),
      vein_features(p)
    ),
    error = function(e) { warning("features failed: ", path, " : ", conditionMessage(e)); NULL }
  )
  out
}

# ---------------------------------------------------------------------
# preflight(): verify the data is reachable BEFORE the slow loop. Stops
# with an actionable message if the folder/file templates don't match
# what's actually on disk -- this is the usual cause of an empty
# features.csv (and the downstream "object 'species' not found" error).
# ---------------------------------------------------------------------
preflight <- function() {
  # 1. Does the dataset directory exist at all?
  if (!dir.exists(DATA_DIR)) {
    stop("DATA_DIR does not exist:\n  ", DATA_DIR,
         "\nFix DATA_DIR in R/config.R (use forward slashes).", call. = FALSE)
  }

  # 2. Check the first expected path for each requested species.
  missing_species <- c()
  example_path <- NULL
  for (id in SPECIES_IDS) {
    p1 <- .leaf_path(id, 1)
    if (is.null(example_path)) example_path <- p1
    if (!file.exists(p1)) missing_species <- c(missing_species, id)
  }

  if (length(missing_species) > 0) {
    # Show what IS on disk so the templates can be corrected.
    dirs <- list.dirs(DATA_DIR, recursive = FALSE, full.names = FALSE)
    dir_sample <- if (length(dirs)) paste(utils::head(dirs, 5), collapse = ", ") else "(none)"
    first_dir <- if (length(dirs)) file.path(DATA_DIR, dirs[1]) else DATA_DIR
    files <- list.files(first_dir)
    file_sample <- if (length(files)) paste(utils::head(files, 5), collapse = ", ") else "(none)"

    stop(
      "Could not find images for species: ", paste(missing_species, collapse = ", "), "\n",
      "Tried path like:\n  ", example_path, "\n\n",
      "What's actually in DATA_DIR:\n",
      "  subfolders : ", dir_sample, "\n",
      "  files in '", if (length(dirs)) dirs[1] else "<DATA_DIR>", "' : ", file_sample, "\n\n",
      "Fix FOLDER_TEMPLATE / FILE_TEMPLATE in R/config.R so they match.\n",
      "  Current FOLDER_TEMPLATE = '", FOLDER_TEMPLATE, "'\n",
      "  Current FILE_TEMPLATE   = '", FILE_TEMPLATE, "'\n",
      "  ({id} = species number, {n} = 3-digit image index)",
      call. = FALSE
    )
  }

  message("preflight OK: found first image for all ",
          length(SPECIES_IDS), " species.")
  invisible(TRUE)
}

run_all <- function() {
  preflight()   # stop early with a clear message if paths are wrong

  rows <- list(); k <- 0
  total <- length(SPECIES_IDS) * N_PER_SPECIES
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  step <- 0
  for (id in SPECIES_IDS) {
    for (n in seq_len(N_PER_SPECIES)) {
      step <- step + 1; setTxtProgressBar(pb, step)
      r <- process_one(id, n)
      if (!is.null(r)) { k <- k + 1; rows[[k]] <- r }
    }
  }
  close(pb)

  # Stop loudly if nothing (or almost nothing) came through, rather than
  # writing an empty CSV that breaks script 06 later.
  if (length(rows) == 0) {
    stop("No images were processed successfully. ",
         "Paths matched in preflight but every image failed to read or ",
         "feature-extract. Re-run with one image to see the warning:\n",
         "  source('R/01_preprocess.R'); ",
         "p <- preprocess_leaf(.leaf_path(SPECIES_IDS[1], 1))",
         call. = FALSE)
  }
  expected <- length(SPECIES_IDS) * N_PER_SPECIES
  if (length(rows) < 0.5 * expected) {
    warning("Only ", length(rows), " of ", expected,
            " images processed. Check earlier warnings for failures.")
  }

  feats <- do.call(rbind, rows)
  out_path <- file.path(OUTPUT_DIR, "features.csv")
  write.csv(feats, out_path, row.names = FALSE)
  message("\nWrote ", nrow(feats), " rows to ", out_path)
  invisible(feats)
}

# Execute when sourced:
features <- run_all()
