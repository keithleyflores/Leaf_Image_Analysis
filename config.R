# =====================================================================
# config.R  --  Paths and run settings.
# =====================================================================

DATA_DIR <- "/Users/keithleyflores/Downloads/leaf_project/leafs"   

OUTPUT_DIR <- "output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# These IDs are the folder numbers 1..15. Species 2 and 8 have long
# petioles -- nice to include at least one of them to show the
# petiole-handling matters, but not required.
SPECIES_IDS <- c(1, 2, 4, 8, 13)

# Human-readable names (from Soderkvist 2001). Used in plots/tables.
SPECIES_NAMES <- c(
  "1"  = "Ulmus carpinifolia",
  "2"  = "Acer platanoides",
  "3"  = "Ulmus glabra",
  "4"  = "Quercus robur",
  "5"  = "Alnus incana",
  "6"  = "Tilia (6)",
  "7"  = "Salix fragilis",
  "8"  = "Populus tremula",
  "9"  = "Corylus avellana",
  "10" = "Sorbus aucuparia",
  "11" = "Prunus padus",
  "12" = "Tilia (12)",
  "13" = "Populus",
  "14" = "Sorbus hybrida",
  "15" = "Fagus silvatica"
)

# Species known to have long petioles (course clarification note).
LONG_PETIOLE_SPECIES <- c(2, 8)

# How image files are named inside each species folder. The default
# matches the dataset convention l{ID}nr{NNN}.tif (e.g. l1nr001.tif).
# {id} is the species number, {n} the 3-digit image index.
FOLDER_TEMPLATE <- "leaf{id}"           # subfolder per species
FILE_TEMPLATE   <- "l{id}nr{n}.tif"     # file inside it
N_PER_SPECIES   <- 75                   # images per species

# Resize longest side to this many px before processing (speeds things
# up a lot; set to NA to use full resolution). 512 is a good balance.
RESIZE_MAX_PX <- 512L

# Grayscale threshold for leaf/background split (course uses 200 on a
# 0-255 scale -> 200/255). Leaf is the DARK region (gray <= threshold).
GRAY_THRESHOLD <- 200L / 255

# Trim long petioles before computing length/width-type features?
TRIM_PETIOLE <- TRUE

# Reproducibility for the sanity-check classifier in script 06.
SEED <- 160

message("config.R loaded. DATA_DIR = ", DATA_DIR,
        " | species = ", paste(SPECIES_IDS, collapse = ", "))

