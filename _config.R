# Defines data paths from environment variables. See README section 4.
#   RICE_ROOT         output tree (02_MART, 05_REPORT); created if absent
#   RICE_DATA         input tree (01_ORG/RWIS_1MIN by default)
#   RICE_META         station metadata folder (default: RICE_DATA/00_META if present, else meta/ in this repository)
#   RICE_RWIS         folder of the 1-min RWIS files, if it is elsewhere
#   RICE_VALID_FILE   1-min records of the independent network (.rds)
#   RICE_VALID_META   station metadata of the independent network (.csv)
#   RICE_GRID_FIELDS  folder of the 500 m objectively analysed fields (R/spatial only)
#   RICE_GRID_LATLON  grid coordinate table (.csv)
# RICE_VALID_ROOT and RICE_GRID_ROOT may be set instead as the parent folders of those four.
# Sourced by the R scripts.

ROOT <- Sys.getenv("RICE_ROOT", unset = getwd())
DATA_ROOT <- Sys.getenv("RICE_DATA", unset = ROOT)
GRID_ROOT <- Sys.getenv("RICE_GRID_ROOT", unset = NA_character_)
VALID_ROOT <- Sys.getenv("RICE_VALID_ROOT", unset = NA_character_)
dir.create(ROOT, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(ROOT)) stop("RICE_ROOT could not be created: ", ROOT)
if (!dir.exists(DATA_ROOT)) stop("RICE_DATA does not exist: ", DATA_ROOT)
for (d in c("02_MART", "05_REPORT/FIGURES", "05_REPORT/LOG"))
  dir.create(file.path(ROOT, d), showWarnings = FALSE, recursive = TRUE)
META_DIR <- Sys.getenv("RICE_META", unset = if (dir.exists(file.path(DATA_ROOT, "00_META")))
  file.path(DATA_ROOT, "00_META") else file.path(CODE, "meta"))
meta_file <- function(name) {
  p <- file.path(META_DIR, name)
  if (!file.exists(p)) p <- file.path(CODE, "meta", name)
  if (!file.exists(p)) stop("station metadata file not found: ", name, "
  Expected in ", META_DIR, " or in meta/ of the repository.")
  p
}
RWIS_DIR <- Sys.getenv("RICE_RWIS", unset = file.path(DATA_ROOT, "01_ORG", "RWIS_1MIN"))
under <- function(root, name) if (is.na(root)) NA_character_ else file.path(root, name)
VALID_FILE <- Sys.getenv("RICE_VALID_FILE", unset = under(VALID_ROOT, "independent_1min.rds"))
VALID_META <- Sys.getenv("RICE_VALID_META", unset = under(VALID_ROOT, "independent_station_metadata.csv"))
GRID_FIELDS <- Sys.getenv("RICE_GRID_FIELDS", unset = under(GRID_ROOT, "fields"))
GRID_LATLON <- Sys.getenv("RICE_GRID_LATLON", unset = under(GRID_ROOT, "grid_latlon.csv"))
need <- function(path, var, what) {
  if (is.na(path))
    stop(what, ": set ", var, " to its location.")
  if (!file.exists(path) && !dir.exists(path))
    stop(what, " not found: ", path, "
  Set ", var, " to its location.")
  path
}
valid_file <- function() need(VALID_FILE, "RICE_VALID_FILE", "1-min records of the independent network")
valid_meta <- function() need(VALID_META, "RICE_VALID_META", "station metadata of the independent network")
grid_fields <- function() need(GRID_FIELDS, "RICE_GRID_FIELDS", "objectively analysed field archive")
grid_latlon <- function() need(GRID_LATLON, "RICE_GRID_LATLON", "grid coordinate table")
rwis_dir <- function() {
  if (!dir.exists(RWIS_DIR))
    stop("1-min RWIS folder not found: ", RWIS_DIR,
         "\n  Set RICE_RWIS to the folder that holds the 1-min files named *rwis_<STN_ID>_*.csv.")
  RWIS_DIR
}
code_file <- function(name) {
  hit <- list.files(CODE, pattern = paste0("^", gsub(".", "[.]", name, fixed = TRUE), "$"),
                    recursive = TRUE, full.names = TRUE)
  if (length(hit) != 1) stop("script not found in the repository: ", name)
  hit
}
