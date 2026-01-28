##### SST OSTIA DOWNLOAD AND PROCESS
library(reticulate)
library(tidyverse)
library(tidync)

## ----- user directories -----
root_dir <- "path/to/project"  # <-- change this for your setup
nc_dir   <- file.path(root_dir, "data", "ostia_nc")
rds_dir  <- file.path(root_dir, "data", "ostia_rds")

dir.create(nc_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

## ----- reticulate / copernicusmarine setup -----
if (!"CopernicusMarine" %in% virtualenv_list()) {
  virtualenv_create(envname = "CopernicusMarine")
  virtualenv_install(
    envname  = "CopernicusMarine",
    packages = c("copernicusmarine")
  )
}

reticulate::use_virtualenv("CopernicusMarine", required = TRUE)
cmt <- reticulate::import("copernicusmarine")

## ----- Adélie Sill–Commonwealth Bay domain -----
lon_min <- 140 # 140°E
lon_max <- 148 # 148°E
lat_min <- -68.5  # 68.5°S
lat_max <- -65.5  # 65.5°S

## ----- target NetCDF paths -----
rep1_nc <- file.path(nc_dir, "sst_REP1_001.nc")
rep2_nc <- file.path(nc_dir, "sst_REP2_001.nc")
nrt_nc  <- file.path(nc_dir, "sst_NRT_001.nc")

## ----- REP1 (1982-01-01 to 2001-12-31) -----
if (!file.exists(rep1_nc)) {
  before_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  
  cmt$subset(
    dataset_id        = "METOFFICE-GLO-SST-L4-REP-OBS-SST",
    variables         = list("analysed_sst"),
    minimum_longitude = lon_min,
    maximum_longitude = lon_max,
    minimum_latitude  = lat_min,
    maximum_latitude  = lat_max,
    start_datetime    = "1982-01-01T00:00:00",
    end_datetime      = "2001-12-31T00:00:00",
    output_directory  = nc_dir
  )
  
  after_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  new_file <- setdiff(after_files, before_files)
  if (length(new_file) == 1) {
    file.rename(new_file, rep1_nc)
  }
}

## ----- REP2 (2002-01-01 to 2022-05-31) -----
if (!file.exists(rep2_nc)) {
  before_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  
  cmt$subset(
    dataset_id        = "METOFFICE-GLO-SST-L4-REP-OBS-SST",
    variables         = list("analysed_sst"),
    minimum_longitude = lon_min,
    maximum_longitude = lon_max,
    minimum_latitude  = lat_min,
    maximum_latitude  = lat_max,
    start_datetime    = "2002-01-01T00:00:00",
    end_datetime      = "2022-05-31T00:00:00",
    output_directory  = nc_dir
  )
  
  after_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  new_file <- setdiff(after_files, before_files)
  if (length(new_file) == 1) {
    file.rename(new_file, rep2_nc)
  }
}

## ----- NRT (2022-06-01 to 2024-12-31) -----
if (!file.exists(nrt_nc)) {
  before_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  
  cmt$subset(
    dataset_id        = "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2",
    variables         = list("analysed_sst"),
    minimum_longitude = lon_min,
    maximum_longitude = lon_max,
    minimum_latitude  = lat_min,
    maximum_latitude  = lat_max,
    start_datetime    = "2022-06-01T00:00:00",
    end_datetime      = "2024-12-31T00:00:00",
    output_directory  = nc_dir
  )
  
  after_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  new_file <- setdiff(after_files, before_files)
  if (length(new_file) == 1) {
    file.rename(new_file, nrt_nc)
  }
}

gc()

## ----- read and standardize REP1 -----
sst_rep1 <- tidync::tidync(rep1_nc) |>
  tidync::hyper_tibble(force = TRUE) |>
  dplyr::mutate(
    time         = as.Date(time, format = "%Y-%m-%d"),
    analysed_sst = analysed_sst - 273.15   # K -> °C
  ) |>
  dplyr::arrange(time) |>
  dplyr::rename(
    t    = time,
    temp = analysed_sst,
    lon  = longitude,
    lat  = latitude
  )

gc()

## ----- read and standardize REP2 -----
sst_rep2 <- tidync::tidync(rep2_nc) |>
  tidync::hyper_tibble(force = TRUE) |>
  dplyr::mutate(
    time         = as.Date(time, format = "%Y-%m-%d"),
    analysed_sst = analysed_sst - 273.15
  ) |>
  dplyr::arrange(time) |>
  dplyr::rename(
    t    = time,
    temp = analysed_sst,
    lon  = longitude,
    lat  = latitude
  )

gc()

## ----- read and standardize NRT -----
sst_nrt <- tidync::tidync(nrt_nc) |>
  tidync::hyper_tibble(force = TRUE) |>
  dplyr::mutate(
    time         = as.Date(time, format = "%Y-%m-%d"),
    analysed_sst = analysed_sst - 273.15
  ) |>
  dplyr::arrange(time) |>
  dplyr::rename(
    t    = time,
    temp = analysed_sst,
    lon  = longitude,
    lat  = latitude
  )

gc()

## ----- combine and save -----
sst_data <- dplyr::bind_rows(sst_rep1, sst_rep2, sst_nrt) |>
  dplyr::arrange(t, lon, lat)

rm(sst_rep1, sst_rep2, sst_nrt)
gc()

saveRDS(
  sst_data,
  file = file.path(rds_dir, "sst_data_ostia_adelie.Rds")
)

##### MCS DETECTION
library(tidyverse)   # data manipulation
library(heatwaveR)   # MHW/MCS detection
library(tidync)      # NetCDF handling (if used upstream)
library(doParallel)  # parallel processing
library(plyr)        # ddply

# sst_data is the combined OSTIA SST object from the previous step
# with columns: t (Date), temp (°C), lon, lat

# cold-spell detection at a single grid cell
MCS_only <- function(df) {
  # climatology and threshold using 1982–2011 baseline
  clim <- ts2clm(
    data              = df,
    climatologyPeriod = c("1982-01-01", "2011-12-31"),
    pctile            = 5,              # 5th percentile threshold
  )
  
  # detect discrete cold-spell events
  event <- detect_event(
    data       = clim,
    coldSpells = TRUE
    # all other arguments left at heatwaveR defaults
  )
  
  # return only the event-level metrics
  event$event
}

# register cores (adjust as needed in practice)
registerDoParallel(cores = 7)

# apply cold-spell detection to each lon–lat pixel
mcs <- plyr::ddply(
  .data      = sst_data,
  .variables = c("lon", "lat"),
  .fun       = MCS_only,
  .parallel  = TRUE
)

# save event-level metrics for later analyses
saveRDS(
  mcs,
  file = file.path(rds_dir, "MCS_events_ostia_adelie.Rds")
)

gc()

##### MCS CUMULATIVE INTENSITY TRENDS
library(tidyverse)
library(lubridate)
library(plyr)
library(doParallel)

## annual cumulative intensity per grid cell and year
mcs_yi_hs <- mcs %>%
  as.data.frame() %>%
  dplyr::select(lon, lat, intensity_cumulative, date_start) %>%
  dplyr::mutate(year = lubridate::year(date_start)) %>%
  dplyr::group_by(lon, lat, year) %>%
  dplyr::summarise(
    intensity_cumulative = sum(intensity_cumulative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::filter(dplyr::between(year, 1982, 2024)) %>%
  tidyr::complete(year = 1982:2024) %>%
  dplyr::mutate(
    intensity_cumulative = ifelse(
      is.na(intensity_cumulative), 0, intensity_cumulative
    )
  ) %>%
  dplyr::ungroup()

## linear trend function for each grid cell
lin_fun <- function(ev) {
  mod1 <- glm(
    intensity_cumulative ~ year,
    family = gaussian(link = "identity"),
    data   = ev
  )
  data.frame(slope = summary(mod1)$coefficients[2, 1])
}

## compute trends per (lon, lat)
registerDoParallel(cores = 7)

mcs_yi_tr <- plyr::ddply(
  .data      = mcs_yi_hs,
  .variables = c("lon", "lat"),
  .fun       = lin_fun,
  .parallel  = TRUE
)

##### MCS ICE FLAGGING
library(tidyverse)
library(heatwaveR)
library(doParallel)

## ----- MCS climatology per pixel (1982–2011 baseline, 5th percentile) -----
MCS_daily <- function(df) {
  # 1) build cold-spell climatology using 5th-percentile threshold
  clim <- ts2clm(
    data              = df,
    climatologyPeriod = c("1982-01-01", "2011-12-31"),
    pctile            = 5,
    coldSpells        = TRUE
  )
  
  # 2) detect discrete cold-spell events
  event <- detect_event(
    data       = clim,
    coldSpells = TRUE
  )
  
  # 3) return the climatology (one row per day with threshold info)
  event$climatology
}

# register cores (adjust as needed)
registerDoParallel(cores = 7)

# apply climatology + event detection to each lon–lat pixel
mcs_clim <- plyr::ddply(
  .data      = sst_data,
  .variables = c("lon", "lat"),
  .fun       = MCS_daily,
  .parallel  = TRUE
)

## ----- “ice-MCS” flagging and summary statistics -----
mcs_ice_summary <- mcs_clim %>%
  # keep only days that are part of a cold-spell event
  dplyr::filter(event) %>%
  
  # flag ice-MCS days when the threshold is below -1.7 °C
  dplyr::mutate(
    icemcs = thresh <= -1.7
  ) %>%
  
  # summarize per lon–lat location
  dplyr::group_by(lon, lat) %>%
  dplyr::summarise(
    totaleventdays = dplyr::n(),                     # all cold-spell days
    totalice       = sum(icemcs, na.rm = TRUE),      # ice-MCS days
    totalnonice    = totaleventdays - totalice,      # non-ice MCS days
    perc_ice       = dplyr::if_else(
      totaleventdays > 0,
      100 * totalice / totaleventdays,
      NA_real_
    ),
    .groups = "drop"
  )

##### ANTLEADS LEAD FREQUENCY
library(terra)
library(dplyr)
library(stringr)

## ----- user directories -----
root_dir <- "path/to/project"  # <-- change this for your setup
ant_dir  <- file.path(root_dir, "data", "AntLeads")  # directory with AntLeads files

latlon_file <- file.path(ant_dir, "AntNSIDC_latlonmap_1km.nc")

## ----- read latitude/longitude map -----
ll <- terra::rast(latlon_file)  # DIMVAR warnings may appear but are expected

## 1) convert to data frame and rename lon/lat
ll_df <- as.data.frame(ll, cells = TRUE, na.rm = TRUE) |>
  dplyr::rename(
    lat = latmap,   # -90 .. -39
    lon = lonmap    # -180 .. 180
  )

## optional quick check
summary(ll_df[, c("lon", "lat")])

## 2) subset to Adélie Sill–Commonwealth Bay region
roi_df <- ll_df |>
  dplyr::filter(
    lon >= 140, lon <= 148,
    lat >= -68.5, lat <= -65.5
  )

roi_cells <- roi_df$cell

## 3) list monthly lead-frequency files (e.g., "200304_LF.nc")
lf_files <- list.files(
  ant_dir,
  pattern = "_LF\\.nc$",
  full.names = TRUE
)

## helper: read one LF file, subset to ROI, compute spatial mean
extract_lf_roi <- function(f) {
  r <- terra::rast(f)  # first layer expected to be meanLF
  
  # extract ROI values and undo ×1000 scaling
  vals <- terra::values(r)[roi_cells] / 1000
  
  lf_mean <- mean(vals, na.rm = TRUE)
  
  # parse year/month from filename "YYYYMM_LF.nc"
  fname <- basename(f)
  ym    <- stringr::str_sub(fname, 1, 6)
  yr    <- as.integer(stringr::str_sub(ym, 1, 4))
  mo    <- as.integer(stringr::str_sub(ym, 5, 6))
  
  tibble::tibble(
    file     = fname,
    year     = yr,
    month    = mo,
    year_dec = yr + (mo - 0.5) / 12,  # mid-month decimal year
    lf_mean  = lf_mean                # mean winter lead frequency in ROI
  )
}

## 4) build monthly lead-frequency time series for the box
antleads_ts <- lf_files |>
  lapply(extract_lf_roi) |>
  dplyr::bind_rows() |>
  dplyr::arrange(year, month)

str(antleads_ts)

##### B-SOSE SEA ICE PROXIES DATA
library(ncdf4)
library(dplyr)
library(lubridate)
library(stringr)
library(tibble)

## ----- user directories and files -----
root_dir <- "path/to/project"  # <-- change to your local project root
sose_dir <- file.path(root_dir, "data", "bsose")

grid_file    <- file.path(sose_dir, "grid.nc")
area_m_0812  <- file.path(sose_dir, "2008to2012_SeaIceArea_m.nc")
heff_m_0812  <- file.path(sose_dir, "2008to2012_SeaIceHeff_m.nc")
empmr_d_0812 <- file.path(sose_dir, "2008to2012_SeaIceEmPmR_d.nc")

fwflx_m_1324 <- file.path(sose_dir, "2013to2024_oceFWflx_m.nc")
area_m_1324  <- file.path(sose_dir, "2013to2024_SeaIceArea_m.nc")
empmr_m_1324 <- file.path(sose_dir, "2013to2024_SeaIceEmPmR_m.nc")
heff_m_1324  <- file.path(sose_dir, "2013to2024_SeaIceHeff_m.nc")

## ----- Adélie Sill–Commonwealth Bay box -----
roi_lonmin <- 140
roi_lonmax <- 148
roi_latmin <- -68.5
roi_latmax <- -65.5

rho_fw      <- 1000   # kg/m^3
sec_per_day <- 86400  # s/day

## ----- ROI mask from B-SOSE grid (XC/YC) -----
ncg <- nc_open(grid_file); on.exit(nc_close(ncg), add = TRUE)

XC <- ncvar_get(ncg, "XC")
YC <- ncvar_get(ncg, "YC")

## wrap longitudes from 0–360 to -180–180 if needed
XC <- ifelse(XC > 180, XC - 360, XC)

roi_mask <- XC >= roi_lonmin & XC <= roi_lonmax &
  YC >= roi_latmin & YC <= roi_latmax

message("ROI mask cells = ", sum(roi_mask, na.rm = TRUE))
if (sum(roi_mask, na.rm = TRUE) == 0) {
  stop("ROI mask empty — check bounds.")
}

## ----- CF time helper -----
cf_time_to_date <- function(time_vals, units_str) {
  units_str <- tolower(units_str)
  origin <- str_match(units_str, "since\\s*([0-9]{4}-[0-9]{2}-[0-9]{2})")[, 2]
  if (is.na(origin)) stop("Couldn't parse CF origin from: ", units_str)
  
  if (str_detect(units_str, "day")) {
    as.Date(time_vals, origin = origin)
  } else if (str_detect(units_str, "hour")) {
    as.POSIXct(
      time_vals * 3600,
      origin = paste0(origin, " 00:00:00"),
      tz = "UTC"
    ) |> as.Date()
  } else if (str_detect(units_str, "second")) {
    as.POSIXct(
      time_vals,
      origin = paste0(origin, " 00:00:00"),
      tz = "UTC"
    ) |> as.Date()
  } else {
    stop("Unsupported CF time units: ", units_str)
  }
}

## ----- robust B-SOSE extractor (monthly or daily) -----
extract_sose_ts <- function(nc_file, var, roi_mask, scale_div = 1,
                            freq = c("monthly", "daily")) {
  
  freq <- match.arg(freq)
  
  nc <- nc_open(nc_file); on.exit(nc_close(nc), add = TRUE)
  
  vnames <- names(nc$var)
  
  ## case-insensitive variable match
  var_idx <- which(tolower(vnames) == tolower(var))
  if (length(var_idx) == 0) {
    var_idx <- grep(tolower(var), tolower(vnames))
  }
  if (length(var_idx) == 0) {
    stop(
      "Variable ", var, " not in file. Vars are: ",
      paste(vnames, collapse = ", ")
    )
  }
  
  var_real <- vnames[var_idx[1]]
  message("Using variable: ", var_real)
  
  arr <- ncvar_get(nc, var_real) / scale_div
  d   <- dim(arr)
  dn  <- names(nc$var[[var_real]]$dim)
  
  ## identify time dimension
  tpos <- which(tolower(dn) %in% c("time", "t", "iter", "time_counter"))
  if (length(tpos) == 0) tpos <- length(d)
  
  ## move time to last dimension
  if (length(d) == 3) {
    perm <- c(setdiff(1:3, tpos), tpos)
    arr  <- aperm(arr, perm)
  } else if (length(d) == 2) {
    arr <- array(arr, dim = c(d, 1L))
  } else {
    stop("Unexpected dims for ", var_real, ": ", paste(d, collapse = "x"))
  }
  
  nt <- dim(arr)[3]
  
  ## ----- time handling -----
  dates <- rep(NA, nt)
  
  tname <- intersect(
    names(nc$dim),
    c("time", "TIME", "t", "iter", "time_counter")
  )
  if (length(tname) > 0) {
    tname     <- tname[1]
    time_vals <- nc$dim[[tname]]$vals
    units_str <- nc$dim[[tname]]$units
    tmp_dates <- cf_time_to_date(time_vals, units_str)
    
    if (length(tmp_dates) == nt) {
      dates <- tmp_dates
    } else {
      dates <- seq.Date(
        tmp_dates[1],
        by = ifelse(freq == "daily", "1 day", "1 month"),
        length.out = nt
      )
    }
  }
  
  if (all(is.na(dates))) {
    yrs_in_name <- str_match(basename(nc_file), "(\\d{4})to(\\d{4})")[, 2:3]
    if (all(!is.na(yrs_in_name))) {
      y0         <- as.integer(yrs_in_name[1])
      start_date <- as.Date(sprintf("%04d-01-01", y0))
      step       <- ifelse(freq == "daily", "1 day", "1 month")
      dates      <- seq.Date(start_date, by = step, length.out = nt)
    } else {
      stop("No CF time and can't infer from name: ", basename(nc_file))
    }
  }
  
  roi_means <- vapply(
    X   = 1:nt,
    FUN = function(k) {
      slice <- arr[, , k]
      mean(slice[roi_mask], na.rm = TRUE)
    },
    FUN.VALUE = numeric(1)
  )
  
  tibble(
    date     = dates,
    year     = year(dates),
    month    = month(dates),
    year_dec = year + (month - 0.5) / 12,
    value    = roi_means,
    var      = var_real
  ) |>
    arrange(date)
}

## ----- extract time series from B-SOSE -----
## 2008–2012 monthly
ts_area_0812 <- extract_sose_ts(area_m_0812, "SIArea", roi_mask, freq = "monthly")
ts_heff_0812 <- extract_sose_ts(heff_m_0812, "SIHeff", roi_mask, freq = "monthly")

## 2008–2012 daily EmPmR -> monthly means
ts_empmr_d_0812 <- extract_sose_ts(empmr_d_0812, "SIEmPmR", roi_mask, freq = "daily")

ts_empmr_0812 <- ts_empmr_d_0812 |>
  group_by(year, month) |>
  summarise(
    value    = mean(value, na.rm = TRUE),
    year_dec = first(year_dec),
    date     = first(date),
    .groups  = "drop"
  ) |>
  mutate(var = "SIEmPmR")

## 2013–2024 monthly
ts_area_1324  <- extract_sose_ts(area_m_1324,  "SIArea",   roi_mask, freq = "monthly")
ts_heff_1324  <- extract_sose_ts(heff_m_1324,  "SIHeff",   roi_mask, freq = "monthly")
ts_empmr_1324 <- extract_sose_ts(empmr_m_1324, "SIEmPmR",  roi_mask, freq = "monthly")
ts_fwflx_1324 <- extract_sose_ts(fwflx_m_1324, "oceFWflx", roi_mask, freq = "monthly")

## ----- combine periods -----
ts_area  <- bind_rows(ts_area_0812,  ts_area_1324)  |> arrange(date)
ts_heff  <- bind_rows(ts_heff_0812,  ts_heff_1324)  |> arrange(date)
ts_empmr <- bind_rows(ts_empmr_0812, ts_empmr_1324) |> arrange(date)
ts_fwflx <- ts_fwflx_1324 |> arrange(date)

## ----- derived proxies -----
## open-water fraction
ts_openwater <- ts_area |>
  mutate(openwater = 1 - value) |>
  select(date, year, month, year_dec, openwater)

## sea-ice volume proxy (area * thickness)
ts_sivol <- ts_area |>
  rename(SIArea = value) |>
  inner_join(
    ts_heff |> rename(SIHeff = value),
    by = c("date", "year", "month", "year_dec")
  ) |>
  mutate(sivol = SIArea * SIHeff) |>
  select(date, year, month, year_dec, sivol)

## freshwater flux from sea ice (m/day)
ts_empmr_mday <- ts_empmr |>
  mutate(fw_m_per_day = (value / rho_fw) * sec_per_day) |>
  select(date, year, month, year_dec, fw_m_per_day)

## ----- winter summaries (Apr–Sep) -----
winter_months <- 4:9

openwater_winter <- ts_openwater |>
  filter(month %in% winter_months) |>
  group_by(year) |>
  summarise(
    openwater_winter = mean(openwater, na.rm = TRUE),
    .groups = "drop"
  )

sivol_winter <- ts_sivol |>
  filter(month %in% winter_months) |>
  group_by(year) |>
  summarise(
    sivol_winter = mean(sivol, na.rm = TRUE),
    .groups = "drop"
  )

fwseaice_winter <- ts_empmr_mday |>
  filter(month %in% winter_months) |>
  group_by(year) |>
  summarise(
    fwseaice_winter_m = sum(fw_m_per_day, na.rm = TRUE),
    .groups = "drop"
  )

fwtotal_winter <- ts_fwflx |>
  mutate(fw_m_per_day = (value / rho_fw) * sec_per_day) |>
  filter(month %in% winter_months) |>
  group_by(year) |>
  summarise(
    fwtotal_winter_m = sum(fw_m_per_day, na.rm = TRUE),
    .groups = "drop"
  )

## ----- final data frames used in plots -----
fig3c_df <- sivol_winter |>
  left_join(openwater_winter, by = "year")

fig3d_df <- fwseaice_winter |>
  left_join(fwtotal_winter, by = "year")

str(fig3c_df)
str(fig3d_df)

##### mCDW DIAGNOSTICS
library(ncdf4)
library(data.table)
library(dplyr)
library(lubridate)
library(gsw)
library(ggplot2)

## ----- user paths -----
root <- "path/to/data"  # <-- change to your local data directory

files <- list(
  ## 2008–2012 (iteration 106)
  th_106 = file.path(root, "bsose_i106_2008to2012_monthly_Theta.nc"),
  sl_106 = file.path(root, "bsose_i106_2008to2012_monthly_Salt.nc"),
  o2_106 = file.path(root, "bsose_i106_2008to2012_monthly_O2.nc"),
  u_106  = file.path(root, "bsose_i106_2008to2012_monthly_Uvel.nc"),
  v_106  = file.path(root, "bsose_i106_2008to2012_monthly_Vvel.nc"),
  
  ## 2013–2023 (iteration 155)
  th_155 = file.path(root, "Theta_bsoseI155_2013to2023_monthly.nc"),
  sl_155 = file.path(root, "Salt_bsoseI155_2013to2023_monthly.nc"),
  o2_155 = file.path(root, "O2_bsoseI155_2013to2023_monthly.nc"),
  u_155  = file.path(root, "Uvel_bsoseI155_2013to2023_monthly.nc"),
  v_155  = file.path(root, "Vvel_bsoseI155_2013to2023_monthly.nc"),
  
  grid   = file.path(root, "grid.nc")
)

## ----- Adélie Sill–Commonwealth Bay box -----
lonmin <- 140
lonmax <- 148
latmin <- -68.5
latmax <- -65.5

## ----- depth settings -----
zmin_m    <- 200   # VmCDW only below this depth
ww_zmax   <- 400   # WW endmember search layer top
mcdw_zmin <- 500   # mCDW endmember deeper than this

## ----- helper functions -----
guess_varname <- function(nc, candidates) {
  vars       <- names(nc$var)
  lower_vars <- tolower(vars)
  lower_cand <- tolower(candidates)
  
  hit <- candidates[lower_cand %in% lower_vars][1]
  if (!is.na(hit)) {
    return(vars[lower_vars == tolower(hit)][1])
  }
  
  stop(
    "Could not find variable among candidates: ",
    paste(candidates, collapse = ", "),
    "\nVars are: ", paste(vars, collapse = ", ")
  )
}

parse_time <- function(nc) {
  possible <- c("time", "TIME", "t", "T", "iter")
  tname    <- possible[possible %in% c(names(nc$dim), names(nc$var))][1]
  if (is.na(tname)) stop("No time dimension/variable found.")
  
  if (tname %in% names(nc$dim)) {
    tvals <- nc$dim[[tname]]$vals
  } else {
    tvals <- ncvar_get(nc, tname)
  }
  
  units <- tryCatch(
    ncatt_get(nc, tname, "units")$value,
    error = function(e) NA_character_
  )
  
  if (!is.na(units) && grepl("since", units)) {
    origin_str <- sub(".*since\\s+", "", units)
    origin     <- as.POSIXct(origin_str, tz = "UTC")
    
    if (grepl("days", units)) {
      dates <- origin + days(tvals)
    } else if (grepl("hours", units)) {
      dates <- origin + hours(tvals)
    } else if (grepl("seconds", units)) {
      dates <- origin + seconds(tvals)
    } else {
      dates <- origin + days(tvals)
    }
    return(as.Date(dates))
  }
  
  fname <- basename(nc$filename)
  yr    <- suppressWarnings(
    as.integer(regmatches(fname, regexpr("[0-9]{4}", fname)))
  )
  if (is.na(yr)) yr <- 2008
  
  seq.Date(
    as.Date(paste0(yr, "-01-15")),
    by         = "1 month",
    length.out = length(tvals)
  )
}

read_slice_rect <- function(nc, vname, i_rng, j_rng, k_rng, t_idx) {
  start <- c(i_rng[1], j_rng[1], k_rng[1], t_idx)
  count <- c(diff(i_rng) + 1, diff(j_rng) + 1, diff(k_rng) + 1, 1)
  ncvar_get(nc, vname, start = start, count = count)
}

## ----- read grid and build mask/indices -----
ncg <- nc_open(files$grid)

XC    <- ncvar_get(ncg, "XC")
YC    <- ncvar_get(ncg, "YC")
RC    <- ncvar_get(ncg, "RC")     # depth centers (negative)
drF   <- ncvar_get(ncg, "DRF")    # layer thickness
RAC   <- ncvar_get(ncg, "RAC")    # cell area at C points
hFacC <- ncvar_get(ncg, "hFacC")
hFacS <- ncvar_get(ncg, "hFacS")  # V faces
hFacW <- ncvar_get(ncg, "hFacW")  # U faces
DXG   <- ncvar_get(ncg, "DXG")    # zonal length at V faces
DYG   <- ncvar_get(ncg, "DYG")    # meridional length at U faces

nc_close(ncg)

mask2d <- (XC >= lonmin & XC <= lonmax & YC >= latmin & YC <= latmax)
ij     <- which(mask2d, arr.ind = TRUE)
i_rng  <- range(ij[, 1])
j_rng  <- range(ij[, 2])

k_idx_all <- which(abs(RC) >= zmin_m)
k_rng     <- range(k_idx_all)
k_idx_ww  <- which(abs(RC) >= zmin_m & abs(RC) <= ww_zmax)

## boundary indices
j_north <- max(ij[, 2])
i_west  <- min(ij[, 1])
i_east  <- max(ij[, 1])

## local boundary indices inside rectangle
jloc_north <- j_north - j_rng[1] + 1
iloc_west  <- i_west  - i_rng[1] + 1
iloc_east  <- i_east  - i_rng[1] + 1

## subset rectangles
XC_rect    <- XC[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2]]
YC_rect    <- YC[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2]]
mask_rect  <- mask2d[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2]]

RAC_rect   <- RAC[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2]]
hFacC_rect <- hFacC[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2], k_rng[1]:k_rng[2]]
hFacS_rect <- hFacS[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2], k_rng[1]:k_rng[2]]
hFacW_rect <- hFacW[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2], k_rng[1]:k_rng[2]]
DXG_rect   <- DXG[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2]]
DYG_rect   <- DYG[i_rng[1]:i_rng[2], j_rng[1]:j_rng[2]]

drF_rect <- drF[k_rng[1]:k_rng[2]]

## precompute cell volumes (m^3)
cellVol_rect <- array(NA_real_, dim = dim(hFacC_rect))
for (kk in seq_along(drF_rect)) {
  cellVol_rect[,, kk] <- RAC_rect * drF_rect[kk] * hFacC_rect[,, kk]
}

## ----- process one B-SOSE period -----
process_period <- function(th_path, sl_path, o2_path,
                           u_path = NULL, v_path = NULL,
                           period_tag = "") {
  
  cat("\n--- Processing period:", period_tag, "---\n")
  
  nc_th <- nc_open(th_path)
  nc_sl <- nc_open(sl_path)
  nc_o2 <- nc_open(o2_path)
  nc_u  <- if (!is.null(u_path)) nc_open(u_path) else NULL
  nc_v  <- if (!is.null(v_path)) nc_open(v_path) else NULL
  
  on.exit({
    nc_close(nc_th); nc_close(nc_sl); nc_close(nc_o2)
    if (!is.null(nc_u)) nc_close(nc_u)
    if (!is.null(nc_v)) nc_close(nc_v)
  }, add = TRUE)
  
  vth <- guess_varname(nc_th, c("Theta", "THETA", "theta", "temp"))
  vsl <- guess_varname(nc_sl, c("Salt", "SALT", "salt", "salinity"))
  vo2 <- guess_varname(nc_o2, c("TRAC03", "trac03", "O2", "oxygen", "o2"))
  
  vu  <- if (!is.null(nc_u)) guess_varname(nc_u, c("Uvel", "UVEL", "uVel", "u")) else NA
  vv  <- if (!is.null(nc_v)) guess_varname(nc_v, c("Vvel", "VVEL", "vVel", "v")) else NA
  
  time_vec   <- parse_time(nc_th)
  years      <- year(time_vec)
  uniq_years <- sort(unique(years))
  
  endm <- data.table(
    year      = uniq_years,
    theta_m   = NA_real_,
    salt_m    = NA_real_,
    o2_m      = NA_real_,
    sigma0_m  = -Inf,
    theta_ww  = NA_real_,
    salt_ww   = NA_real_,
    o2_ww     = NA_real_,
    CT_ww     = Inf
  )
  
  ## ----- pass 1: derive endmembers -----
  cat("Pass 1: deriving endmembers per year...\n")
  
  for (tt in seq_along(time_vec)) {
    
    yr <- years[tt]
    
    th <- read_slice_rect(nc_th, vth, i_rng, j_rng, k_rng, tt)
    sl <- read_slice_rect(nc_sl, vsl, i_rng, j_rng, k_rng, tt)
    o2 <- read_slice_rect(nc_o2, vo2, i_rng, j_rng, k_rng, tt)
    
    th[!mask_rect] <- NA
    sl[!mask_rect] <- NA
    o2[!mask_rect] <- NA
    
    i_len <- dim(th)[1]
    j_len <- dim(th)[2]
    k_len <- dim(th)[3]
    
    lon_vec <- rep(as.vector(XC_rect), times = k_len)
    lat_vec <- rep(as.vector(YC_rect), times = k_len)
    z_vec   <- rep(RC[k_rng[1]:k_rng[2]], each = i_len * j_len)
    
    th_vec <- as.vector(th)
    sl_vec <- as.vector(sl)
    o2_vec <- as.vector(o2)
    
    good <- is.finite(th_vec) & is.finite(sl_vec) & is.finite(o2_vec)
    if (!any(good)) next
    
    lon_g <- lon_vec[good]
    lat_g <- lat_vec[good]
    z_g   <- z_vec[good]
    th_g  <- th_vec[good]
    sl_g  <- sl_vec[good]
    o2_g  <- o2_vec[good]
    
    p_g  <- gsw_p_from_z(z_g, lat_g)
    SA_g <- gsw_SA_from_SP(sl_g, p_g, lon_g, lat_g)
    CT_g <- gsw_CT_from_pt(SA_g, th_g)
    sig0 <- gsw_sigma0(SA_g, CT_g)
    
    irow <- which(endm$year == yr)
    
    ## mCDW: densest deep tail
    deep_layer <- abs(z_g) >= mcdw_zmin
    if (any(deep_layer)) {
      sig_deep <- sig0[deep_layer]
      q_hi     <- quantile(sig_deep, 0.995, na.rm = TRUE)
      hi_set   <- deep_layer & (sig0 >= q_hi)
      if (any(hi_set)) {
        endm$sigma0_m[irow] <- max(sig0[hi_set], na.rm = TRUE)
        endm$theta_m[irow]  <- median(th_g[hi_set], na.rm = TRUE)
        endm$salt_m[irow]   <- median(sl_g[hi_set], na.rm = TRUE)
        endm$o2_m[irow]     <- median(o2_g[hi_set], na.rm = TRUE)
      }
    }
    
    ## WW: cold tail 200–400 m
    ww_layer <- abs(z_g) >= zmin_m & abs(z_g) <= ww_zmax
    if (any(ww_layer)) {
      CTw    <- CT_g[ww_layer]
      q_lo   <- quantile(CTw, 0.05, na.rm = TRUE)
      lo_set <- ww_layer & (CT_g <= q_lo)
      if (any(lo_set)) {
        endm$CT_ww[irow]    <- min(CT_g[lo_set], na.rm = TRUE)
        endm$theta_ww[irow] <- median(th_g[lo_set], na.rm = TRUE)
        endm$salt_ww[irow]  <- median(sl_g[lo_set], na.rm = TRUE)
        endm$o2_ww[irow]    <- median(o2_g[lo_set], na.rm = TRUE)
      }
    }
    
    if (tt %% 12 == 0) cat("  ", period_tag, ":", yr, "done\n")
  }
  
  print(endm)
  
  ## ----- repair bad endmember years -----
  endm[, dtheta := theta_m - theta_ww]
  endm[, dsalt  := salt_m  - salt_ww]
  
  bad <- which(
    !is.finite(endm$dtheta) | !is.finite(endm$dsalt) |
      abs(endm$dtheta) < 0.05 | abs(endm$dsalt) < 0.02
  )
  
  if (length(bad) > 0) {
    cat("Repairing bad endmember years:", endm$year[bad], "\n")
    good_idx <- setdiff(seq_len(nrow(endm)), bad)
    
    for (ii in bad) {
      yg      <- endm$year[good_idx]
      nearest <- good_idx[which.min(abs(yg - endm$year[ii]))]
      
      endm$theta_m[ii]  <- endm$theta_m[nearest]
      endm$salt_m[ii]   <- endm$salt_m[nearest]
      endm$o2_m[ii]     <- endm$o2_m[nearest]
      endm$sigma0_m[ii] <- endm$sigma0_m[nearest]
      
      endm$theta_ww[ii] <- endm$theta_ww[nearest]
      endm$salt_ww[ii]  <- endm$salt_ww[nearest]
      endm$o2_ww[ii]    <- endm$o2_ww[nearest]
      endm$CT_ww[ii]    <- endm$CT_ww[nearest]
    }
    
    endm[, dtheta := theta_m - theta_ww]
    endm[, dsalt  := salt_m  - salt_ww]
  }
  
  ## ----- pass 2: monthly VmCDW, volume, transport -----
  cat("Pass 2: computing monthly VmCDW, volume, transport...\n")
  
  out_monthly <- data.table(
    time        = time_vec,
    year        = years,
    month       = month(time_vec),
    H_mCDW_m    = NA_real_,
    V_mCDW_m3   = NA_real_,
    S_mCDW_gkg  = NA_real_,
    Qn_in_m3s   = NA_real_,
    Qe_in_m3s   = NA_real_,
    Qw_in_m3s   = NA_real_,
    Qtot_in_m3s = NA_real_
  )
  
  for (tt in seq_along(time_vec)) {
    
    yr <- years[tt]
    em <- endm[year == yr]
    
    th <- read_slice_rect(nc_th, vth, i_rng, j_rng, k_rng, tt)
    sl <- read_slice_rect(nc_sl, vsl, i_rng, j_rng, k_rng, tt)
    
    th[!mask_rect] <- NA
    sl[!mask_rect] <- NA
    
    Vm_th <- (th - em$theta_ww) / (em$theta_m - em$theta_ww)
    Vm_sl <- (sl - em$salt_ww)  / (em$salt_m  - em$salt_ww)
    Vm    <- (Vm_th + Vm_sl) / 2
    Vm    <- pmin(pmax(Vm, 0), 1)
    
    Vm_vol <- Vm * cellVol_rect
    denom  <- sum(Vm_vol, na.rm = TRUE)
    
    S_mCDW <- if (is.finite(denom) && denom > 0) {
      sum(sl * Vm_vol, na.rm = TRUE) / denom
    } else {
      NA_real_
    }
    out_monthly$S_mCDW_gkg[tt] <- S_mCDW
    
    thick_cell <- array(0, dim = dim(Vm)[1:2])
    vol_cell   <- array(0, dim = dim(Vm)[1:2])
    
    for (kk in seq_len(dim(Vm)[3])) {
      thick_cell <- thick_cell + Vm[,, kk] * drF_rect[kk] * hFacC_rect[,, kk]
      vol_cell   <- vol_cell   + Vm[,, kk] * cellVol_rect[,, kk]
    }
    
    Hbar <- sum(thick_cell * RAC_rect, na.rm = TRUE) /
      sum(RAC_rect[mask_rect], na.rm = TRUE)
    Vtot <- sum(vol_cell, na.rm = TRUE)
    
    out_monthly$H_mCDW_m[tt]  <- Hbar
    out_monthly$V_mCDW_m3[tt] <- Vtot
    
    ## north edge (Vvel)
    Qn_in <- 0
    if (!is.null(nc_v)) {
      vsec <- ncvar_get(
        nc_v, vv,
        start = c(i_rng[1], j_north, k_rng[1], tt),
        count = c(diff(i_rng) + 1, 1, diff(k_rng) + 1, 1)
      )
      vsec <- drop(vsec)
      
      Vm_face <- Vm[, jloc_north, ]
      hS_face <- hFacS_rect[, jloc_north, ]
      dx_face <- DXG_rect[, jloc_north]
      
      area_face <- vsec * NA_real_
      for (kk in seq_len(ncol(vsec))) {
        area_face[, kk] <- dx_face * drF_rect[kk] * hS_face[, kk]
      }
      
      Qn_in <- sum(pmax(-vsec, 0) * Vm_face * area_face, na.rm = TRUE)
      out_monthly$Qn_in_m3s[tt] <- Qn_in
    }
    
    ## east/west edges (Uvel)
    Qe_in <- 0
    Qw_in <- 0
    if (!is.null(nc_u)) {
      ## east boundary: inflow = westward (u < 0)
      usec_e <- ncvar_get(
        nc_u, vu,
        start = c(i_east, j_rng[1], k_rng[1], tt),
        count = c(1, diff(j_rng) + 1, diff(k_rng) + 1, 1)
      )
      usec_e <- drop(usec_e)
      
      Vm_e  <- Vm[iloc_east, , ]
      hW_e  <- hFacW_rect[iloc_east, , ]
      dy_e  <- DYG_rect[iloc_east, ]
      
      area_e <- usec_e * NA_real_
      for (kk in seq_len(ncol(usec_e))) {
        area_e[, kk] <- dy_e * drF_rect[kk] * hW_e[, kk]
      }
      
      Qe_in <- sum(pmax(-usec_e, 0) * Vm_e * area_e, na.rm = TRUE)
      out_monthly$Qe_in_m3s[tt] <- Qe_in
      
      ## west boundary: inflow = eastward (u > 0)
      usec_w <- ncvar_get(
        nc_u, vu,
        start = c(i_west, j_rng[1], k_rng[1], tt),
        count = c(1, diff(j_rng) + 1, diff(k_rng) + 1, 1)
      )
      usec_w <- drop(usec_w)
      
      Vm_w  <- Vm[iloc_west, , ]
      hW_w  <- hFacW_rect[iloc_west, , ]
      dy_w  <- DYG_rect[iloc_west, ]
      
      area_w <- usec_w * NA_real_
      for (kk in seq_len(ncol(usec_w))) {
        area_w[, kk] <- dy_w * drF_rect[kk] * hW_w[, kk]
      }
      
      Qw_in <- sum(pmax(usec_w, 0) * Vm_w * area_w, na.rm = TRUE)
      out_monthly$Qw_in_m3s[tt] <- Qw_in
    }
    
    out_monthly$Qtot_in_m3s[tt] <- Qn_in + Qe_in + Qw_in
    
    if (tt %% 12 == 0) {
      cat("  ", period_tag, ":", yr, "month", month(time_vec[tt]), "done\n")
    }
  }
  
  list(
    endmembers = endm,
    monthly    = out_monthly
  )
}

## ----- run both periods and form annual means -----
res_106 <- process_period(
  files$th_106, files$sl_106, files$o2_106,
  u_path = files$u_106, v_path = files$v_106,
  period_tag = "i106_2008-2012"
)

res_155 <- process_period(
  files$th_155, files$sl_155, files$o2_155,
  u_path = files$u_155, v_path = files$v_155,
  period_tag = "i155_2013-2023"
)

monthly_all <- rbindlist(list(res_106$monthly, res_155$monthly)) |>
  arrange(time)

endm_all <- rbindlist(list(
  res_106$endmembers |> mutate(period = "i106"),
  res_155$endmembers |> mutate(period = "i155")
))

annual_all <- monthly_all |>
  group_by(year) |>
  summarise(
    H_mCDW_m   = mean(H_mCDW_m,   na.rm = TRUE),
    V_mCDW_m3  = mean(V_mCDW_m3,  na.rm = TRUE),
    S_mCDW_gkg = mean(S_mCDW_gkg, na.rm = TRUE),
    Qn_in_Sv   = mean(Qn_in_m3s,   na.rm = TRUE) / 1e6,
    Qe_in_Sv   = mean(Qe_in_m3s,   na.rm = TRUE) / 1e6,
    Qw_in_Sv   = mean(Qw_in_m3s,   na.rm = TRUE) / 1e6,
    Qtot_in_Sv = mean(Qtot_in_m3s, na.rm = TRUE) / 1e6,
    .groups    = "drop"
  )

saveRDS(monthly_all, file.path(root, "adelie_mCDW_monthly_withQ_3edges.rds"))
saveRDS(annual_all,  file.path(root, "adelie_mCDW_annual_withQ_3edges.rds"))
saveRDS(endm_all,    file.path(root, "adelie_endmembers_annual.rds"))

##### OTHER PARAMETERS
## ----- packages -----
library(terra)      # rasters / NetCDF / ASC / GRIB
library(sf)         # shapefiles
library(dplyr)
library(lubridate)
library(ggplot2)
library(pals)       # palettes
library(cowplot)    # for side-by-side panels

## ----- paths -----
data_dir <- "path/to/data"  # <-- change to your local data directory

bathy_file <- file.path(data_dir, "02_gvdem500v3.asc")
u_file     <- file.path(data_dir, "03_east_seavelocity.nc")
v_file     <- file.path(data_dir, "04_north_seavelocity.nc")
ice_file   <- file.path(data_dir, "05_icebergs", "Icebergs_20251114.shp")
b9b_csv    <- file.path(data_dir, "06_b09b.csv")

## ERA5 10 m wind (GRIB)
wu_file <- file.path(data_dir, "07_10m_wind_u.grib")
wv_file <- file.path(data_dir, "08_10m_wind_v.grib")

## ----- region of interest (Adélie–George V Land shelf) -----
roi <- ext(140, 150, -68.5, -64.5)   # [xmin, xmax, ymin, ymax]

## ----- bathymetry (gvdem) -----
bathy     <- rast(bathy_file)
bathy_roi <- crop(bathy, roi)

bathy_df <- as.data.frame(bathy_roi, xy = TRUE, na.rm = FALSE)
names(bathy_df)[3] <- "z"

## depth binning (0–1000 m in 100 m steps, then 2000, 4000)
depth_breaks <- c(0, seq(100, 1000, 100), 2000, 4000)

bathy_df <- bathy_df |>
  mutate(
    depth_pos = ifelse(z < 0, -z, NA_real_),
    depth_pos = pmin(pmax(depth_pos, 0), 3700),
    depth_band = cut(
      depth_pos,
      breaks = depth_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  )

depth_levels <- levels(bathy_df$depth_band)
bathy_cols   <- hcl.colors(length(depth_levels), "Purp", rev = TRUE)
names(bathy_cols) <- depth_levels

contour_breaks <- depth_breaks

## ----- ocean velocity (mean u, v) -----
u <- rast(u_file) |> crop(roi)
v <- rast(v_file) |> crop(roi)

u_mean <- app(u, mean, na.rm = TRUE)
v_mean <- app(v, mean, na.rm = TRUE)

vel_df <- as.data.frame(c(u_mean, v_mean), xy = TRUE, na.rm = TRUE)
names(vel_df)[3:4] <- c("u", "v")

## subsample for clarity
vel_df <- vel_df[seq(1, nrow(vel_df), by = 6), ]

vel_df2 <- vel_df |>
  mutate(
    speed = sqrt(u^2 + v^2),
    u_dir = if_else(speed == 0, 0, u / speed),
    v_dir = if_else(speed == 0, 0, v / speed)
  )

base_len <- 0.25

vel_df2 <- vel_df2 |>
  mutate(
    s_norm = ifelse(
      max(speed, na.rm = TRUE) > 0,
      speed / max(speed, na.rm = TRUE),
      0
    ),
    len  = (0.3 + 0.7 * s_norm) * base_len,
    xend = x + u_dir * len,
    yend = y + v_dir * len
  )

## ----- 10 m wind (ERA5) -----
wu <- rast(wu_file) |> crop(roi)
wv <- rast(wv_file) |> crop(roi)

wu_mean <- app(wu, mean, na.rm = TRUE)
wv_mean <- app(wv, mean, na.rm = TRUE)

wind_df <- as.data.frame(c(wu_mean, wv_mean), xy = TRUE, na.rm = TRUE)
names(wind_df)[3:4] <- c("u", "v")

## subsample
wind_df <- wind_df[seq(1, nrow(wind_df), by = 6), ]

wind_df2 <- wind_df |>
  mutate(
    speed = sqrt(u^2 + v^2),
    u_dir = if_else(speed == 0, 0, u / speed),
    v_dir = if_else(speed == 0, 0, v / speed)
  )

wind_df2 <- wind_df2 |>
  mutate(
    s_norm = ifelse(
      max(speed, na.rm = TRUE) > 0,
      speed / max(speed, na.rm = TRUE),
      0
    ),
    len  = (0.3 + 0.7 * s_norm) * base_len,
    xend = x + u_dir * len,
    yend = y + v_dir * len
  )

##### ROBUSTNESS
library(tidyverse)   # data manipulation
library(heatwaveR)   # MHW/MCS detection
library(tidync)      # NetCDF handling (if used upstream)
library(doParallel)  # parallel processing
library(plyr)        # ddply

# sst_data is the combined OSTIA SST object from the previous step
# with columns: t (Date), temp (°C), lon, lat

# cold-spell detection at a single grid cell
MCS_only2 <- function(df) {
  # climatology and threshold using 1982–2011 baseline
  clim <- ts2clm(
    data              = df,
    climatologyPeriod = c("1991-01-01", "2020-12-31"),
    pctile            = 5              # 5th percentile threshold
  )
  
  # detect discrete cold-spell events
  event <- detect_event(
    data       = clim,
    coldSpells = TRUE
    # all other arguments left at heatwaveR defaults
  )
  
  # return only the event-level metrics
  event$event
}

# register cores (adjust as needed in practice)
registerDoParallel(cores = 7)

# apply cold-spell detection to each lon–lat pixel
mcs2 <- plyr::ddply(
  .data      = sst_data,
  .variables = c("lon", "lat"),
  .fun       = MCS_only,
  .parallel  = TRUE
)

# save event-level metrics for later analyses
saveRDS(
  mcs2,
  file = file.path(rds_dir, "MCS_events_ostia_adelie2.Rds")
)

gc()

library(tidyverse)
library(lubridate)
library(plyr)
library(doParallel)

## annual cumulative intensity per grid cell and year
mcs_yi_hs2 <- mcs2 %>%
  as.data.frame() %>%
  dplyr::select(lon, lat, intensity_cumulative, date_start) %>%
  dplyr::mutate(year = lubridate::year(date_start)) %>%
  dplyr::group_by(lon, lat, year) %>%
  dplyr::summarise(
    intensity_cumulative = sum(intensity_cumulative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::filter(dplyr::between(year, 1982, 2024)) %>%
  tidyr::complete(year = 1982:2024) %>%
  dplyr::mutate(
    intensity_cumulative = ifelse(
      is.na(intensity_cumulative), 0, intensity_cumulative
    )
  ) %>%
  dplyr::ungroup()

## linear trend function for each grid cell
lin_fun <- function(ev) {
  mod1 <- glm(
    intensity_cumulative ~ year,
    family = gaussian(link = "identity"),
    data   = ev
  )
  data.frame(slope = summary(mod1)$coefficients[2, 1])
}

## compute trends per (lon, lat)
registerDoParallel(cores = 7)

mcs_yi_tr2 <- plyr::ddply(
  .data      = mcs_yi_hs2,
  .variables = c("lon", "lat"),
  .fun       = lin_fun,
  .parallel  = TRUE
)

mcs_df2 <- mcs_yi_tr2 %>%
  dplyr::mutate(lon = as.numeric(lon),
                lat = as.numeric(lat))











