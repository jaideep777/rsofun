params <-
  list(output_dir = ".")

#' 
## --------------------------------------------------------------------
library(tidyverse)
library(reshape2)
library(FluxDataKit)
library(lubridate)
library(ncdf4)
library(here)
source(here("vignettes","read_meta_fdk.R"))

# load stuff
# use for cloud coverage and le 
zenodo_driver = readRDS(here("data","filtered_driver.rds"))

valid_years = read.csv(here("ancillary_data","valid_years_final.csv"), header = T, fileEncoding = "UTF-16")

load(here("data","p_model_drivers.rda"))


nc = nc_open(here("ancillary_data","cwdx80.nc"))
lons = ncvar_get(nc, "lon")
lats = ncvar_get(nc, "lat")
S80 = ncvar_get(nc, "cwdx80")

args = commandArgs(trailingOnly=TRUE)


# insert your path here
lsm_path = paste0(here("data-raw","LSM"))
csv_path = paste0(here("data-raw","CSV"))

files_csv = list.files(csv_path)
files_lsm = list.files(lsm_path)

####-----
# filter sites
####-----
# Get filename for HH data for matching site
files_csv = files_csv[grep("HH", files_csv)] 

# remove sites with cropland, wetland and bad le_corr
sites = sapply(files_csv, function(x){substr(x,5,10)})

keep <- fdk_site_info|>
  filter(igbp_land_use != "CRO" & igbp_land_use != "WET")

sites = sites[which(sites %in% keep$sitename)]

keep <- fdk_site_fullyearsequence |>
  filter(drop_lecorr != TRUE)

sites = sites[which(sites %in% keep$sitename)]

zenodo_driver = zenodo_driver[which(zenodo_driver$sitename %in% sites),]

driver = NULL

####-----
# for cycle
####-----

#from here starts for cycle
num = 0
for (i in sites[1:20]){ #change to all
  num = num + 1
  message("site number ",num)
  
  meta <- suppressWarnings(
    try(
      read_meta_fdk(
        site = i,
        path = lsm_path,
        meta_data = T
      )
    )
  )
  
  
  # select years with good le_corr
  ystart = fdk_site_fullyearsequence[fdk_site_fullyearsequence$sitename == i,]$year_start_lecorr
  yend = fdk_site_fullyearsequence[fdk_site_fullyearsequence$sitename == i,]$year_end_lecorr

  file_csv = files_csv[grep(i,files_csv)]
  
  message("- reading FLUXNET format halfhourly data")
  hhdf <- readr::read_csv(paste0(csv_path,"/",file_csv))
  
  # Add date and time columns to hhdf for easier further processing.
  # ---------------------------------------------------------
  hhdf =
    hhdf |>
    mutate(time = lubridate::as_datetime(as.character(TIMESTAMP_START), tz = "GMT", format="%Y%m%d%H%M")) |>
    mutate(date = lubridate::as_date(time))
  
  # slice by years
  hhdf = hhdf |> filter(lubridate::year(date) %in% ystart:yend)
  
  
  message("- Add SW_OUT=NA if not present")
  if (!("SW_OUT" %in% colnames(hhdf))) {
    hhdf$SW_OUT = NA
  }
  
  
  # Aggregate to daily 24-hr means  ----------------------------------------------------------
  message("- downsampling FLUXNET format - 24 hr means")
  ddf_24hr_mean <-
    try(
      hhdf |>
        group_by(date) |>
        select(-TIMESTAMP_START, -TIMESTAMP_END) |>
        summarize_all(.funs = mean)
    )
  
  
  # Aggregate around daily maximum ppfd for acclimating model
  # ---------------------------------------------------------
  test.3day = hhdf %>% filter(date >= as_date(paste0(floor((ystart+yend)/2),"-06-01")) &
                                date <= as_date(paste0(floor((ystart+yend)/2),"-06-03")) )
  
  aggregate_daily_3hr_maxima = function(df){
    # Get the time at which SW_IN is maximum
    maxppfd <- df %>% filter(SW_IN_F_MDS == max(SW_IN_F_MDS))
    max_t <- maxppfd$time[1]
    
    # Select times that lie in 3-hr interval around max_t
    df_aroundmax <- df %>% filter(time < (max_t + 1.5*3600) &
                                    time > (max_t - 1.5*3600) )
    
    # take mean of selected entries
    df_mean <- df_aroundmax |>
      select(-TIMESTAMP_START, -TIMESTAMP_END) |>
      summarize_all(.funs = mean)
    
    df_mean
  }
  
  
  #'
  ## --------------------------------------------------------------------
  # Apply 3hr maxima aggregation to all data
  # ----------------------------------------
  message("- downsampling FLUXNET format - daily 3-hr means around max ppfd")
  ddf_3hr_maxima <- hhdf |>
    group_by(date) |>
    do(aggregate_daily_3hr_maxima(.)) |>
    ungroup()
  
  #'
  ## --------------------------------------------------------------------
  aggregate_daily_daylength = function(df){
    # Get the time at which SW_IN > 0
    pos_ppfd <- df %>% filter(SW_IN_F_MDS > 10)
    # if SW_IN is unavailable in that year calc daylength based on NETRAD
    if (nrow(pos_ppfd) < 2){
      pos_ppfd <- df %>% filter(NETRAD > 25)
    }
    
    tmax <- max(pos_ppfd$time)
    tmin <- min(pos_ppfd$time)
    
    # Select times that lie in 3-hr interval around max_t
    df_aroundmax <- df %>% filter(time <= tmax &
                                    time >= tmin )
    
    # take mean of selected entries
    df_mean <- df_aroundmax |>
      select(-TIMESTAMP_START, -TIMESTAMP_END) |>
      summarize_all(.funs = mean) |>
      mutate(daylength = difftime(tmax, tmin, units="hours") |> as.numeric())
    
    df_mean
  }
  
  #'
  ## --------------------------------------------------------------------
  # Apply daytime mean aggregation to all data
  # ------------------------------------------
  message("- downsampling FLUXNET format - daytime means")
  ddf_daytime_mean <- hhdf |>
    group_by(date) |>
    do(aggregate_daily_daylength(.)) |>
    ungroup()
  
  
  #'
  ## --------------------------------------------------------------------
  # Calculate daily tmax and tmin from hh data
  # ------------------------------------------
  tmaxmin <-
    hhdf |>
    group_by(date) |>
    summarize(
      tmax = max(TA_F_MDS),
      tmin = min(TA_F_MDS)
    )
  
  #'
  ## --------------------------------------------------------------------
  # Creating driver object  ------------------------------------------------------
  message("- compiling drivers")
  
  ddf_24hr_mean <- ddf_24hr_mean |> filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2))
  ddf_3hr_maxima <- ddf_3hr_maxima |> filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2))
  ddf_daytime_mean <- ddf_daytime_mean |> filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2))
  
  ccov = zenodo_driver[zenodo_driver$sitename == i,4][[1]][[1]]$ccov
  le = zenodo_driver[zenodo_driver$sitename == i,4][[1]][[1]]$le
  
  ccov = ccov[1:dim(ddf_24hr_mean)[1]]
  le = le[1:dim(ddf_24hr_mean)[1]]
  
  
  site_lon = meta[[1]]$longitude
  site_lat = meta[[1]]$latitude
  
  lonid = which(lons > site_lon)[1]-1
  latid = which(lats > site_lat)[1]-1
  n = 1
  S80_slice = S80[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site = mean(as.numeric(S80_slice, na.rm=T))
  whc_site_sd = sd(as.numeric(S80_slice, na.rm=T))
  
  
  p_hydro_drivers <- p_model_drivers
  
  p_hydro_drivers$sitename[[1]] = i
  
  
  p_hydro_drivers$site_info[[1]] =
    tibble(
      lon=meta[[1]]$longitude,
      lat=meta[[1]]$latitude,
      elv = meta[[1]]$elevation,
      canopy_height = ifelse(is.na(meta[[1]]$canopy_height), yes = 20, meta[[1]]$canopy_height),
      reference_height = ifelse(is.na(meta[[1]]$reference_height), yes = 22, meta[[1]]$reference_height),
      whc = whc_site,
      whc_sd = whc_site_sd,
      IGBP_veg_short = meta[[1]]$IGBP_veg_short
    )
  
  
  kfFEC = 2.04

  p_hydro_drivers$forcing <-
    ddf_24hr_mean |>
    left_join(tmaxmin) |>
    group_by(date) |>
    summarize(
      date = date,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = P_F * 48 /(60 * 60 * 24), # P_F [mm timestep-1] * 48 [timesteps day-1] / 86400 [secs day-1 ]
      tmin = tmin, # TMIN_F_MDS,
      tmax = tmax, # TMAX_F_MDS,
      fapar = FPAR,
      co2 = CO2_F_MDS
    ) |>
    mutate(ccov = ccov)  |>
    mutate(le = le)  |>
    
    list()
  
  p_hydro_drivers$forcing_acclim <-
    ddf_3hr_maxima |>
    left_join(tmaxmin) |>
    group_by(date) |>
    summarize(
      date = date,
      time = time,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = NA, # P_F * 48 / (60 * 60 * 24),
      tmin = tmin, # TMIN_F_MDS,
      tmax = tmax, # TMAX_F_MDS,
      fapar = FPAR,
      co2 = CO2_F_MDS
    ) |>
    mutate(ccov = ccov)  |>
    list()
  
  driver = rbind(driver,p_hydro_drivers)
}

####-----
# save
####-----

# write all drivers to file
# apply compression to minimize space
message(paste0("save driver in ",here("data","driver.rda")))
save(driver,
     file = here("data","driver.rda")
)
