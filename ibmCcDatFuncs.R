# This code was written by Isaac Brito-Morales (i.britomorales@uq.edu.au)
# Please do not distribute this code without permission.
# NO GUARANTEES THAT CODE IS CORRECT
# Caveat Emptor!
# Function's arguments
# ipath: directory where the netCDF files are located
# opath: directory to allocate the new regrided netCDF files
# resolution = resolution for the regrid process
regrid <- function(ipath, opath, resolution) {
  ####################################################################################
  ####### Defining the main packages
  ####################################################################################
  # List of pacakges that we will be used
  list.of.packages <- c("doParallel", "parallel", "stringr", "data.table")
  # If is not installed, install the pacakge
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  # Load packages
  lapply(list.of.packages, require, character.only = TRUE)
  ####################################################################################
  ####### Getting the path and directories for the files
  ####################################################################################
  # Establish the find bash command
  line1 <- paste(noquote("find"), noquote(ipath), "-type", "f", "-name",
                 noquote("*.nc"), "-exec", "ls", "-l", "{}")
  line2 <- paste0("\\", ";")
  line3 <- paste(line1, line2)
  # Getting a list of directories for every netCDF file
  dir_files <- system(line3, intern = TRUE)
  dir_nc <- strsplit(x = dir_files, split = " ")
  nc_list <- lapply(dir_nc, function(x){f1 <- tail(x, n = 1)})
  # Cleaning the directories to get a final vector of directories
  final_nc <- lapply(nc_list, function(x) {
    c1 <- str_split(unlist(x), pattern = "//")
    c2 <- paste(c1[[1]][1], c1[[1]][2], sep = "/")})
  files.nc <- unlist(final_nc)
  ####################################################################################
  ####### Starting the regrid process
  ####################################################################################
  # Resolution
  if(resolution == "1") {
    grd <- "r360x180"
  } else if(resolution == "0.5") {
    grd <- "r720x360"
  } else if(resolution == "0.25") {
    grd <- "r1440x720"
  }
  # Parallel looop
  UseCores <- 3 # we can change this number
  cl <- makeCluster(UseCores)
  registerDoParallel(cl)
  foreach(j = 1:length(files.nc), .packages = c("stringr")) %dopar% {
    # Trying to auto the name for every model
    var_obj <- system(paste("cdo -showname", files.nc[j]), intern = TRUE)
    var_all <- str_replace_all(string = var_obj, pattern = " ", replacement = "_")
    var <- tail(unlist(strsplit(var_all, split = "_")), n = 1)
    # Running CDO regrid
    system(paste(paste("cdo -remapbil,", grd, ",", sep = ""),
                 paste("-selname",var, sep = ","), files.nc[j],
                 paste0(opath, basename(files.nc[j])), sep = (" "))) # -P 2
  }
  stopCluster(cl)
}


# EXAMPLE Run the regrid R function
# regrid(ipath = "/data/ClimateModels/",
#        opath = "/data/ClimateModelsRegrid/",
#        resolution = "0.5")


########## ********** ########## 

# Arguments
# ipath: directory where the netCDF files are located
# opath1: directory to allocate the split files
# opath2: directory to allocate the vertical average files
olayer <- function(ipath, opath1, opath2) {
  ####################################################################################
  ####### Defining the main packages (tryining to auto this)
  ####################################################################################
  # List of pacakges that we will be used
  list.of.packages <- c("doParallel", "parallel", "stringr", "data.table")
  # If is not installed, install the pacakge
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  # Load packages
  lapply(list.of.packages, require, character.only = TRUE)
  ####################################################################################
  ####### Getting the path and directories for the files
  ####################################################################################
  # Establish the find bash command
  line1 <- paste(noquote("find"), noquote(ipath), "-type", "f", "-name",
                 noquote("*.nc"), "-exec", "ls", "-l", "{}")
  line2 <- paste0("\\", ";")
  line3 <- paste(line1, line2)
  # Getting a list of directories for every netCDF file
  dir_files <- system(line3, intern = TRUE)
  dir_nc <- strsplit(x = dir_files, split = " ")
  nc_list <- lapply(dir_nc, function(x){f1 <- tail(x, n = 1)})
  # Cleaning the directories to get a final vector of directories
  final_nc <- lapply(nc_list, function(x) {
    c1 <- str_split(unlist(x), pattern = "//")
    c2 <- paste(c1[[1]][1], c1[[1]][2], sep = "/")})
  files.nc <- unlist(final_nc)
  ####################################################################################
  ####### Filtering by layers and generating new netCDF files with outputs
  ####################################################################################
  # Parallel looop
  cl <- makeCluster(3)
  registerDoParallel(cl)
  foreach(j = 1:length(files.nc), .packages = c("stringr")) %dopar% {
    # Trying to auto the name for every model
    var_obj <- system(paste("cdo -showname", files.nc[j]), intern = TRUE)
    var_all <- str_replace_all(string = var_obj, pattern = " ", replacement = "_")
    var <- tail(unlist(strsplit(var_all, split = "_")), n = 1)
    # Defining depths
    levels <- as.vector(system(paste("cdo showlevel", files.nc[j]), intern = TRUE))
    lev <- unlist(strsplit(levels, split = " "))
    depths <- unique(lev[lev != ""])
    # Some can come in cm
    if(depths[1] >= 50) {
      sf <- depths[as.numeric(depths) <= 500]
      ep <- depths[as.numeric(depths) >= 0 & as.numeric(depths) <= 20000]
      mp <- depths[as.numeric(depths) > 20000 & as.numeric(depths) <= 100000]
      bap <- depths[as.numeric(depths) > 100000]
    } else {
      sf <- depths[as.numeric(depths) <= 5]
      ep <- depths[as.numeric(depths) >= 0 & as.numeric(depths) <= 200]
      mp <- depths[as.numeric(depths) > 200 & as.numeric(depths) <= 1000]
      bap <- depths[as.numeric(depths) > 1000]
    }
    # Running CDO
    # Surface
    system(paste(paste("cdo -L -sellevel,",
                       paste0(sf, collapse = ","), ",", sep = ""),
                 paste("-selname,", var, sep = ""), files.nc[j],
                 paste0(opath1, "01-sf_", basename(files.nc[j]))))
    # Epipelagic
    system(paste(paste("cdo -L -sellevel,",
                       paste0(ep, collapse = ","), ",", sep = ""),
                 paste("-selname,", var, sep = ""), files.nc[j],
                 paste0(opath1, "02-ep_", basename(files.nc[j]))))
    # Mesopelagic
    system(paste(paste("cdo -L -sellevel,",
                       paste0(mp, collapse = ","), ",", sep = ""),
                 paste("-selname,", var, sep = ""), files.nc[j],
                 paste0(opath1, "03-mp_", basename(files.nc[j]))))
    # Bathypelagic
    system(paste(paste("cdo -L -sellevel,",
                       paste0(bap, collapse = ","), ",", sep = ""),
                 paste("-selname,", var, sep = ""), files.nc[j],
                 paste0(opath1, "04-bap_", basename(files.nc[j]))))
  }
  stopCluster(cl)
  ####################################################################################
  ####### Getting the path and directories for the "split by depth" files
  ####################################################################################
  # Establish the find bash command
  line1.1 <- paste(noquote("find"), noquote(opath1), "-type", "f", "-name",
                   noquote("*.nc"), "-exec", "ls", "-l", "{}")
  line2.1 <- paste0("\\", ";")
  line3.1 <- paste(line1.1, line2.1)
  # Getting a list of directories for every netCDF file
  dir_files.2 <- system(line3.1, intern = TRUE)
  dir_nc.2 <- strsplit(x = dir_files.2, split = " ")
  nc_list.2 <- lapply(dir_nc.2, function(x){f1 <- tail(x, n = 1)})
  # Cleaning the directories to get a final vector of directories
  final_nc.2 <- lapply(nc_list.2, function(x) {
    c1 <- str_split(unlist(x), pattern = "//")
    c2 <- paste(c1[[1]][1], c1[[1]][2], sep = "/")})
  files.nc.2 <- unlist(final_nc.2)
  ####################################################################################
  ####### Filtering by layers and generating the "weighted-average depth layer"
  ####################################################################################
  # Parallel looop
  cl <- makeCluster(3)
  registerDoParallel(cl)
  foreach(i = 1:length(files.nc.2), .packages = c("stringr")) %dopar% {
    # Running CDO
    system(paste(paste("cdo -L vertmean", sep = ""), files.nc.2[i],
                 paste0(opath2, basename(files.nc.2[i])), sep = (" ")))
  }
  stopCluster(cl)
}


# EXAMPLE Running the ocean depth layer R function will:
#   olayer(ipath = "/data/ClimateModelsRegrid/",
#          opath1 = "/data/ClimateModelsRegridLayer/",
#          opath2 = "/data/ClimateModelsRegridLayerMean/")

########## ********** ########## 

merge_files <- function(ipath, opath1) {
  ####################################################################################
  ####### Defining the main packages (tryining to auto this)
  ####################################################################################
  # List of pacakges that we will be used
  list.of.packages <- c("doParallel", "parallel", "stringr", "data.table")
  # If is not installed, install the pacakge
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  # Load packages
  lapply(list.of.packages, require, character.only = TRUE)
  ####################################################################################
  ####### Getting the path and directories for the files
  ####################################################################################
  # Establish the find bash command
  line1 <- paste(noquote("find"), noquote(ipath), "-type", "f", "-name",
                 noquote("*.nc"), "-exec", "ls", "-l", "{}")
  line2 <- paste0("\\", ";")
  line3 <- paste(line1, line2)
  # Getting a list of directories for every netCDF file
  dir_files <- system(line3, intern = TRUE)
  dir_nc <- strsplit(x = dir_files, split = " ")
  nc_list <- lapply(dir_nc, function(x){f1 <- tail(x, n = 1)})
  # Cleaning the directories to get a final vector of directories
  final_nc <- lapply(nc_list, function(x) {
    c1 <- str_split(unlist(x), pattern = "//")
    c2 <- paste(c1[[1]][1], c1[[1]][2], sep = "/")})
  files.nc <- unlist(final_nc)
  ####################################################################################
  ####### Filtering by layers and generating new netCDF files with outputs
  ####################################################################################
  # Filtering (not dplyr!) by ocean layers
  sf <- files.nc[str_detect(string = basename(files.nc), pattern = "01-sf*") == TRUE]
  ep <- files.nc[str_detect(string = basename(files.nc), pattern = "02-ep*") == TRUE]
  mp <- files.nc[str_detect(string = basename(files.nc), pattern = "03-mp*") == TRUE]
  bap <- files.nc[str_detect(string = basename(files.nc), pattern = "04-bap*") == TRUE]
  # Defining how many models are per ocean layer
  model_list_sf <- lapply(sf, function(x)
  {d1 <- unlist(strsplit(x = basename(x), split = "_"))[4]})
  model_list_ep <- lapply(ep, function(x)
  {d1 <- unlist(strsplit(x = basename(x), split = "_"))[4]})
  model_list_mp <- lapply(bap, function(x)
  {d1 <- unlist(strsplit(x = basename(x), split = "_"))[4]})
  model_list_bap <- lapply(bap, function(x)
  {d1 <- unlist(strsplit(x = basename(x), split = "_"))[4]})
  models <- unique(unlist(c(model_list_sf, model_list_ep, model_list_mp, model_list_bap)))
  # Parallel looop
  cl <- makeCluster(3)
  registerDoParallel(cl)
  foreach(i = 1:length(models), .packages = c("stringr")) %dopar% {
    f1 <- ep[str_detect(string = basename(sf), pattern = models[i]) == TRUE]
    system(paste(paste("cdo -L mergetime", paste0(f1, collapse = " "), sep = " "),
                 paste0(opath1, paste(unlist(strsplit(basename(f1[1]), "_"))[c(1:7)],
                                      collapse = "_"), ".nc"), sep = (" ")))
    f2 <- ep[str_detect(string = basename(ep), pattern = models[i]) == TRUE]
    system(paste(paste("cdo -L mergetime", paste0(f2, collapse = " "), sep = " "),
                 paste0(opath1, paste(unlist(strsplit(basename(f2[1]), "_"))[c(1:7)],
                                      collapse = "_"), ".nc"), sep = (" ")))
    f3 <- mp[str_detect(string = basename(mp), pattern = models[i]) == TRUE]
    system(paste(paste("cdo -L mergetime", paste0(f3, collapse = " "), sep = " "),
                 paste0(opath1, paste(unlist(strsplit(basename(f3[1]), "_"))[c(1:7)],
                                      collapse = "_"), ".nc"), sep = (" ")))
    f4 <- bap[str_detect(string = basename(bap), pattern = models[i]) == TRUE]
    system(paste(paste("cdo -L mergetime", paste0(f4, collapse = " "), sep = " "),
                 paste0(opath1, paste(unlist(strsplit(basename(f4[1]), "_"))[c(1:7)],
                                      collapse = "_"), ".nc"), sep = (" ")))
  }
  stopCluster(cl)
}

# EXAMPLE Running the merge function
# merge_files(ipath = "/Users/bri273/Desktop/CDO/models_regrid_vertmean/",
#             opath1 = "/Users/bri273/Desktop/CDO/models_regrid_zmerge/")

