### mapping functions for site assessments
### written by Jessica Couture, Research Scientist
### Marine SPARC, Conservation International


##### Libraries needed
library(tidyverse)
library(mregions)
library(googledrive)
library(sf)
library(kableExtra)
library(patchwork)
library(ncdf4)
library(viridis)
library(weathermetrics)
library(chron)
library(lattice)
library(lubridate)


##################################################################
########################### Extent Map ###########################
##################################################################
# mpaShp = sf object of MPA shape file
# siteName = string representing the site to be mapped
# mrSearch = string to search marine regions DB for EEZ shp file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# zoomLat = latitudinal buffer for the entire extent map in decimal degrees (size depends on context extent map scale)
# zoomLon = longitudinal buffer for the entire extent map in decimal degrees (size depends on context extent map scale)

# output: ggplot map object with MPA shapefile, plot area and surrounding countries plotted


extMapFunc<-function(mpaShp,buff,zoomLat,zoomLon,corr=c(NA,"lonFx")){
  
  if(is.na(corr)) {
   
     mpa<-mpaShp%>%
      st_set_crs(4326)#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    
  } else {
    mpa<-mpaShp%>%
      st_set_crs(4326)%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaShp)[1])%>%
      st_set_crs(4326)
  }
  
  siteBB<-st_bbox(mpa)
  
  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  } else {
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,
                                              180+(180+siteBB$xmin),
                                              siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,
                                              180+(180+siteBB$xmax),
                                              siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  }
  
  bbDf<-data.frame(lon=c(as.numeric(ploArea$xmin),
                         as.numeric(ploArea$xmax),
                         as.numeric(ploArea$xmax),
                         as.numeric(ploArea$xmin)),
                   lat=c(as.numeric(ploArea$ymin),
                         as.numeric(ploArea$ymin),
                         as.numeric(ploArea$ymax),
                         as.numeric(ploArea$ymax)))%>%
   sfheaders::sf_polygon()%>%
    st_set_crs(4326)
     # st_as_sf(.,
    #          coords=c("lon","lat"),
    #          crs=st_crs(mpa))%>%
    # mutate(id=1)%>%
    # group_by(id)%>%
    # # summarise()%>%
    # st_convex_hull()
  # st_cast("POLYGON")
  
  
  
  if(is.na(corr)) {
    totArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-zoomLon),
                       xmax=as.numeric(siteBB$xmax+zoomLon),
                       ymin=as.numeric(siteBB$ymin-zoomLat),
                       ymax=as.numeric(siteBB$ymax+zoomLat)),crs = st_crs(mpa))
  } else {
    totArea<-st_bbox(c(xmin=ifelse(as.numeric(siteBB$xmin)<0,
                                   180+(180+as.numeric(siteBB$xmin)),
                                   as.numeric(siteBB$xmin))-zoomLon,
                       xmax=ifelse(as.numeric(siteBB$xmax)<0,
                                   180+(180+as.numeric(siteBB$xmax)),
                                   as.numeric(siteBB$xmax))+zoomLon,
                       ymin=as.numeric(siteBB$ymin)-zoomLat,
                       ymax=as.numeric(siteBB$ymax)+zoomLat),
                     crs = st_crs(mpa))
  }
  
  ### crop base layers from world map
  if(is.na(corr)){
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      group_by(group) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast(to="POLYGON")%>%
      st_crop(.,totArea) 
    
  } else {
    
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      filter(lon2>totArea$xmin & lon2<totArea$xmax,
             lat>totArea$ymin & lat<totArea$ymax)%>%
      # st_drop_geometry()%>%
      group_by(subregion)%>%
      mutate(dup=n()>3)%>%
      filter(dup==TRUE)%>%
      ungroup()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id="group")%>%
      # st_make_valid()%>%
      st_set_crs(4326)
    
  }
  
  ### pull EEZ shp file from the marine regions data
  # having a lot of trouble with the mregions package, removing for now
  
  # mrEezName<-mr_geo_code(place = mrSearch, like = TRUE, fuzzy = FALSE)%>%
  #   filter(placeType=="EEZ")%>%
  #   select(accepted)%>%
  #   as.numeric()
  # 
  # 
  # if(is.na(corr)) {
  #   
  #   siteEez<-mr_shp(key="MarineRegions:eez",filter=8429)%>%
  #     st_as_sf(.,
  #              coords=c("long","lat"),
  #              crs=4326)
  #   
  # } else {
  #   siteEez<-mr_shp(key="MarineRegions:eez",filter=mrEezName)%>%
  #     st_as_sf(.,
  #              coords=c("long","lat"),
  #              crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
  #     st_cast(to="POINT")%>% # convert polygons to points
  #     mutate(lon=st_coordinates(geometry)[,"X"] ,
  #            lat=st_coordinates(geometry)[,"Y"])%>%
  #     mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
  #                           TRUE ~ lon))%>%
  #     st_drop_geometry()%>%
  #     sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(siteEez)[1])
  # }

  extMap<-ggplot()+
    geom_sf(data=geoShp)+
    # geom_sf(data=siteEez,fill="grey33",color=NA,alpha=0.3)+
    geom_sf(data=mpa,fill="deepskyblue3",alpha=0.3,color="deepskyblue4",linetype=2)+
    geom_sf(data=bbDf,fill=NA,color="darkorange3")+
    annotate(geom="text", x=ploArea$xmin+5, y=ploArea$ymin-1, label="Area assessed",
             color="darkorange3")+
    scale_x_continuous(expand = c(0,0),limits = c(totArea$xmin,totArea$xmax),breaks=seq(round(totArea$xmin,digits = 0.1),round(totArea$xmax,digits = 0.1),by=10))+
    scale_y_continuous(expand = c(0,0),limits = c(totArea$ymin,totArea$ymax))+
    theme_bw()+
    labs(title = "Extent map",x="longitude",y="latitude")
  
  return(extMap)
  
}



##################################################################
##################### Temperatures at depth ######################
########################### Static ###############################

# mpaShp = sf object of MPA shape file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# targDepth = one of "epi" for epipelagic, "meso" for mesoplagic, "bathy" for bathypelagic to indicate which dataset to pull from the drive
# YEAR = year in the future for analysis from 2020 - 2100, default is 2050
# corr = indicates whether longitude (lonFx) needs to be corrected for crossing the 180 meridian, if no use NA

# output: ggplot map object with MPA shapefile and projected temperature and indicated depth range in 2050



temp3dPlots<-function(mpaShp,buff,targDepth=c("epi","meso","bathy"),YEAR = 2050,corr=c(NA, "lonFx")){
  
  
  ### set bounding box for target area
 if(is.na(corr)) {
   
   mpa<-mpaShp%>%
    st_set_crs(4326)
   
 } else {
    
   mpa<-mpaShp%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaShp)[1])%>%
      st_set_crs(4326)
  }
  
  siteBB<-st_bbox(mpa)
  
  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  } else {
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,180+(180+siteBB$xmin),siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,180+(180+siteBB$xmax),siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  }
  
  
  ### crop base layers from world map
  
  if(is.na(corr)) {
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      group_by(group) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast("POLYGON")%>%
      st_crop(.,ploArea)
    
  } else {
    
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      filter(lon2>ploArea$xmin & lon2<ploArea$xmax,
             lat>ploArea$ymin & lat<ploArea$ymax)%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id="group")%>%
      st_make_valid()
  }
  

  
  ### pull in ensamble data based on target depth region
  
  datFile<-"tempThetaoProj.rds" #where to put it
  
  driveFile<-ifelse(targDepth=="epi","02-ep_thetao_AEMean_ssp585_r1i1p1f1_2015-2100.rds",
                    ifelse(targDepth=="meso","03-mp_thetao_AEMean_ssp585_r1i1p1f1_2015-2100.rds",
                           "04-bap_thetao_AEMean_ssp585_r1i1p1f1_2015-2100.rds")) 
  
  drive_download(driveFile, #raster file
                 path=datFile,
                 overwrite = TRUE)
  
  rastT<-read_rds(datFile) # raster file
  
  
  ### delete large file once the object has been created
  if (file.exists(datFile)) {
    unlink(datFile)
    cat(paste("The file has been deleted: ",datFile))
  }
  
  ### Isolate corresponding raster layers based on year defined
  indxN<-1:12 #2015, now
  indxF<-((YEAR-2015)*12+1):((YEAR-2015)*12+12) # pull data for the future
  
  ### pull 2015 and future (YEAR) data layers
  
  rastYrN<-rastT[[indxN]] ### pull just the data for 2015
  rastYrF<-rastT[[indxF]] ### pull just the data for given year in the future
  
  ### calculate average annual temps from 12mo stack (returns one raster for each year)
  
  rastAnnN<-raster::overlay(rastYrN,fun=mean)
  rastAnnF<-raster::overlay(rastYrF,fun=mean)
  
  ### calculate change in temp
  
  dTempRast<-overlay(rastAnnF,rastAnnN,fun=function(x,y) x-y)
  
  ### convert raster to dataframe for ggplotting
  
  if(is.na(corr)) {
    
    dTempSpDf<-crop(dTempRast,ploArea)%>% # crop to intended area
      as.data.frame(.,xy=TRUE)%>% ## convert to df
      rename(lon=x,
             lat=y)
    
  } else {
    
    dTempSpDf<-dTempRast%>%
      as.data.frame(.,xy=TRUE)%>%
      mutate(lon=case_when(x<0 ~ 180+(180+x),
                           TRUE ~ x))%>%
      rename(lat=y)%>%
      filter(lon>ploArea["xmin"],
             lon<ploArea["xmax"],
             lat>ploArea["ymin"],
             lat<ploArea["ymax"])
    
  }
  

  
  #### plot
  
  subtit<-paste("2015 -",YEAR,sep=" ")
  
  deltaTemp <- ggplot()+
    geom_raster(data=dTempSpDf,aes(x=lon,y=lat,fill=layer))+
    geom_sf(data = geoShp,fill="grey",color="black")+
    geom_sf(data = mpa,fill=NA, color="white")+
    scale_fill_gradient2(low="navyblue",mid="white",high="firebrick")+ #limits=c(round(min(dTempRast$ ??)-5,digits=-1)l,round(max(dTempRast$ ??),digits=-1))
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
          plot.subtitle = element_text(size=13, color = "blue4",hjust = 0.5))+
    labs(title="Change in ocean temperature",subtitle = subtit)
  
  return(deltaTemp)
  
}


####################################################################
############ Climate projections from NOAA CMIP6 data ##############
######################### SST, SSS, O2 ############################

# mpaShp = sf object of MPA shape file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# siteName = short code for site name to be used in the name of the output file (string with no spaces)
# var = one of "tos" (SST), "sos" (SSSal), "o2" ([oxygen]) variables to plot.SST and SSS data are by definition only calculated at one depth (0m) and oxygen concentration is assessed at shallow and deep 
# corr = indicates whether longitude (lonFx) needs to be corrected for crossing the 180 meridian, if no use NA

# output: returns 1 or 2 plots of the projected metrics around the input MPA site


noaaCMIP6<-function(mpaShp,
                    buff,
                    siteName,
                    var=c("tos","sos","o2"),
                    corr=c(NA,"lonFx")) {

  ### pull in global data from the drive based on indicated variable (var)
  datFile<-"data/tempDat/noaaTmp.nc"
  
  driveFile<-if(var=="tos") {
    "tos_Omon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc"
  } else {
    if(var=="sos"){
      "sos_Omon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc"
    } else {
      "o2_Oyr_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_2015-2100.nc"
    }
  }
  
  drive_download(driveFile,
                 path=datFile,
                 overwrite = TRUE)
  
  ncDat<-nc_open(datFile)
  
  
  ### pull variables from nc file: ncDat
  lon<-ncvar_get(ncDat,"longitude")
  lat<-ncvar_get(ncDat,"latitude")
  datVar<-ncvar_get(ncDat,var) # pull variable data
  timeBnd<-ncvar_get(ncDat,"time_bnds")
  
  ttl<-ncatt_get(ncDat,var,"long_name")
  ttlNm<-as.character(ttl$value)
  
  ### delete large file from local
  
  if (file.exists(datFile)) {
    unlink(datFile)
    cat(paste("The file has been deleted: ",datFile))
  }
  
  ### convert arrays to DF
  
  if(var %in% c("tos","sos")){
    
    newYr<-seq(1,ncol(timeBnd),by=12)
    aveDf<-NULL
    datesDf<-tibble(date=as_date(timeBnd[1,],origin="1850-01-01"))%>%
      mutate(year=year(date))
    
    
    for(i in 1:length(newYr)){
      
      yrDf<-expand.grid(lat=lat[1,],lon=lon[,1])%>%
        mutate(lon=ifelse(lon>180,-(360-lon),lon))
      
      
      for(j in 1:12){ # iterate over time
        k=newYr[i]+j-1
        
        monDf <- expand.grid(lat=lat[1,],lon=lon[,1])%>%
          mutate(lon=ifelse(lon>180,-(360-lon),lon),
                 varVal=as.vector(t(datVar[,,k])))
        
        yrDf = yrDf %>% bind_cols(monDf%>%select(varVal))    
        
      }
      yrDf<-yrDf%>%
        mutate(aveVar=rowMeans(.[-1:-2],na.rm = TRUE))%>%
        select(lat,lon,aveVar)%>%
        mutate(year=2014+i)#datesDf[k,"year"])
      
      aveDf<-aveDf%>%
        bind_rows(.,yrDf)
      
      i=i+1
    }
  } else {
    
    o2Dfshallow<-NULL
    
    datesDf<-tibble(date=as_date(timeBnd[1,],origin="1850-01-01"))%>%
      mutate(year=year(date))
    
    for (i in 1:nrow(datesDf)){ # iterate over years
      
      tempDf <- expand.grid(lat=lat[1,],lon=lon[,1])%>%
        mutate(lon=ifelse(lon>180,-(360-lon),lon),
               o2=as.vector(t(datVar[,,3,i])))%>% # 3 is depth of 25m
        mutate(year=as.numeric(datesDf[i,"year"]))
      
      o2Dfshallow = o2Dfshallow %>% bind_rows(tempDf)
    }
    
    ## Deeper: >425m
    
    o2DfDeep<-NULL
    
    for (i in 1:nrow(datesDf)){ # iterate over years
      
      tempDf <- expand.grid(lat=lat[1,],lon=lon[,1])%>%
        mutate(lon=ifelse(lon>180,-(360-lon),lon),
               o2=as.vector(t(datVar[,,26,i])))%>% # 26 is depth of 427.3156m
        mutate(year=as.numeric(datesDf[i,"year"]))
      
      o2DfDeep = o2DfDeep %>% bind_rows(tempDf)
    }
  }
  
  ### calculate changes to 2060 and 2100
  if(var %in% c("tos","sos")) {
    
    varDiffs<-aveDf%>%
      filter(year==2020)%>%
      rename(var2020=aveVar)%>%
      bind_cols(aveDf%>%
                  filter(year==2060)%>%
                  select(var2060=aveVar),
                aveDf%>%
                  filter(year==2100)%>%
                  select(var2100=aveVar))%>%
      mutate(diff2060=var2060-var2020,
             diff2100=var2100-var2020)
  } else {
    
    o2DiffsSh<-o2Dfshallow%>%
      filter(year==2020)%>%
      rename(o22020=o2)%>%
      bind_cols(o2Dfshallow%>%
                  filter(year==2060)%>%
                  select(o22060=o2),
                o2Dfshallow%>%
                  filter(year==2100)%>%
                  select(o22100=o2))%>%
      mutate(diff2060=o22060-o22020,
             diff2100=o22100-o22020)
    
    o2DiffsDp<-o2DfDeep%>%
      filter(year==2020)%>%
      rename(o22020=o2)%>%
      bind_cols(o2DfDeep%>%
                  filter(year==2060)%>%
                  select(o22060=o2),
                o2DfDeep%>%
                  filter(year==2100)%>%
                  select(o22100=o2))%>%
      mutate(diff2060=o22060-o22020,
             diff2100=o22100-o22020)
    
  }
  
  ### set up plot area and shp files
  
  if(is.na(corr)) {
    
    mpa<-mpaShp%>%
      st_set_crs(4326)
    
    siteBB<-st_bbox(mpa)
    
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  } else {
    
    mpa<-mpaShp%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaShp)[1])%>%
      st_set_crs(4326)
    
    siteBB<-st_bbox(mpa)
    
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,180+(180+siteBB$xmin),siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,180+(180+siteBB$xmax),siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  }
  
  ### Clip var data to site
  
  if(is.na(corr)) {
    if(var %in% c("tos","sos")){
    
    varSite<-varDiffs%>%
      filter(lon>as.numeric(ploArea$xmin),
             lon<as.numeric(ploArea$xmax),
             lat>as.numeric(ploArea$ymin),
             lat<as.numeric(ploArea$ymax))
    } else {
      o2SiteSh<-o2DiffsSh%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)
      
      o2SiteDp<-o2DiffsDp%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)
    }
  } else {
    
    if(var %in% c("tos","sos")){
      
      varSite<-varDiffs%>%
        mutate(lon=case_when(lon < 0 ~ 180+(180+lon),
                             TRUE ~ lon))%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))
    } else {
      o2SiteSh<-o2DiffsSh%>%
        mutate(lon=case_when(lon < 0 ~ 180+(180+lon),
                             TRUE ~ lon))%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)
      
      o2SiteDp<-o2DiffsDp%>%
        mutate(lon=case_when(lon < 0 ~ 180+(180+lon),
                             TRUE ~ lon))%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)
    }
  }
  
  ### crop base layers from world map
  
  if(is.na(corr)) {
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      group_by(group) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast("POLYGON")%>%
      st_crop(.,ploArea)
    
  } else {
    
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      filter(lon2>ploArea$xmin & lon2<ploArea$xmax,
             lat>ploArea$ymin & lat<ploArea$ymax)%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id="group")%>%
      st_make_valid()%>%
      st_set_crs(4326)
  }
  
  ## get eez
  
  # mrEezName<-mr_geo_code(place = mrSearch, like = TRUE, fuzzy = FALSE)%>%filter(placeType=="EEZ")%>%select(accepted)%>%as.numeric()
  # 
  # siteEez<-mr_shp(key="MarineRegions:eez",filter=mrEezName)%>%
  #   st_as_sf(.,
  #            crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  # 
  # if(corr=="lonFx") {
  #   siteEez<-siteEez%>%
  #     st_as_sf(.,
  #              coords=c("long","lat"),
  #              crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
  #     st_cast(to="POINT")%>% # convert polygons to points
  #     mutate(lon=st_coordinates(geometry)[,"X"] ,
  #            lat=st_coordinates(geometry)[,"Y"])%>%
  #     mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
  #                           TRUE ~ lon))%>%
  #     st_drop_geometry()%>%
  #     sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(siteEez)[1])
  # }
  
  ### plot data
  if(var %in% c("tos","sos")) {
    
    if(any(varSite$diff2060 > 0) & any(varSite < 0)) {
      sclFllGrd<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                      name=ifelse(var=="tos","degC","ppt"))
    } else {
      
      if(all(varSite$diff2060>=0)) {
        sclFllGrd<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                        name=ifelse(var=="tos","degC","ppt"))
      } else {
        sclFllGrd<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                        name=ifelse(var=="tos","degC","ppt"))
      }
    }
    
    sitePlt<-ggplot()+
      geom_tile(data=varSite,aes(x=lon,y=lat,fill=diff2060),width=1,height=1)+
      geom_sf(data=mpa,color="grey21",fill=NA)+
      geom_sf(data=geoShp)+
      # geom_sf(data=mrRev)+
      sclFllGrd+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
      theme(legend.direction = "vertical",legend.position = "right")+
      ggtitle(paste("Projected Change in Sea Surface",ifelse(var=="tos","Temperature","Salinity")),subtitle = "2020 to 2060")+
      theme_bw()
    
    return(sitePlt)
  } else {
    
    if(any(o2SiteSh$diff2060 > 0) & any(o2SiteSh < 0)) {
      sclFllGrdSh<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                      name="% change in\n[oxygen]")
    } else {
      
      if(all(o2SiteSh$diff2060>=0)) {
        sclFllGrdSh<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                       name="% change in\n[oxygen]")
      } else {
        sclFllGrdSh<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                       name="% change in\n[oxygen]")
      }
    }
    
    o2DiffSh<-ggplot()+
      geom_tile(data=o2SiteSh,aes(x=lon,y=lat,fill=o2perc60),width=1,height=1)+
      geom_sf(data=mpa,color="white",fill=NA)+
      geom_sf(data=geoShp)+
      # geom_sf(data=mrRev)+
      sclFllGrdSh+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
      theme(legend.direction = "vertical",legend.position = "right")+
      ggtitle("Projected Change in Oxygen Concentration: shallow (25m)",subtitle = "2020 to 2060")+
      theme_bw()
    
    if(any(o2SiteDp$diff2060 > 0) & any(o2SiteDp < 0)) {
      sclFllGrdDp<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                        name="% change in\n[oxygen]")
    } else {
      
      if(all(o2SiteDp$diff2060>=0)) {
        sclFllGrdDp<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                         name="% change in\n[oxygen]")
      } else {
        sclFllGrdDp<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                         name="% change in\n[oxygen]")
      }
    }
    
    o2DiffDp<-ggplot()+
      geom_tile(data=o2SiteDp,aes(x=lon,y=lat,fill=o2perc60),width=1,height=1)+
      geom_sf(data=mpa,color="white",fill=NA)+
      geom_sf(data=geoShp)+
      # geom_sf(data=mrRev)+
      sclFllGrdDp+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
      theme(legend.direction = "vertical",legend.position = "right")+
      ggtitle("Projected Change in Oxygen Concentration: Deep (>425m)",subtitle = "2020 to 2060")+
      theme_bw()
    
    return(list(o2DiffSh,o2DiffDp))
  }
  
}



###################################################################
############ pulling local species range data from ###############
########################## Aquamaps ##############################

# mpaShp = sf object of MPA shape file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# siteName = short code for site name to be used in the name of the output file (string with no spaces)
# nowFut = one of "now" for data from 2019 or "fut" for data projections, to indicate which dataset to pull from 
# corr = indicates whether longitude (lonFx) needs to be corrected for crossing the 180 meridian, if no use NA

# output: new csv file with species that intsect with the indicated site area (shapefile + buffer)
# !!! NOT coded for sites that cross the 180 meridian (like the Lau Seascape, see LauCC for this code)

aqmapClip<-function(mpaShp,
                    buff,
                    siteName,
                    nowFut=c("now","fut"),
                    corr=c(NA,"lonFx")) {
  
  mpa<-mpaShp%>%
    st_set_crs("EPSG 4326")#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  ### indicate which data file
  if(nowFut=="now"){
    aqData<-"hcaf_species_native-001.csv"
  } else {
    aqData<-"rcp85_hcaf_species_native_2050.csv"
  }
  
  ### download the data
  datFile<-"aquamapsDat.csv" #where to put it

  drive_download(aqData, #csv file
                     path=datFile,
                     overwrite = TRUE)

  aqDat<-read_csv(datFile,
                  col_select = c(SpeciesID,CenterLat,CenterLong,Probability))

  ### delete large file once the object has been created
  
  if (file.exists(datFile)) {
   unlink(datFile)
   cat(paste("The file has been deleted: ",datFile))
  }
  
  ### establish site coordinates
  siteBB<-st_bbox(mpaShp)
  
  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  } else {
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,180+(180+siteBB$xmin),siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,180+(180+siteBB$xmax),siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  }
  
  ### Clip aquamaps data to site
  
  if(is.na(corr)) {
    
    aqSiteDat<-aqDat%>%
        filter(CenterLong>as.numeric(ploArea$xmin),
               CenterLong<as.numeric(ploArea$xmax),
               CenterLat>as.numeric(ploArea$ymin),
               CenterLat<as.numeric(ploArea$ymax))%>%
        st_as_sf(coords=c("CenterLong","CenterLat"),
                 crs=4326)
    
  } else {
    
    aqSiteDat<-aqDat%>%
      mutate(lon=case_when(CenterLong < 0 ~ 180+(180+CenterLong),
                            TRUE ~ CenterLong))%>%
      rename(lat=CenterLat)%>%
      filter(lon>ploArea$xmin & lon<ploArea$xmax,
             lat>ploArea$ymin & lat<ploArea$ymax)%>%
        st_as_sf(coords=c("lon","lat"),
                 crs=4326)
    
  }
  
  fileNm<-paste("data/aqClips/",siteName,"Aq",nowFut,".shp",sep="")
  
  ### write to new file
  st_write(aqSite,fileNm)
  
  # print(fileNm)
  # return(aqSiteDat)
  
}

##################################################################
################## Changes to Species Richness ###################
########################## Aquamaps ##############################

# ***** requires files created by the aqmapClip function above. So must be run AFTER that function. *****

# mpaShp = sf object of MPA shape file
# siteName = siteName as used in the aqMapClip function to create the site specific data file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# thresh = threshold between 0-1 to apply to data (usually one of 0.4, 0.5, 0.6, 0.7,default is 0.3)
# YEAR = year in the future for analysis from 2020 - 2100, default is 2050
# corr = indicates whether longitude (lonFx) needs to be corrected for crossing the 180 meridian, if no use NA

# output ggplot map object with MPA shapefile, plot area and surrounding countries plotted
# !!! NOT coded for sites that cross the 180 meridian (like the Lau Seascape, see LauCC for this code)


dSppRich<-function(mpaShp,siteName,buff,thresh=0.4,corr=c(NA, "lonFx")){
  
  ## pull in & format maps data
  if(is.na(corr)) {
    
    mpa<-mpaShp%>%
      st_set_crs(4326)#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    
  } else {
    mpa<-mpaShp%>%
      st_set_crs(4326)%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaShp)[1])%>%
      st_set_crs(4326)
  }
  
  siteBB<-st_bbox(mpa)
  
  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  } else {
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,
                                              180+(180+siteBB$xmin),
                                              siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,
                                              180+(180+siteBB$xmax),
                                              siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  }
  
  aqSpp<-read_csv("data/speciesoccursum.csv")
  
  ### current data
  nowFile<-paste("data/aqClips/",siteName,"Aqnow",".shp",sep="")
  if (file.exists(nowFile)) {
    
    siteNow<-read_sf(nowFile)
    
  } else {
    stop(paste("The file does not exist: ",nowFile,". Create this file with 'aqmapClip()' before re-running"))
  }
  
  colnames(siteNow)<-c("SpeciesID","probability","geometry")
  
  ### future data
  futFile<-paste("data/aqClips/",siteName,"Aqfut",".shp",sep="")
  
  if (file.exists(futFile)) {
    
    siteFut<-read_sf(futFile)
    
  } else {
    stop(paste("The ",futFile, " file does not exist. Create this file with 'aqmapClip()' before re-running"))
  }
  
  colnames(siteFut)<-c("SpeciesID","probability","geometry")
  
  ### Calculate species richness by cell (# of species)
  
  nowSR<-siteNow%>%
    filter(probability >= thresh)%>%
    mutate(geomID = as.character(geometry))%>%
    group_by(geomID)%>%
    summarise(sppRich=n())
  
  futSR<-siteFut%>%
    filter(probability >= thresh)%>%
    mutate(geomID=as.character(geometry))%>%
    group_by(geomID)%>%
    summarise(sppRich50=n())
  
  dSppRich<-nowSR%>%
    st_join(.,futSR)%>%
    mutate(srDiff=sppRich50-sppRich,
           lon=st_coordinates(.)[,1],
           lat=st_coordinates(.)[,2])
  
  ### plot changes in species richness
  
  ## get countries map
  
  ### crop base layers from world map
  
  if(is.na(corr)) {
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      group_by(group) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast("POLYGON")%>%
      st_crop(.,ploArea)
    
  } else {
    
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      filter(lon2>ploArea$xmin & lon2<ploArea$xmax,
             lat>ploArea$ymin & lat<ploArea$ymax)%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id="group")%>%
      st_make_valid()
  }
  
  ## get eez
  
  # mrEezName<-mr_geo_code(place = mrSearch, like = TRUE, fuzzy = FALSE)%>%filter(placeType=="EEZ")%>%select(accepted)%>%as.numeric()
  # 
  # siteEez<-mr_shp(key="MarineRegions:eez",filter=mrEezName)%>%
  #   st_as_sf(.,
  #            crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  # 
  # if(corr=="lonFx") {
  #   siteEez<-siteEez%>%
  #     st_as_sf(.,
  #              coords=c("long","lat"),
  #              crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
  #     st_cast(to="POINT")%>% # convert polygons to points
  #     mutate(lon=st_coordinates(geometry)[,"X"] ,
  #            lat=st_coordinates(geometry)[,"Y"])%>%
  #     mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
  #                           TRUE ~ lon))%>%
  #     st_drop_geometry()%>%
  #     sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(siteEez)[1])
  # }
  
  ## Plot layers
  
  dSRplot<-ggplot()+
    geom_tile(data=dSppRich,aes(x=lon,y=lat,fill=srDiff),width=0.5,height=0.5)+
    scale_fill_gradient2(low="navyblue",mid="white",high="firebrick",midpoint = 0,name="Change in\nnumber of species")+
    # geom_sf(siteEez,fill=NA,color="grey70")+
    geom_sf(data = geoShp,fill="grey",color="black")+
    geom_sf(data = mpa,fill=NA, color="deepskyblue4",linetype=2)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    labs(title="Changes in species richness (2019 to 2050)", subtitle = "Data source: Aquamaps",x="longitude",y="latitude")+
    theme_bw()
  
  return(dSRplot)
  
  
  }