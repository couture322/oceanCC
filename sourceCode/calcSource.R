### mapping functions for site assessments
### written by Jessica Couture, Research Scientist
### Marine SPARC, Conservation International


##### Libraries needed
library(tidyverse)
library(mregions)

##################################################################
########################### Extent Map ###########################
##################################################################
# mpaPoly = sf object of MPA shape file
# siteName = string representing the site to be mapped
# mrSearch = string to search marine regions DB for EEZ shp file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# zoomLat = latitudinal buffer for the entire extent map in decimal degrees (size depends on context extent map scale)
# zoomLon = longitudinal buffer for the entire extent map in decimal degrees (size depends on context extent map scale)

# output ggplot map object with MPA shapefile, plot area and surrounding countries plotted
# !!! NOT coded for sites that cross the 180 meridian (like the Lau Seascape, see LauCC for this code)

extMapFunc<-function(mpaPoly,siteName,mrSearch,buff,zoomLat,zoomLon,corr=c(NA,"lonFx")){
  
  mpa<-mpaPoly%>%
    st_set_crs("EPSG 4326")#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  if(corr=="lonFx") {
    mpa<-mpa%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaPoly)[1])
    # st_crs(crs=4326)#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  }
  
  siteBB<-st_bbox(mpa)
  
  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  } else {if(corr=="lonFx") {
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
  }
  
  bbDf<-data.frame(lon=c(as.numeric(ploArea$xmin),
                         as.numeric(ploArea$xmax),
                         as.numeric(ploArea$xmax),
                         as.numeric(ploArea$xmin)),
                   lat=c(as.numeric(ploArea$ymin),
                         as.numeric(ploArea$ymin),
                         as.numeric(ploArea$ymax),
                         as.numeric(ploArea$ymax)))%>%
    st_as_sf(.,
             coords=c("lon","lat"))%>%
    mutate(id=1)%>%
    group_by(id)%>%
    summarise()%>%
    st_convex_hull()
  # st_cast("POLYGON")
  
  
  
  if(is.na(corr)) {
    totArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-zoomLon),
                       xmax=as.numeric(siteBB$xmax+zoomLon),
                       ymin=as.numeric(siteBB$ymin-zoomLat),
                       ymax=as.numeric(siteBB$ymax+zoomLat)),crs = st_crs(mpa))
  } else {if(corr=="lonFx") {
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
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id="group")%>%
      st_make_valid()
    
  }
  
  
  mrEezName<-mr_geo_code(place = mrSearch, like = TRUE, fuzzy = FALSE)%>%filter(placeType=="EEZ")%>%select(accepted)%>%as.numeric()
  
  siteEez<-mr_shp(key="MarineRegions:eez",filter=mrEezName)%>%
    st_as_sf(.,
             crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  if(corr=="lonFx") {
    siteEez<-siteEez%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(siteEez)[1])
  }
  
  extMap<-ggplot()+
    geom_sf(data=geoShp,fill="grey")+
    geom_sf(data=siteEez,fill="grey33",color=NA,alpha=0.3)+
    geom_sf(data=mpa,fill="deepskyblue3",alpha=0.3,color="deepskyblue4",linetype=2)+
    geom_sf(data=bbDf,fill=NA,color="darkorange3",linetype=2)+
    scale_x_continuous(expand = c(0,0),limits = c(totArea$xmin,totArea$xmax))+
    scale_y_continuous(expand = c(0,0),limits = c(totArea$ymin,totArea$ymax))+
    theme_bw()+
    labs(subtitle = "Extent map",title=siteName,x="longitude",y="latitude")
  
  return(extMap)
  
}



##################################################################
##################### Temperatures at depth ######################
########################### Static ###############################

# mpaPoly = sf object of MPA shape file
# buff = distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
# targDepth = one of "epi" for epipelagic, "meso" for mesoplagic, "bathy" for bathypelagic to indicate which dataset to pull from the drive
# YEAR = year in the future for analysis from 2020 - 2100, default is 2050
# corr = indicates whether longitude (lonFx) needs to be corrected for crossing the 180 meridian, if no use NA

# output ggplot map object with MPA shapefile, plot area and surrounding countries plotted
# !!! NOT coded for sites that cross the 180 meridian (like the Lau Seascape, see LauCC for this code)


temp3dPlots<-function(mpaPoly,buff,targDepth=c("epi","meso","bathy"),YEAR = 2050,corr=c(NA, "lonFx")){
  
  
  ### set bounding box for target area
  mpa<-mpaPoly%>%
    st_set_crs("EPSG 4326")
  
  if(corr=="lonFx") {
    mpa<-mpa%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaPoly)[1])
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
  
  geoShp<-map_data("world")%>%
    st_as_sf(.,
             coords=c("long","lat"),
             crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
    group_by(group) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON")%>%
    st_crop(.,ploArea)
  
  # print(head(geoShp))
  
  if(corr=="lonFx") {
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
  
  if(corr=="lonFx") {
    
    dTempSpDf<-dTempRast%>%
      as.data.frame(.,xy=TRUE)%>%
      mutate(lon=case_when(x<0 ~ 180+(180+x),
                           TRUE ~ x))%>%
      rename(lat=y)%>%
      filter(lon>ploArea["xmin"],
             lon<ploArea["xmax"],
             lat>ploArea["ymin"],
             lat<ploArea["ymax"])
    
  } else {
    
    dTempSpDf<-crop(dTempRast,ploArea)%>% # crop to intended area
      as.data.frame(.,xy=TRUE)%>% ## convert to df
      rename(lon=x,
             lat=y)
    
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
