####Housekeeping####
rm(list=ls())
setwd('/Users/aaron/OneDrive - Imperial College London/University/Project/data/') #CHANGE THIS

####Loading Libraries####
#library("maps")
#library("mapdata")
library("ggplot2")
library("ggmap")
library(sf)
library(rgdal)
library(broom)
library(readxl)

meta <- read_excel("FieldDataSheet_environment_map.xlsx")
names(meta) <- tolower(names(meta))
# create the land-use variable from the transect variable
meta$landuse <- substring(meta$transect, 1, 2)

# correct any inconsistencies in names convention of the landuse
meta <- meta %>% mutate(landuse = str_replace(landuse, "Wf", "WF")) %>% 
  mutate(landuse = str_replace(landuse, "Df", "DF")) %>% 
  mutate(landuse = str_replace(landuse, "Cf", "CF")) %>% 
  mutate(landuse = str_replace(landuse, "Ug", "UG"))
register_google(key="AIzaSyDfB1YygYEdyBi08YBZkVupiP3obA-qZmQ")

# llmeans1 <- meta %>% filter(landuse == "CF")
# llmeans1 <- sapply(llmeans1[18:17], mean)

llmeans2 <- meta %>% filter(landuse != "CF")
llmeans2 <- sapply(llmeans2[18:17], mean)

# create a meta file without CF in it 
meta2 <- meta %>% group_by(transect) %>% 
                  filter(landuse != "CF") %>% 
                  summarise(long. = mean(long.), lat. = mean(lat.))

# create a subset only of CF transects
meta4 <- meta %>% group_by(transect) %>% 
  filter(landuse == "CF") %>% 
  summarise(long. = mean(long.), lat. = mean(lat.))

# merge the the meta file and data_final from final R script 
meta3 <- merge(meta2, data_final)
llmeans <- llmeans2

# get sample map
sample_map <- get_map(location=llmeans,
                     maptype="satellite",
                     source="google",
                     zoom=15) #Change this bit to get the right zoom

sat_image1 <- ggmap(sample_map) +
geom_point(data=meta3, aes(x=long., y=lat., color=landuse, shape = landuse), size=2.8)+ #Change this line for if you want lines between the starts and ends I guess
  #You could also colour-code by transect or transect type (I put some nice columns in the excel sheet for this)
  labs(y="Latitude (ºN)", x="Longitude (ºE)") #+
  #geom_polygon(data = shp_tran)

# con forest
llmeans3 <- meta %>% filter(landuse == "CF")
llmeans3 <- sapply(llmeans3[18:17], mean)

# create a meta file
meta4 <- meta %>% group_by(transect) %>% 
  filter(landuse == "CF") %>% 
  summarise(long. = mean(long.), lat. = mean(lat.))

# merge the meta 4 file and data_final (has to be loaded from other script)
meta5 <- merge(meta4, data_final)

llmeans <- llmeans3

# get the sample map
sample_map <- get_map(location=llmeans,
                      maptype="satellite",
                      source="google",
                      zoom=15) #Change this bit to get the right zoom

sat_image2 <- ggmap(sample_map) +
  geom_point(data=meta5, aes(x=long., y=lat., color=landuse, shape = landuse), size=2.8)+ #Change this line for if you want lines between the starts and ends I guess
  #You could also colour-code by transect or transect type (I put some nice columns in the excel sheet for this)
  labs(y="Latitude (ºN)", x="Longitude (ºE") #+
#geom_polygon(data = shp_tran)


library(cowplot)
library(grDevices)
# plot the maps in one figure
rud_both <- ggdraw() +
  draw_plot(sat_image1, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(sat_image2, x = 0.5, y = 0, width = 0.5, height = 1) 

rud_both
