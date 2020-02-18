
# packages
library(raster)
library(sf)
library(viridis)
library(gridExtra)
library(ggplot2)
library(dplyr)

rasterOptions(maxmemory = 40e9, chunksize = 4e9, tmpdir="F:/Temp/") # Define this based on your computer settings.

### Download data

# Set working directory to source file location first

# Download the winter flux synthesis data to the data folder from here: https://doi.org/10.3334/ORNLDAAC/1692
# (or any other synthesis data set that you want to visualize here)
d <- read.csv("../data/nongrowing_season_CO2_flux.csv")

# Convert to spatial points data frame in WGS 84 projection
sp <- d
class(sp)
coordinates(sp) <- ~longitude + latitude
proj4string(sp) = CRS("+init=epsg:4326") 


# Download the study domain shapefile. I use the permafrost distribution based on Brown et al. 2002. I downloaded the file from here https://nsidc.org/data/ggd318, and transformed it to WGS 84 projection in ArcMap because R creates artificial boundaries to the file for some reason.
# (or any other study domain that you want to use here)
perma <- shapefile("../data/perma_dissolved_wgs.shp")
# Crop to the northern regions only
perma <- crop(perma, extent(-180, 180, 48, 90))

# Download a world map for visualizations
data("wrld_simpl", package = "maptools")        
wm <- crop(wrld_simpl, extent(-180, 180, 48, 90))   

# Let's use the data sets that are available via raster package to characterize the environmental space
?getData

# Download WordClim rasters. First, current climate
climate <- getData("worldclim",var="bio",res=10)
# Select annual precipitation and temperature (see here: http://www.worldclim.org/bioclim)
climate <- climate[[c(1,12)]]
names(climate) <- c("Temp","Prec")

# Future climate
fut_climate <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70)
fut_climate <- fut_climate[[c(1,12)]]
names(fut_climate) <- c("Fut_Temp","Fut_Prec")

# Create a raster stack 
climate_all <- stack(climate, fut_climate)

# Crop the stack to the study domain
r2 <- crop(climate_all, extent(perma), "rr2.tif", overwrite=TRUE)
plot(r2)
r3 <- mask(r2, perma, "rr3.tif", overwrite=TRUE)
plot(r3)
names(r3) <- names(climate_all)

# Extract climate data to winter flux sites
clim <- as.data.frame(extract(climate_all, sp))


### 1. Plot the observations on a map

plot(perma)
plot(sp, add=T, p=1, col="red")


### 2. Plot the density of observations on a map

# We first create an empty raster that is used to show the number of observations per raster cell
# Add a raster from the stack
r <- r3[[1]] 
# resample to ~100 km so that the number of observations is shown per 100 x 100 pixel
a = raster(extent(r)); res(a)=c(1, 1) 
    
# The add the number of observations to each pixel
a[] <- 0
tab <- table(cellFromXY(a, sp))
a[as.numeric(names(tab))] <- tab
a[a==0] <- NA

# Set the breaks for the plot 
# breaks for 1 x 1 deg
start <- minValue(a)
breaks <- c(start, 15,30, 45, 300)
labels <- c(start, 15,30, 45, 300)
limits <- c(start, 300)

# Plot
rasterVis::gplot(a) +
  geom_polygon(data = wm, aes(x = long, y = lat, group = group), fill="white", color="black", size=0.001) + 
  geom_polygon(data = perma, aes(x = long, y = lat, group = group), fill = "grey", colour = "black", alpha = 0.9) +
  coord_map("stereographic", orientation=c(90, 0, 0), ylim=50)  + theme_void() +
  geom_tile(aes(fill=value), alpha=0.8) + scale_fill_viridis(option="plasma", direction=-1, na.value="transparent", limits=limits, breaks=breaks, labels=labels)+ theme(legend.position="none") + ggtitle("Density of observations") + labs(fill="") + theme(plot.title=element_text(size=18, face="bold"))



### 3. Observations across countries (or some other domains)
countries <- extract(wm, sp)
countries <- countries[!is.na(countries$NAME), ]

# Relative number
countries_count <- countries %>% group_by(NAME) %>% summarise(n=n())
p1 <- ggplot(countries_count, aes(NAME, n)) + geom_bar(stat="identity") + ggtitle("Number of observations in each domain")


# The spatial extent of the same domains
wm$area_sqkm  <- area(wm) / 1000000
countries_count <- merge(countries_count, wm, by="NAME")
p2 <- ggplot(countries_count, aes(NAME, area_sqkm)) + geom_bar(stat="identity") + ggtitle("The area of each domain")

# Plot
grid.arrange(p1, p2)



### 4. Scatterplots of climatic gradients


# A basic scatterplot

climate_df <- as.data.frame(r3)
names(climate_df) <- names(clim)

# Current climate
p <- ggplot(climate_df, aes(Temp, Prec)) + geom_point()
p + geom_point(data=clim, col="red")

# Future climate
p <- ggplot(climate_df, aes(Fut_Temp, Fut_Prec)) + geom_point()
p + geom_point(data=clim, col="red")



# Heatmap of 2d bin counts across climatic gradients

# set the limits
minx <- min(climate_df$Temp, na.rm=TRUE)
maxx <- max(climate_df$Temp, na.rm=TRUE)

miny <- min(climate_df$Prec, na.rm=TRUE)
maxy <- max(climate_df$Prec, na.rm=TRUE)

p1 <- ggplot(climate_df, aes(Temp, Prec) ) +
  geom_bin2d() +
 ggtitle("Whole region") + xlim(minx, maxx) + ylim(miny, maxy)

p2 <- ggplot(clim, aes(Temp, Prec) ) +
  geom_bin2d() +
  ggtitle("Sites") + xlim(minx, maxx) + ylim(miny, maxy)


grid.arrange(p1, p2)



# Smoother density estimates (i.e. a smoothed version of the histogram) for climatic gradients

climate_df$class <- "Whole"
clim$class <- "Sites"

climate_merged <- rbind(climate_df, clim)

p1 <- ggplot(climate_merged) + geom_density(aes(Temp, col=factor(class)))

p2 <- ggplot(climate_merged) + geom_density(aes(Prec, col=factor(class)))

grid.arrange(p1, p2)





### 4. Where are the gaps spatially?

# Based on min and max
minx <- min(clim$Temp, na.rm=TRUE)
maxx <- max(clim$Temp, na.rm=TRUE)

miny <- min(clim$Prec, na.rm=TRUE)
maxy <- max(clim$Prec, na.rm=TRUE)

# Separately
plot(r3[[1]]<minx | r3[[1]]>maxx)
plot(r3[[2]]<miny | r3[[2]]>maxy)

# Together
plot(r3[[2]]<miny | r3[[2]]>maxy | r3[[1]]<minx | r3[[1]]>maxx)


# Based on quantiles
minx <- quantile(clim$Temp, c(0.05), na.rm=TRUE)
maxx <- quantile(clim$Temp, c(0.95), na.rm=TRUE)

miny <- quantile(clim$Prec, c(0.05), na.rm=TRUE)
maxy <- quantile(clim$Prec, c(0.95), na.rm=TRUE)

plot(r3[[1]]<minx | r3[[1]]>maxx)
plot(r3[[2]]<miny | r3[[2]]>maxy)
plot(r3[[2]]<miny | r3[[2]]>maxy | r3[[1]]<minx | r3[[1]]>maxx)



### 5. A more analytical approach to explore gaps across space (see Virkkala et al. 2019 ERL)


# Create new "pseudo" sites that do not have any measurements
new_sites <- sampleRandom(r3, size=1100, sp=T) 
plot(new_sites)

# Extract climate data to these sites
climate_df_newsites <- data.frame(extract(r3, new_sites))
names(climate_df_newsites) <- names(clim)[1:4]
# Create a site identificator (no site)
climate_df_newsites$site <- 0

# Create a site identificator for the existing site data
clim$site <- 1

# Merge the site and pseudo site data frames
modeldata <- rbind(subset(clim, select=-c(class)), climate_df_newsites)


# Create the model
library(gbm)

# Parameters - these can be tuned depending on how general versus specific you want the model to be. See https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2656.2008.01390.x
Interaction_Depth <- 5
Numb_Trees <- 100 
Min_Obs <- 20 

# Train the model
m1 <- gbm(site ~ Temp + Prec, data=modeldata, distribution="bernoulli", n.trees=Numb_Trees, n.minobsinnode = Min_Obs, interaction.depth = Interaction_Depth)
# Estimate the optimal number of boosting iterations
trees_opt <- gbm.perf(m1)
# Predict with the model
pred <- predict(r3, m1,  type="response",  progress="text", n.trees=trees_opt)

# Plot the prediction. High values = high probability that there is a site in those environmental conditions
plot(pred)







