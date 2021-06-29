###################################################################
## create raster with spectra from different numbers of species
##  - - random patterns - -
###################################################################

require(raster)
# load data
setwd("E:/KIT_Forschung/1_SVH_paper/Teja_SVH")
data <- read.csv("canopy_kit_plant_functional_gradient.csv")
setwd("E:/KIT_Forschung/1_SVH_paper/5_sims/sim_ras")


# get data for one time point
data_sm <- data[2:22,]
# extract reflectances
refl <- data_sm[,17:2117]
# reduce number of bands to only 1 our of 20 bands
subs <- seq(1,2100,20)
refl_sm2 <- refl[,subs]
# kick out bands affected by atmopshere and noise at the end of the SWIR range
refl_sm <- refl_sm2[,-c(71:79, 101:105)]
# check whether spectra look ok
plot(seq(1,91,1), refl_sm[1,])
# we have spectra with 91 bands left


# define function to create raster with
# each pixel representing a spectrum of 
# a species

createSimRaster <- function(dim_x, dimy_y, species){
  
  # create empty array (100 x 100 pixels and 2100 bands - this matches the
  # spectra available in our base data)

  empty <- array(NA, dim=c(dim_x,dim_y,91))
  
  
  # create 10000 random numbers to assign each "pixel" of array to 
  # to one of the 23 species

  # create required combinations of pixel positions
  x <- seq(1, dim_x, 1)
  y <- seq(1, dim_y, 1)
  d1 <- expand.grid(x = x, y = y)
  
  # now iterate through pixel locations and randomly insert spectra
  for (i in 1:nrow(d1)){
    # get spectra to insert
    spec <- as.numeric(refl_sm[species[i],])
    # insert spectra
    empty[d1[i,1], d1[i,2], ] <- spec
    #print(i)
    
  }
  
  # convert array to raster-brick
  full <- brick(empty)
  # plot raster
  plotRGB(full, r=45, g=30, b=15, stretch="hist")
  # return raster stack
  full
}


## now create some rasters for different amounts of species

# set size of raster
dim_x = 100
dim_y = 100


set.seed(21)

nrit = 100

# now create 10 rasters of 23 species

spec22_res <- list()

set.seed(22)

for (i22 in 1:nrit){
  
  # prepare random species vector with desired number of species
  # in this case we consider all 23 species of time point 1
  species22 <- sample(seq(1,22,1), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species22)
  spec22_res[[i22]] <- img
  print(i22)
  
}

save(spec22_res, file = "spec22_100ras.RData")

# now create 10 rasters of 4 species

spec4_res <- list()

set.seed(22)

for (i4 in 1:nrit){
  
  species4 <- sample(sample(seq(1,23,1),4, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species4)
  spec4_res[[i4]] <- img
  print(i4)
  
}

save(spec4_res, file = "spec4_100ras.RData")


# now create 10 rasters of 8 species

spec8_res <- list()

set.seed(22)

for (i8 in 1:nrit){
  
  species8 <- sample(sample(seq(1,23,1),8, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species8)
  spec8_res[[i8]] <- img
  print(i8)
  
}

save(spec8_res, file = "spec8_100ras.RData")

# now create 10 rasters of 12 species

spec12_res <- list()

set.seed(22)

for (i12 in 1:nrit){
  
  species12 <- sample(sample(seq(1,23,1),12, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species12)
  spec12_res[[i12]] <- img
  print(i12)
  
}

save(spec12_res, file = "spec12_100ras.RData")


# now create 10 rasters of 16 species

spec16_res <- list()

set.seed(22)

for (i16 in 1:nrit){
  
  species16 <- sample(sample(seq(1,23,1),16, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species16)
  spec16_res[[i16]] <- img
  print(i16)
  
}

save(spec16_res, file = "spec16_100ras.RData")


# now create 10 rasters of 20 species

spec20_res <- list()

set.seed(22)

for (i20 in 1:nrit){
  
  species20 <- sample(sample(seq(1,23,1),20, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species20)
  spec20_res[[i20]] <- img
  print(i20)
  
}

save(spec20_res, file = "spec20_100ras.RData")






