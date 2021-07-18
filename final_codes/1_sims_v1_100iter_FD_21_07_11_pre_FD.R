###################################################################
## create raster with spectra from different numbers of species
##  - - random patterns - -
###################################################################


require(raster)
require(data.table)
require(matrixStats)
require(FD)
# load data
setwd("E:/KIT_Forschung/1_SVH_paper/127_IAVS_special_issue_spectral_diversity/Teja_SVH")
data <- read.csv("canopy_kit_plant_functional_gradient.csv")
setwd("E:/KIT_Forschung/1_SVH_paper/5_sims_new/sim_ras_07_12")


# get data for one time point
data_sm <- data[2:22,]
colnames(data_sm[,1:16])
# extract reflectances
refl <- data_sm[,17:2117]
water_abs = c(1400:1500,1800:2000,2400:2500)
bands = seq(400,2500, by = 1)
bands = bands[-which(bands %in% water_abs==TRUE)]
dummy = bands-399

refl2 <- refl[,dummy]

# reduce number of bands to only 1 out of 20 bands
subs <- seq(1,1698,20)
refl_sm <- refl2[,subs]
# check whether spectra look ok
plot(seq(1,85,1), refl_sm[1,])
# we have spectra with 85 bands left


# define function to create raster with
# each pixel representing a spectrum of 
# a species

createSimRaster <- function(dim_x, dimy_y, species){
  
  # create empty array (100 x 100 pixels and 2100 bands - this matches the
  # spectra available in our base data)

  empty <- array(NA, dim=c(dim_x,dim_y,85))
  
  
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


nrit = 100

# now create 10 rasters of 23 species

spec21_res <- list()
spec21_spec <- list()
spec21_cov <- list()

set.seed(8)

for (i21 in 1:nrit){
  
  # prepare random species vector with desired number of species
  # in this case we consider all 23 species of time point 1
  species21 <- sample(seq(1,21,1), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species21)
  spec21_res[[i21]] <- img
  print(i21)
  
  # calculate and save cover values and species of each image
  cov <- as.data.frame(table(species21))
  cov[,1] <- as.numeric(as.character(cov[,1]))
  cov[,2] <- as.numeric(cov[,2])
  
  #colnames(cov)
  alls <- colSums(cov)[2]
  
  cov_out <- as.data.frame(t(cov$Freq/alls))
  colnames(cov_out) <- data_sm$species[cov$species21]
  
  spec21_cov[[i21]] <- cov_out
  spec21_spec[[i21]]<- species21

}

save(spec21_res, file = "spec21_100ras.RData")
save(spec21_cov, file = "spec21_100cov.RData")
save(spec21_spec, file = "spec21_100spec.RData")

# now create 10 rasters of 4 species

spec4_res <- list()
spec4_spec <- list()
spec4_cov <- list()

set.seed(8)

for (i4 in 1:nrit){
  
  species4 <- sample(sample(seq(1,21,1),4, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species4)
  spec4_res[[i4]] <- img
  print(i4)
  
  # calculate FD of each image
  # get cover values of each species
  cov4 <- as.data.frame(table(species4))
  cov4[,1] <- as.numeric(as.character(cov4[,1]))
  cov4[,2] <- as.numeric(cov4[,2])
  
  colnames(cov)
  alls4 <- colSums(cov4)[2]
  
  cov_out4 <- as.data.frame(t(cov4$Freq/alls4))
  colnames(cov_out4) <- data_sm$species[cov4$species4]

  spec4_cov[[i4]] <- cov_out4
  spec4_spec[[i4]]<- species4
  rm(img)
  
}

save(spec4_res, file = "spec4_100ras.RData")
save(spec4_cov, file = "spec4_100cov.RData")
save(spec4_spec, file = "spec4_100spec.RData")

# now create 10 rasters of 8 species

spec8_res <- list()
spec8_spec <- list()
spec8_cov <- list()

set.seed(8)

for (i8 in 1:nrit){
  
  species8 <- sample(sample(seq(1,21,1),8, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species8)
  spec8_res[[i8]] <- img
  print(i8)
  
  # calculate FD of each image
  # get cover values of each species
  cov8 <- as.data.frame(table(species8))
  cov8[,1] <- as.numeric(as.character(cov8[,1]))
  cov8[,2] <- as.numeric(cov8[,2])
  
  colnames(cov)
  alls8 <- colSums(cov8)[2]
  
  cov_out8 <- as.data.frame(t(cov8$Freq/alls8))
  colnames(cov_out8) <- data_sm$species[cov8$species8]
  
  spec8_cov[[i8]] <- cov_out8
  spec8_spec[[i8]]<- species8
  rm(img)
  
}

save(spec8_res, file = "spec8_100ras.RData")
save(spec8_cov, file = "spec8_100cov.RData")
save(spec8_spec, file = "spec8_100spec.RData")

# now create 10 rasters of 12 species

spec12_res <- list()
spec12_spec <- list()
spec12_cov <- list()

set.seed(22)

for (i12 in 1:nrit){
  
  species12 <- sample(sample(seq(1,21,1),12, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species12)
  spec12_res[[i12]] <- img
  print(i12)
  
  # calculate FD of each image
  # get cover values of each species
  cov12 <- as.data.frame(table(species12))
  cov12[,1] <- as.numeric(as.character(cov12[,1]))
  cov12[,2] <- as.numeric(cov12[,2])
  
  colnames(cov)
  alls12 <- colSums(cov12)[2]
  
  cov_out12 <- as.data.frame(t(cov12$Freq/alls12))
  colnames(cov_out12) <- data_sm$species[cov12$species12]
  
  spec12_cov[[i12]] <- cov_out12
  spec12_spec[[i12]]<- species12
  rm(img)
  
}

save(spec12_res, file = "spec12_100ras.RData")
save(spec12_cov, file = "spec12_100cov.RData")
save(spec12_spec, file = "spec12_100spec.RData")

# now create 10 rasters of 16 species

spec16_res <- list()
spec16_spec <- list()
spec16_cov <- list()

set.seed(22)

for (i16 in 1:nrit){
  
  species16 <- sample(sample(seq(1,21,1),16, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species16)
  spec16_res[[i16]] <- img
  print(i16)
  
  # calculate FD of each image
  # get cover values of each species
  cov16 <- as.data.frame(table(species16))
  cov16[,1] <- as.numeric(as.character(cov16[,1]))
  cov16[,2] <- as.numeric(cov16[,2])
  
  colnames(cov)
  alls16 <- colSums(cov16)[2]
  
  cov_out16 <- as.data.frame(t(cov16$Freq/alls16))
  colnames(cov_out16) <- data_sm$species[cov16$species16]
  
  spec16_cov[[i16]] <- cov_out16
  spec16_spec[[i16]]<- species16
  rm(img)
  
}

save(spec16_res, file = "spec16_100ras.RData")
save(spec16_cov, file = "spec16_100cov.RData")
save(spec16_spec, file = "spec16_100spec.RData")

# now create 10 rasters of 20 species

spec20_res <- list()
spec20_spec <- list()
spec20_cov <- list()

set.seed(22)

for (i20 in 1:nrit){
  
  species20 <- sample(sample(seq(1,21,1),20, replace=F), dim_x*dim_y, replace=T)
  img <- createSimRaster(dim_x, dim_y, species20)
  spec20_res[[i20]] <- img
  print(i20)
  
  # calculate FD of each image
  # get cover values of each species
  cov20 <- as.data.frame(table(species20))
  cov20[,1] <- as.numeric(as.character(cov20[,1]))
  cov20[,2] <- as.numeric(cov20[,2])
  
  colnames(cov)
  alls20 <- colSums(cov20)[2]
  
  cov_out20 <- as.data.frame(t(cov20$Freq/alls20))
  colnames(cov_out20) <- data_sm$species[cov20$species20]
  
  spec20_cov[[i20]] <- cov_out20
  spec20_spec[[i20]]<- species20
  rm(img)
  
}

save(spec20_res, file = "spec20_100ras.RData")
save(spec20_cov, file = "spec20_100cov.RData")
save(spec20_spec, file = "spec20_100spec.RData")





