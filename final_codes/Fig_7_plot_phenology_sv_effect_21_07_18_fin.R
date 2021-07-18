require(raster)
require(RStoolbox)

###########################################
##### function to create synthetic raster
###########################################

createSimRaster <- function(dim_x, dimy_y, species, refl_sm){
  
  # create empty array (100 x 100 pixels and 2100 bands - this matches the
  # spectra available in our base data)
  
  empty <- array(NA, dim=c(dim_x,dim_y,85))
  
  
  # create 10000 random numbers to assign each "pixel" of array to 
  # to one of the 23 species
  
  # create required combinations of pixel positions
  x <- seq(1, dim_x, 1)
  y <- seq(1, dim_y, 1)
  gr1 <- expand.grid(x = x, y = y)
  
  # now iterate through pixel locations and randomly insert spectra
  for (i99 in 1:nrow(gr1)){
    # get spectra to insert / select it based on species name of current pixel
    corspec <- names[species[i99]]
    spec <- as.numeric(refl_sm[refl_sm$species==corspec, 2:86])
    # insert spectra
    empty[gr1[i99,1], gr1[i99,2], ] <- spec
    #print(i)
    
  }
  
  # convert array to raster-brick
  full <- brick(empty)
  # plot raster
  plotRGB(full, r=45, g=30, b=15, stretch="hist")
  # return raster stack
  full
}

###########################################
##### function to calculate spectral variation
###########################################

get_specVarPCA <- function(img){
  
  #calculate PCA------
  pca_raster<-rasterPCA(img)
  class(pca_raster)
  
  image_1 <- pca_raster$map$PC1
  results <-mean(dist(getValues(image_1)[!is.na(getValues(image_1))]), na.rm=T)
  image_2<-pca_raster$map$PC2
  results2<-mean(dist(getValues(image_2)[!is.na(getValues(image_2))]), na.rm=T)
  val_1<- getValues(image_1)[!is.na(getValues(image_1))]
  val_2<- getValues(image_2)[!is.na(getValues(image_2))]
  dist_both<-sqrt((dist(val_1)^2)+(dist(val_2)^2))
  results3 <-mean(dist_both, na.rm=T)
  
  # return results
  print(results3)
  results3
  
}



###########################################
# create synthetic images across time
###########################################

#load base data
setwd("D:/1_SVH_paper_29_06_21/6_phenology_sv_effect")
data <- read.table("pheno_spec_var_data_ready.txt", sep="\t", header=T)

# get names from first time point (this is the reference for the already
# created species-compositions to create synthetic images used for Figure 3, 
# which will be re-used here)
names <- data$species[1:21]




########################################################################
## now create some rasters for 4 species images across time
########################################################################

#load species data for synthetic images with 4 species
setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_13/sim_ras/sim_rasters")
load("spec4_100spec.RData")

# set output directorx
setwd("D:/1_SVH_paper_29_06_21/6_phenology_sv_effect/rasters")

# set size of raster
dim_x = 100
dim_y = 100

# create empty list to store results
imgs_4spec <- list()

for (i41 in 1:8){
  
  # take first set of species (vector with length 10000)
  d1 <- spec4_spec[[i41]]
  # get unique values (these are 4)
  unspec <- as.numeric(unique(d1))
  # and corresponding species names
  curnam <- names[unspec]
  # reduce the full dataset to all rows with the four species
  data_cur <- data[data$species%in%curnam,]
  
  # check how many of the species are available per time-stamp
  compl1 <- table(data_cur$timestamp1.1.363.)
  # identify only time-stamps where all four species have been measured
  z <- as.numeric(which(compl1 == 4))
  
  # now create synthetic images
  imgs_temp <- list()
  
  if (z>0) {
  
    for (i4 in 1:length(z)){
      
      # select spectra of species of current time step
      data_calc <- data[data$timestamp1.1.363.==z[i4],18:2118]
      # remove noisy bands
      water_abs = c(1400:1500,1800:2000,2400:2500)
      bands = seq(400,2500, by = 1)
      bands = bands[-which(bands %in% water_abs==TRUE)]
      dummy = bands-399
      refl2 <- data_calc[,dummy]
      
      # subset spectral bands to reduce size of image
      subs <- seq(1,1698,20)
      refl_sm2 <- refl2[,subs]
      # add species-name to each spectra to allow to select
      # the correct spectra for each synthetic image and time step
      refl_sm <- cbind(data[data$timestamp1.1.363.==z[i4],2], refl_sm2)
      colnames(refl_sm)[1] <- "species"
      
      # call the function to create the synthetic image
      img <- createSimRaster(dim_x, dim_y, d1, refl_sm)
      # save image of current time step into list
      imgs_temp[[i4]] <- img
      print(i4)
        
    }
    
    # save temporal image stack into list
    imgs_4spec[[i41]] <- imgs_temp
    
  } else {next}


}

#########################################################
## and calculate corresponding spectral variation values
#########################################################

# create empty list to store results
spec_var_time_all4 <- list()

# loop through time-stacks of synthetic images
for (i42 in 1:length(imgs_4spec)){
  
    # take time stack of first synthetic image
    imglist <- imgs_4spec[[i42]]
    
    # create empty list to store results
    sv_time <- list()
    
    # loop through time-steps ot current stack
    for (i43 in 1:length(imglist)){
      
      # take time step
      img43 <- imglist[[i43]]  
      # calculate spectral variation
      specvar <- get_specVarPCA(img43)
      # store spectral variation to list
      sv_time[[i43]] <- specvar
      print(i43)
      
      }
      
    spec_var_time_all4[[i42]]<-sv_time
    
}

save(spec_var_time_all4, file="4_spec_sv_over_time.RData")

dev.off()



########################################################################
## now create some rasters for 8 species images across time
########################################################################

#load species data for synthetic images with 4 species
setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_13/sim_ras/sim_rasters")
load("spec8_100spec.RData")

setwd("D:/1_SVH_paper_29_06_21/6_phenology_sv_effect/rasters")

# set size of raster
dim_x = 100
dim_y = 100


imgs_8spec <- list()

for (i81 in 1:8){
  
  d1 <- spec8_spec[[i81]]
  unspec <- as.numeric(unique(d1))
  curnam <- names[unspec]
  
  data_cur <- data[data$species%in%curnam,]
  
  compl1 <- table(data_cur$timestamp1.1.363.)
  z <- as.numeric(which(compl1 == 8))
  
  imgs_temp <- list()
  
  if (z>0) {
    
    for (i8 in 1:length(z)){
      
      data_calc <- data[data$timestamp1.1.363.==z[i8],18:2118]
      water_abs = c(1400:1500,1800:2000,2400:2500)
      bands = seq(400,2500, by = 1)
      bands = bands[-which(bands %in% water_abs==TRUE)]
      dummy = bands-399
      
      refl2 <- data_calc[,dummy]
      
      
      subs <- seq(1,1698,20)
      refl_sm2 <- refl2[,subs]
      # kick out bands affected by atmopshere and noise at the end of the SWIR range
      refl_sm <- cbind(data[data$timestamp1.1.363.==z[i8],2], refl_sm2)
      colnames(refl_sm)[1] <- "species"
      
      img <- createSimRaster(dim_x, dim_y, d1, refl_sm)
      imgs_temp[[i8]] <- img
      print(i8)
      
    }
    
    imgs_8spec[[i81]] <- imgs_temp
    
  } else {next}
  
  
}

#########################################################
## and calculate corresponding spectral variation values
#########################################################

spec_var_time_all8 <- list()

for (i82 in 1:length(imgs_8spec)){
  
  imglist <- imgs_8spec[[i82]]
  
  sv_time <- list()
  
  for (i83 in 1:length(imglist)){
    
    img83 <- imglist[[i83]]  
    specvar <- get_specVarPCA(img83)
    sv_time[[i83]] <- specvar
    print(i83)
    
  }
  
  spec_var_time_all8[[i82]]<-sv_time
  
}


save(spec_var_time_all8, "8_spec_sv_over_time.RData")

########################################################################
## now create some rasters for 12 species images across time
########################################################################

#load species data for synthetic images with 4 species
setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_13/sim_ras/sim_rasters")
load("spec12_100spec.RData")

setwd("D:/1_SVH_paper_29_06_21/6_phenology_sv_effect/rasters")

# set size of raster
dim_x = 100
dim_y = 100


imgs_12spec <- list()

for (i121 in 1:8){
  
  d1 <- spec12_spec[[i121]]
  unspec <- as.numeric(unique(d1))
  curnam <- names[unspec]
  
  data_cur <- data[data$species%in%curnam,]
  
  compl1 <- table(data_cur$timestamp1.1.363.)
  z <- as.numeric(which(compl1 == 12))
  
  imgs_temp <- list()
  
  if (z>0) {
    
    for (i12 in 1:length(z)){
      
      data_calc <- data[data$timestamp1.1.363.==z[i12],18:2118]
      water_abs = c(1400:1500,1800:2000,2400:2500)
      bands = seq(400,2500, by = 1)
      bands = bands[-which(bands %in% water_abs==TRUE)]
      dummy = bands-399
      
      refl2 <- data_calc[,dummy]
      
      
      subs <- seq(1,1698,20)
      refl_sm2 <- refl2[,subs]
      # kick out bands affected by atmopshere and noise at the end of the SWIR range
      refl_sm <- cbind(data[data$timestamp1.1.363.==z[i12],2], refl_sm2)
      colnames(refl_sm)[1] <- "species"
      
      img <- createSimRaster(dim_x, dim_y, d1, refl_sm)
      imgs_temp[[i12]] <- img
      print(i12)
      
    }
    
    imgs_12spec[[i121]] <- imgs_temp
    
  } else {next}
  
  
}

#########################################################
## and calculate corresponding spectral variation values
#########################################################

spec_var_time_all12 <- list()

for (i122 in 1:length(imgs_12spec)){
  
  imglist <- imgs_12spec[[i122]]
  
  sv_time <- list()
  
  for (i123 in 1:length(imglist)){
    
    img123 <- imglist[[i123]]  
    specvar <- get_specVarPCA(img123)
    sv_time[[i123]] <- specvar
    print(i123)
    
  }
  
  spec_var_time_all12[[i122]]<-sv_time
  
}


save(spec_var_time_all12, "12_spec_sv_over_time.RData")


########################################################################
## now create some rasters for 16 species images across time
########################################################################

#load species data for synthetic images with 4 species
setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_13/sim_ras/sim_rasters")
load("spec16_100spec.RData")

setwd("D:/1_SVH_paper_29_06_21/6_phenology_sv_effect/rasters")

# set size of raster
dim_x = 100
dim_y = 100


imgs_16spec <- list()

for (i161 in 1:8){
  
  d1 <- spec16_spec[[i161]]
  unspec <- as.numeric(unique(d1))
  curnam <- names[unspec]
  
  data_cur <- data[data$species%in%curnam,]
  
  compl1 <- table(data_cur$timestamp1.1.363.)
  z <- as.numeric(which(compl1 == 16))
  
  imgs_temp <- list()
  
  if (z>0) {
    
    for (i16 in 1:length(z)){
      
      data_calc <- data[data$timestamp1.1.363.==z[i16],18:2118]
      water_abs = c(1400:1500,1800:2000,2400:2500)
      bands = seq(400,2500, by = 1)
      bands = bands[-which(bands %in% water_abs==TRUE)]
      dummy = bands-399
      
      refl2 <- data_calc[,dummy]
      
      
      subs <- seq(1,1698,20)
      refl_sm2 <- refl2[,subs]
      # kick out bands affected by atmopshere and noise at the end of the SWIR range
      refl_sm <- cbind(data[data$timestamp1.1.363.==z[i16],2], refl_sm2)
      colnames(refl_sm)[1] <- "species"
      
      img <- createSimRaster(dim_x, dim_y, d1, refl_sm)
      imgs_temp[[i16]] <- img
      print(i16)
      
    }
    
    imgs_16spec[[i161]] <- imgs_temp
    
  } else {next}
  
  
}

#########################################################
## and calculate corresponding spectral variation values
#########################################################

spec_var_time_all16 <- list()

for (i162 in 1:length(imgs_16spec)){
  
  imglist <- imgs_16spec[[i162]]
  
  sv_time <- list()
  
  for (i163 in 1:length(imglist)){
    
    img163 <- imglist[[i163]]  
    specvar <- get_specVarPCA(img163)
    sv_time[[i163]] <- specvar
    print(i163)
    
  }
  
  spec_var_time_all16[[i162]]<-sv_time
  
}


save(spec_var_time_all6, "16_spec_sv_over_time.RData")


########################################################################
## now create some rasters for 20 species images across time
########################################################################

#load species data for synthetic images with 4 species
setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_13/sim_ras/sim_rasters")
load("spec20_100spec.RData")

setwd("D:/1_SVH_paper_29_06_21/6_phenology_sv_effect/rasters")

# set size of raster
dim_x = 100
dim_y = 100


imgs_20spec <- list()

for (i201 in 1:8){
  
  d1 <- spec20_spec[[i201]]
  unspec <- as.numeric(unique(d1))
  curnam <- names[unspec]
  
  data_cur <- data[data$species%in%curnam,]
  
  compl1 <- table(data_cur$timestamp1.1.363.)
  z <- as.numeric(which(compl1 == 20))
  
  imgs_temp <- list()
  
  if (z>0) {
    
    for (i20 in 1:length(z)){
      
      data_calc <- data[data$timestamp1.1.363.==z[i20],18:2118]
      water_abs = c(1400:1500,1800:2000,2400:2500)
      bands = seq(400,2500, by = 1)
      bands = bands[-which(bands %in% water_abs==TRUE)]
      dummy = bands-399
      
      refl2 <- data_calc[,dummy]
      
      
      subs <- seq(1,2098,20)
      refl_sm2 <- refl2[,subs]
      # kick out bands affected by atmopshere and noise at the end of the SWIR range
      refl_sm <- cbind(data[data$timestamp1.1.363.==z[i20],2], refl_sm2)
      colnames(refl_sm)[1] <- "species"
      
      img <- createSimRaster(dim_x, dim_y, d1, refl_sm)
      imgs_temp[[i20]] <- img
      print(i20)
      
    }
    
    imgs_20spec[[i201]] <- imgs_temp
    
  } else {next}
  
  
}

#########################################################
## and calculate corresponding spectral variation values
#########################################################

spec_var_time_all20 <- list()

for (i202 in 1:length(imgs_20spec)){
  
  imglist <- imgs_20spec[[i202]]
  
  sv_time <- list()
  
  for (i203 in 1:length(imglist)){
    
    img203 <- imglist[[i203]]  
    specvar <- get_specVarPCA(img203)
    sv_time[[i203]] <- specvar
    print(i203)
    
  }
  
  spec_var_time_all20[[i202]]<-sv_time
  
}

save(spec_var_time_all20, "20_spec_sv_over_time.RData")





############################################################################
## now plot spectral variation over time for each of the 8 x 4 datasets
############################################################################

# load colorBrewer to create colors
library(RColorBrewer)
colgrad2 = brewer.pal(8, "Set1")


#par(mfrow=c(1,4))
dev.off()


png(filename = "spectral_var_phen.png", width=2400, height=1400, res=300)

layout(matrix(c(1,1,1,2,2,3,4), nrow = 1, ncol = 7, byrow = TRUE))
par(mar=c(5,6,4,1))
# create empty plot as basic
plot(seq(0,14,1),cex.axis=1.3, seq(0,1.4,0.1), main="4 species", type="l", xlim=c(1,16), ylim=c(0,1.2), col="white", ylab="Spectral variation", xlab="Timestep", cex.lab=2)

# loop through the 8 synthetic images
for (u in 1:8){
  
  # take first set of species (vector with length 10000)
  d1 <- spec4_spec[[u]]
  # get unique values (these are 4)
  unspec <- as.numeric(unique(d1))
  # and corresponding species names
  curnam <- names[unspec]
  # reduce the full dataset to all rows with the four species
  data_cur <- data[data$species%in%curnam,]
  
  # check how many of the species are available per time-stamp
  compl1 <- table(data_cur$timestamp1.1.363.)
  # identify only time-stamps where all four species have been measured
  z <- as.numeric(which(compl1 == 4))
  par(new=T)
  # plot time-stamps against spectral variation
  plot(z, spec_var_time_all4[[u]], type="l", xlim=c(1,16), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
  par(new=T)
  plot(z, spec_var_time_all4[[u]], pch=0, xlim=c(1,16), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
}


par(mar=c(5,1,4,1))
plot(seq(0,14,1), seq(0,1.4,0.1), main="8 species", type="l", axes=F, xlim=c(1,12),ylab="", xlab="", ylim=c(0,1.2), col="white")
axis(1, cex.axis=1.3)
box()

for (u in 1:8){
  
  # take first set of species (vector with length 10000)
  d1 <- spec8_spec[[u]]
  # get unique values (these are 4)
  unspec <- as.numeric(unique(d1))
  # and corresponding species names
  curnam <- names[unspec]
  # reduce the full dataset to all rows with the four species
  data_cur <- data[data$species%in%curnam,]
  
  # check how many of the species are available per time-stamp
  compl1 <- table(data_cur$timestamp1.1.363.)
  # identify only time-stamps where all four species have been measured
  z <- as.numeric(which(compl1 == 8))
  par(new=T)
  # plot time-stamps against spectral variation
  plot(z, spec_var_time_all8[[u]], type="l", xlim=c(1,12), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
  par(new=T)
  plot(z, spec_var_time_all8[[u]], pch=0, xlim=c(1,12), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
}

par(mar=c(5,1,4,1))
plot(seq(0,14,1), seq(0,1.4,0.1), main="12 species", axes=F, type="l", xlim=c(1,3), ylab="", xlab="",ylim=c(0,1.2), col="white")
axis(1, at=1:3, labels=c(1,2,3), cex.axis=1.3)
box()

for (u in 1:8){
  
  # take first set of species (vector with length 10000)
  d1 <- spec12_spec[[u]]
  # get unique values (these are 4)
  unspec <- as.numeric(unique(d1))
  # and corresponding species names
  curnam <- names[unspec]
  # reduce the full dataset to all rows with the four species
  data_cur <- data[data$species%in%curnam,]
  
  # check how many of the species are available per time-stamp
  compl1 <- table(data_cur$timestamp1.1.363.)
  # identify only time-stamps where all four species have been measured
  z <- as.numeric(which(compl1 == 12))
  par(new=T)
  # plot time-stamps against spectral variation
  plot(z, spec_var_time_all12[[u]], type="l", xlim=c(1,3), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
  par(new=T)
  plot(z, spec_var_time_all12[[u]], pch=0, xlim=c(1,3), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
}

par(mar=c(5,1,4,1))
plot(seq(0,14,1), seq(0,1.4,0.1), main="16 species", type="l", xlim=c(1,3), ylab="", xlab="", axes=F, ylim=c(0,1.2), col="white")
axis(1, at=1:3, labels=c(1,2,3), cex.axis=1.3)
box()

for (u in 1:8){
  
  # take first set of species (vector with length 10000)
  d1 <- spec16_spec[[u]]
  # get unique values (these are 4)
  unspec <- as.numeric(unique(d1))
  # and corresponding species names
  curnam <- names[unspec]
  # reduce the full dataset to all rows with the four species
  data_cur <- data[data$species%in%curnam,]
  
  # check how many of the species are available per time-stamp
  compl1 <- table(data_cur$timestamp1.1.363.)
  # identify only time-stamps where all four species have been measured
  z <- as.numeric(which(compl1 == 16))
  par(new=T)
  # plot time-stamps against spectral variation
  plot(z, spec_var_time_all16[[u]], type="l", xlim=c(1,3), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
  par(new=T)
  plot(z, spec_var_time_all16[[u]], pch=0, xlim=c(1,3), ylim=c(0,1.2), axes=F, ann=F, col=colgrad2[u])
}
dev.off()
