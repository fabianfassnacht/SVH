require(rgdal)
require(RStoolbox)
require(raster)

# load raster data

setwd("E:/KIT_Forschung/1_SVH_paper/5_sims_new/sim_ras_07_15_fd")

load("spec4_100ras.RData")
load("spec8_100ras.RData")
load("spec12_100ras.RData")
load("spec16_100ras.RData")
load("spec20_100ras.RData")
load("spec21_100ras.RData")



# function to calculate spectral variation from first
# two PCA axes

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

# now calculate spectral variation from all 10 rasters of each number of species

### 4 species

niter = 100

sv_4spec <- list()

for (i in 1:niter){
  
  img <- spec4_res[[i]]
  img <- aggregate(img, 2)
  specvar <- get_specVarPCA(img)
  sv_4spec[[i]] <- specvar
  print(i)
  
}

save(sv_4spec, file="sv_4spec_100iter_agg2_v2.RData")

### 8 species

sv_8spec <- list()

for (i in 1:niter){
  
  img <- spec8_res[[i]]
  img <- aggregate(img, 2)
  specvar <- get_specVarPCA(img)
  sv_8spec[[i]] <- specvar
  print(i)
  
}

save(sv_8spec, file="sv_8spec_100iter_agg2_v.RData")
### 12 species

sv_12spec <- list()

for (i in 1:niter){
  
  img <- spec12_res[[i]]
  img <- aggregate(img, 2)
  specvar <- get_specVarPCA(img)
  sv_12spec[[i]] <- specvar
  print(i)
  
}

save(sv_12spec, file="sv_12spec_100iter_agg2_v.RData")

### 16 species

sv_16spec <- list()

for (i in 1:niter){
  
  img <- spec16_res[[i]]
  img <- aggregate(img, 2)
  specvar <- get_specVarPCA(img)
  sv_16spec[[i]] <- specvar
  print(i)
  
}

save(sv_16spec, file="sv_16spec_100iter_agg2_v.RData")

### 20 species

sv_20spec <- list()

for (i in 1:niter){
  
  img <- spec20_res[[i]]
  img <- aggregate(img, 2)
  specvar <- get_specVarPCA(img)
  sv_20spec[[i]] <- specvar
  print(i)
  
}

save(sv_20spec, file="sv_20spec_100iter_agg2_v.RData")

### 22 species

sv_21spec <- list()

for (i in 1:niter){
  
  img <- spec21_res[[i]]
  img <- aggregate(img, 2)
  specvar <- get_specVarPCA(img)
  sv_21spec[[i]] <- specvar
  print(i)
  
}

save(sv_21spec, file="sv_21spec_100iter_agg2_v.RData")



sv_4 <- unlist(sv_4spec)
sv_8 <- unlist(sv_8spec)
sv_12 <- unlist(sv_12spec)
sv_16 <- unlist(sv_16spec)
sv_20 <- unlist(sv_20spec)
sv_22 <- unlist(sv_22spec)


res_1 <- data.frame(sv_4, sv_8, sv_12, sv_16, sv_20, sv_22)

boxplot(res_1)



