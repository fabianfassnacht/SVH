require(rgdal)
require(RStoolbox)
require(raster)

# function to calculate spectral variation from first
# two PCA axes

get_specVarPCA <- function(img, grid, cm){
  
  #calculate PCA------
  
  pca_raster<-rasterPCA(img[[-c(9:10)]])
  class(pca_raster)
  #plot(pca_raster$map$PC1)
  #summary(pca_raster$model)
  
  values_PCA<-as.data.frame(raster::extract(pca_raster$map, c(1:ncell(pca_raster$map))))
  #hist(na.omit(values_PCA$PC1))
  
  
  #for each subplot
  results<-data.frame(plot=c(1:7), MPD_PC1=NA, MPD_PC2=NA, MPD_PCboth=NA, spat_res=rep(cm))
  
  for(i in 1:7){
    subplot<-grid[i,]
    image_1<-mask(crop(pca_raster$map$PC1, subplot), subplot)
    results[results$plot==i,"MPD_PC1"]<-mean(dist(getValues(image_1)[!is.na(getValues(image_1))]), na.rm=T)
    # results[results$plot==i,"sd_PC1"]<-sd(getValues(image_1)[!is.na(getValues(image_1))], na.rm=T)
    image_2<-mask(crop(pca_raster$map$PC2, subplot), subplot)
    results[results$plot==i,"MPD_PC2"]<-mean(dist(getValues(image_2)[!is.na(getValues(image_2))]), na.rm=T)
    # results[results$plot==i,"sd_PC2"]<-sd(getValues(image_2)[!is.na(getValues(image_2))], na.rm=T)
    val_1<- getValues(image_1)[!is.na(getValues(image_1))]
    val_2<- getValues(image_2)[!is.na(getValues(image_2))]
    dist_both<-sqrt((dist(val_1)^2)+(dist(val_2)^2))
    results[results$plot==i,"MPD_PCboth"]<-mean(dist_both, na.rm=T)
    print(i)
    
  }
  results
}


#Load data---------------
setwd("D:/1_SVH_paper_24_06_21/Garchinger_Heide")
net <- readOGR(".", "heide_32632_v1")

img_3 <-stack("D:/1_SVH_paper_24_06_21/Garchinger_Heide/subset_0_of_S2B_MSIL1C_20190830T102029_N0208_R065_T32UPU_20190830T130621_resampled.tif") #RasterStack #PC

sv3 <- get_specVarPCA(img_3, net, 1000)

sv3
