require(rgdal)
require(RStoolbox)
require(raster)

# function to calculate spectral variation from first
# two PCA axes

get_specVarPCA <- function(img, grid, cm){
  
  #calculate PCA------
  
  pca_raster<-rasterPCA(img[[1:3]])
  class(pca_raster)
  #plot(pca_raster$map$PC1)
  #summary(pca_raster$model)
  
  values_PCA<-as.data.frame(raster::extract(pca_raster$map, c(1:ncell(pca_raster$map))))
  #hist(na.omit(values_PCA$PC1))
  
  
  #for each subplot
  results<-data.frame(plot=c(1:8), MPD_PC1=NA, MPD_PC2=NA, MPD_PCboth=NA, spat_res=rep(cm))
  
  for(i in 1:8){
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
setwd("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/")
net <- readOGR(".", "grid")

img_384 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_384.tif") #RasterStack #PC
img_192 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_192.tif") #RasterStack #PC
img_96 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_96.tif") #RasterStack #PC
img_48 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_48.tif") #RasterStack #PC
img_24 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_24.tif") #RasterStack #PC
img_12 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_12.tif") #RasterStack #PC
img_6 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/uav_6.tif") #RasterStack #PC
img_3 <-stack("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/UAV_image_1.tif") #RasterStack #PC

sv3 <- get_specVarPCA(img_3, net, 3)
sv6 <- get_specVarPCA(img_6, net, 6)
sv12 <- get_specVarPCA(img_12, net, 12)
sv24 <- get_specVarPCA(img_24, net, 24)
sv48 <- get_specVarPCA(img_48, net, 48)
sv96 <- get_specVarPCA(img_96, net, 96)
sv192 <- get_specVarPCA(img_192, net, 192)
sv384 <- get_specVarPCA(img_384, net, 384)


full_data <- rbind(sv3, sv6, sv12, sv24, sv48, sv96, sv192, sv384)
full_data <- rbind(sv12, sv24, sv48, sv96, sv192, sv384)
cols <- c("#b35806", "#e08214", "#fdb863", "#fee0b6", "#d8daeb", "#b2abd2", "#8073ac", "#542788")

save(full_data, file = "full_data_sv_12_384.RData")


png(width = 1000, height=800, filename = "spec_var_scale_effect.png", res=150)
par(mar=c(4,4,4,4))
plot(1:8, seq(1,80,10), ylim=c(0,100), xlim=c(1,400), col="white", axes=F, xlab="pixel size in cm", ylab="spec. var. (mean dist. pixels / PCA axes 1+2)")

for (i2 in 1:8){
  
  sub_data <- full_data[full_data$plot==i2,]
  par(new=T)
  plot(sub_data$spat_res, sub_data$MPD_PCboth, axes=F, ann=F, type="l", col=cols[i2], lwd=3, ylim=c(0,100), xlim=c(1,400))
  
}
box()
axis(1)
axis(2)
dev.off()
