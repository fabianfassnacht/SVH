# function to aggregate synthetic rasters
aggloop <- function (inlist, fac) {

  out1 <- list()
  for (i in 1:100){
    
    r <- inlist[[i]]
    r1 <- aggregate(r, fac)
    out1[[i]] <- r1
    
  }
  out1
}

# function to assign extents for mosaicing
assign_ex <- function(inlist, exes){
  
  out2 <- list()
  for (i2 in 1:500){
    
    r <- inlist[[i2]]
    extent(r) <- extent(as.numeric(exes[i2,]))
    out2[[i2]] <- r 
    
  }
  out2
}

# function to get number of kmeans classes


getsvpca <- function(kmap, exes){
  
  out3 <-list()  
  
  for (i3 in 1:500){
    
    m1 <- crop(kmap, extent(as.numeric(exes[i3,])))
    
    image_1 <- m1$PC1
    results <-mean(dist(getValues(image_1)[!is.na(getValues(image_1))]), na.rm=T)
    image_2<-m1$PC2
    results2<-mean(dist(getValues(image_2)[!is.na(getValues(image_2))]), na.rm=T)
    val_1<- getValues(image_1)[!is.na(getValues(image_1))]
    val_2<- getValues(image_2)[!is.na(getValues(image_2))]
    dist_both<-sqrt((dist(val_1)^2)+(dist(val_2)^2))
    results3 <-mean(dist_both, na.rm=T)
    
    # return results
    print(results3)
    out3[[i3]] <- results3
    
  }
  out3
}

require(raster)
require(RStoolbox)

# load synthetic raster files
setwd("E:/Fabian_specvar/SVH")
load("spec4_100ras.RData")
load("spec8_100ras.RData")
load("spec12_100ras.RData")
load("spec16_100ras.RData")
load("spec20_100ras.RData")


# prepare extents to mosaic all raster next to each other
exes <- data.frame(seq(0,49900,100), seq(100,50000,100),rep(0,500), rep(100,500))
colnames(exes) <- c("xmin", "xmax", "ymin", "ymax")

# run a loop to assign the extents to each raster and store
# raster objects with extents to new list

agg2_all_sp <- c(spec4_res, spec8_res, spec12_res, spec16_res, spec20_res)
agg2_all_sp_ex <- assign_ex(agg2_all_sp, exes)



# mosaic all 100 rasters to one file
agg2_all_sp_ex$fun <- "mean"
r_full_all_sp_ex <- do.call(mosaic, agg2_all_sp_ex)

setwd("E:/KIT_Forschung/1_SVH_paper_03_10_21/2_manuscript/Revisions_AVS/1_code")
writeRaster(r_full_all_sp_ex, filename="mosaiced_all_species_agg2_for_kmeans.tif", format="GTiff")

summary(r_full_all_sp_ex)

table(is.na(values(r_full_all_sp_ex)))

# calculate PCA of mosaic
pca_raster<-rasterPCA(r_full_all_sp_ex, nSamples=100)

# take principal component raster
kmap <- pca_raster$map

# check spectral variation for first two components of each individual image
sv_pca <- getsvpca(kmap, exes)

save(sv_pca, file="pca_sv.RData" )

setwd("E:/KIT_Forschung/1_SVH_paper_03_10_21/2_manuscript/Revisions_AVS/1_code")
png(filename = "kmeans_boxplots_agg10_all_species.png", height=1000, width=750, res=100)
par(mfrow=c(2,1))

boxplot(unlist(sv_pca)[1:100], 
        unlist(sv_pca)[101:200],
        unlist(sv_pca)[201:300],
        unlist(sv_pca)[301:400],
        unlist(sv_pca)[401:500], axes=F, ylab="number of kmeans clusters", ylim=c(0,1.2))
axis(1, at = 1:5, labels=c("4 species", "8 species", "12 species", "16 species", "20 species"))
axis(2)
box()

dev.off()
