require(raster)
require(lattice)
require(gridExtra)
options(scipen=999)

setwd("E:/KIT_Forschung/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale")

img <- stack("UAV_image_1.tif")
# img2 <- aggregate(img, fact=2, fun=mean)
# img3 <- aggregate(img, fact=4, fun=mean)
# img4 <- aggregate(img, fact=8, fun=mean)
# img5 <- aggregate(img, fact=16, fun=mean)
# img6 <- aggregate(img, fact=32, fun=mean)
# img7 <- aggregate(img, fact=64, fun=mean)
# img8 <- aggregate(img, fact=128, fun=mean)
# 
# writeRaster(img2, filename = "uav_6.tif", format="GTiff")
# writeRaster(img3, filename = "uav_12.tif", format="GTiff")
# writeRaster(img4, filename = "uav_24.tif", format="GTiff")
# writeRaster(img5, filename = "uav_48.tif", format="GTiff")
# writeRaster(img6, filename = "uav_92.tif", format="GTiff")
# writeRaster(img7, filename = "uav_192.tif", format="GTiff")
# writeRaster(img8, filename = "uav_384.tif", format="GTiff")

img2 <- stack("uav_6.tif")
img3 <- stack("uav_12.tif")
img4 <- stack("uav_24.tif")
img5 <- stack("uav_48.tif")
img6 <- stack("uav_96.tif")
img7 <- stack("uav_192.tif")
img8 <- stack("uav_384.tif")

img_v <- values(img[[1]])
img2_v <- values(img2[[1]])
img3_v <- values(img3[[1]])
img4_v <- values(img4[[1]])
img5_v <- values(img5[[1]])
img6_v <- values(img6[[1]])
img7_v <- values(img7[[1]])
img8_v <- values(img8[[1]])


p1 <- histogram(img_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="% of all values", xlab="DN blue band - 3 cm pixel", breaks=22)
p2 <- histogram(img2_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="", xlab="6 cm pixel", breaks=22)
p3 <- histogram(img3_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="", xlab="12 cm pixel", breaks=22)
p4 <- histogram(img4_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="", xlab="24 cm pixel", breaks=22)
p5 <- histogram(img5_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="% of all value", xlab="48 cm pixel", breaks=22)
p6 <- histogram(img6_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="", xlab="96 cm pixel", breaks=22)
p7 <- histogram(img7_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="", xlab="192 cm pixel", breaks=22)
p8 <- histogram(img8_v, col="grey", xlim=c(0,255), ylim=c(0,18), ylab="", xlab="384 cm pixel", breaks=22)


# Plot prints
dev.off()
png(width=1400, height=1000, res=150, filename="Fig_1_hist_scale.png")
print(p1, split = c(1, 1, 4, 2), more = TRUE)
print(p2, split = c(2, 1, 4, 2), more = TRUE)
print(p3, split = c(3, 1, 4, 2), more = TRUE)
print(p4, split = c(4, 1, 4, 2), more = TRUE)
print(p5, split = c(1, 2, 4, 2), more = TRUE)
print(p6, split = c(2, 2, 4, 2), more = TRUE)
print(p7, split = c(3, 2, 4, 2), more = TRUE)
print(p8, split = c(4, 2, 4, 2), more = FALSE)
dev.off()



png(width=1400, height=1000, res=150, filename="Fig_1_hist_scale_rgb.png")

par(mfrow=c(1,3), oma=c(0,0,0,0), mar=c(1,1,1,1))

plotRGB(img, r=3, g=2, b=1, stretch="hist", main="pixel size 3 cm", margins=T)
plotRGB(img5, r=3, g=2, b=1, stretch="hist", main="pixel size 48 cm", margins=T)
plotRGB(img8, r=3, g=2, b=1, stretch="hist", main="pixel size 384 cm", margins=T)

dev.off()

