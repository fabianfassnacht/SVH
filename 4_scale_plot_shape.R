require(rgdal)
setwd("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale")

grid_new <- readOGR(".", "grid")
setwd("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/1_Figure_1_scale/data")
load("full_data_sv_12_384_new.RData")
# grid_old <- readOGR(".", "grid")

# plot(grid_new, col=grid_new@data$id)
# text(grid_new@data$left, grid_new@data$top, grid_new@data$id)
# plot(grid_old, col=grid_new@data$id)
# text(grid_old@data$left, grid_old@data$top, grid_old@data$id)
# grid_new@data$id
# grid_old@data$id

cols <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666")
cols2 <- c("#f0027f", "#386cb0", "#666666", "#bf5b17","#beaed4","#7fc97f", "#ffff99","#fdc086")



png(width = 1000, height=800, filename = "spec_var_scale_effect_new.png", res=150)
par(mar=c(4,4,4,4))
plot(1:8, seq(1,80,10), ylim=c(0,100), xlim=c(1,400), col="white", axes=F, xlab="pixel size in cm", ylab="Spec. Var.", cex.lab=1.2)

for (i2 in 1:8){
  
  sub_data <- full_data[full_data$plot==i2,]
  par(new=T)
  plot(sub_data$spat_res, sub_data$MPD_PCboth, axes=F, ann=F, type="l", col=cols[i2], lwd=3, ylim=c(0,100), xlim=c(1,400))
  
}
box()
axis(1, cex=1.2)
axis(2, cex=1.2)
legend("bottomleft", col=cols2, legend=c("Pl. 1", "Pl. 2","Pl. 3",
                                        "Pl. 4","Pl. 5","Pl. 6",
                                        "Pl. 7","Pl. 8"), pch=15, ncol=2)

dev.off()