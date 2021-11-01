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

getkmcl <- function(kmap, exes){

  out3 <-list()  
  
  for (i3 in 1:500){
    
    m1 <- crop(kmap, extent(as.numeric(exes[i3,])))
    l <- length(table(values(m1)))
    out3[[i3]] <- l
    
  }
  out3
}

require(raster)
require(RStoolbox)

# load synthetic raster files
setwd("E:/KIT_Forschung/1_SVH_paper_18_07_21/5_simulation_graph_new_21_07_17_fd/sim_ras_07_15_fd")
load("spec4_100ras.RData")
load("spec8_100ras.RData")
load("spec12_100ras.RData")
load("spec16_100ras.RData")
load("spec20_100ras.RData")

# aggregate to 50 x 50 pixels
agg2_4sp <- aggloop(spec4_res, 2) 
agg2_8sp <- aggloop(spec8_res, 2) 
agg2_12sp <- aggloop(spec12_res, 2) 
agg2_16sp <- aggloop(spec16_res, 2) 
agg2_20sp <- aggloop(spec20_res, 2) 

# prepare extents to mosaic all raster next to each other
exes <- data.frame(seq(0,24950,50), seq(50,25000,50),rep(0,500), rep(50,500))
colnames(exes) <- c("xmin", "xmax", "ymin", "ymax")

# run a loop to assign the extents to each raster and store
# raster objects with extents to new list

agg2_all_sp <- c(agg2_4sp, agg2_8sp, agg2_12sp, agg2_16sp, agg2_20sp)

agg2_all_sp_ex <- assign_ex(agg2_all_sp, exes)



# mosaic all 100 rasters to one file
agg2_all_sp_ex$fun <- "mean"

r_full_all_sp_ex <- do.call(mosaic, agg2_all_sp_ex)
setwd("E:/KIT_Forschung/1_SVH_paper_03_10_21/2_manuscript/Revisions_AVS/1_code")
save(r_full_all_sp_ex, file="mosaiced_all_species_agg2_for_kmeans.Rdata")

# now run kemans algorithm on the new raster
set.seed(12)
#unC_all_sp_20cl <- unsuperClass(r_full_all_sp_ex, nSamples = 10000, nClasses = 20, nStarts = 5)
#km_m_all_sp_20cl <- unC_all_sp_20cl$map
unC_all_sp_100cl <- unsuperClass(r_full_all_sp_ex, nSamples = 10000, nClasses = 100, nStarts = 5)
km_m_all_sp_100cl <- unC_all_sp_100cl$map

# check number of classes for each kmeans raster file

save(km_m_all_sp_100cl, file="km_m_all_sp_100cl.RData")

#cl_4sp_20cl <- getkmcl(km_m_all_sp_20cl, exes)
cl_4sp_100cl <- getkmcl(km_m_all_sp_100cl, exes)

save(cl_4sp_100cl, file="cl_4sp_100cl.RData")

setwd("E:/KIT_Forschung/1_SVH_paper_03_10_21/2_manuscript/Revisions_AVS/1_code")
png(filename = "kmeans_boxplots_agg2_all_species.png", height=1000, width=750, res=100)
par(mfrow=c(2,1))

boxplot(unlist(cl_4sp_100cl)[1:100], 
        unlist(cl_4sp_100cl)[101:200],
        unlist(cl_4sp_100cl)[201:300],
        unlist(cl_4sp_100cl)[301:400],
        unlist(cl_4sp_100cl)[401:500], axes=F, ylab="number of kmeans clusters")
axis(1, at = 1:5, labels=c("4 species", "8 species", "12 species", "16 species", "20 species"))
axis(2)
box()

dev.off()
