setwd("D:/1_SVH_paper_24_06_21/127_IAVS_special_issue_spectral_diversity/139_SVH_concept/Teja_SVH")

# load spectra and trait data
vars <- read.csv("canopy_kit_plant_functional_gradient.csv")

# check dates
vars$timestamp

# separate dataset into two annual datasets
vars_2016 <- vars[1:364,]
vars_2017 <- vars[365:659,]


colnames(vars)

## 1-17 traits
## 18-2117 reflectance


# 3 herbs (forb) / 3 grasses (graminoid) => 3 spectral signatures over for three time periods and 2 years

# grasses

# Calamagrostis epigejos
# Nardus stricta
# Agrostis capillaris 


# herbs

# Aegopodium podagraria
# Geum urbanum
# Festuca ovina

# both

# Campanula rotundifolia 
# Festuca ovina
# Lamium purpureum 
# Poa annua
# Trifolium pratense
# Stellaria media 




gr1_2016 <- vars_2016[vars_2016$species == "Calamagrostis epigejos",]
gr2_2016 <- vars_2016[vars_2016$species == "Nardus stricta",]
gr3_2016 <- vars_2016[vars_2016$species == "Agrostis capillaris",]
he1_2016 <- vars_2016[vars_2016$species == "Aegopodium podagraria",]
he2_2016 <- vars_2016[vars_2016$species == "Geum urbanum",]
he3_2016 <- vars_2016[vars_2016$species == "Festuca ovina",]


gr1_2016_refl <- gr1_2016[,18:2117]
gr2_2016_refl <- gr2_2016[-3,18:2117]
gr3_2016_refl <- gr3_2016[,18:2117]

he1_2016_refl <- he1_2016[,18:2117]
he2_2016_refl <- he2_2016[,18:2117]
he3_2016_refl <- he3_2016[,18:2117]

nm = 5
water_abs = c(1400:1500,1800:2000,2400:2500)
bands = seq(400,2500, by = nm)
bands = bands[-which(bands %in% water_abs==TRUE)]
dummy = bands-399


library(RColorBrewer)
colgrad = brewer.pal(9, "Purples")
colgrad2 = brewer.pal(6, "Set1")

###############################################
######### Figure intra-seasonal variability
###############################################

png(filename = "seasonal_var_species_new.png", width=2500, height=1500, res=300)

par(mfrow=c(2,3), oma=c(1,4,1,1), mar=c(4,4,2,2))

## grass 1

id = c(seq(2,8,1),8,8,9,9,8,8,rev(seq(2,8,1)))
plot(bands, gr1_2016_refl[1,dummy], main="Calamagrostis epigejos", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:17){
  par(new=T)
  plot(bands, gr1_2016_refl[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

## grass 2

id = c(seq(2,8,1),8,9,9,9,8,rev(seq(2,8,1)))
plot(bands, gr2_2016_refl[1,dummy], main="Nardus stricta", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:19){
  par(new=T)
  plot(bands, gr2_2016_refl[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

## grass 3

id = c(2,seq(2,8,1),9,9,9,rev(seq(2,8,1)),2)
plot(bands, gr3_2016_refl[1,dummy], main="Agrostis capillaris", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:17){
  par(new=T)
  plot(bands, gr3_2016_refl[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}




## herb 1

id = c(seq(2,8,1),8,8,9,9,8,8,rev(seq(2,8,1)))
plot(bands, he1_2016_refl[1,dummy], main="Aegopodium podagraria", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:20){
  par(new=T)
  plot(bands, he1_2016_refl[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

## herb 2

id = c(seq(2,8,1),8,8,9,9,8,8,rev(seq(2,8,1)))
plot(bands, he2_2016_refl[1,dummy], main="Geum urbanum", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:20){
  par(new=T)
  plot(bands, he2_2016_refl[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

## herb 3

id = c(seq(2,8,1),8,8,9,9,8,8,rev(seq(2,8,1)))
plot(bands, he3_2016_refl[1,dummy], main="Festuca ovina", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:20){
  par(new=T)
  plot(bands, he3_2016_refl[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

dev.off()



###############################################
######### Figure intra-seasonal variability
###############################################


gr1_2016$timestamp
gr2_2016$timestamp
gr3_2016$timestamp

he1_2016$timestamp
he2_2016$timestamp
he3_2016$timestamp


# dataset 31st of may

ds_31may <- rbind(gr1_2016_refl[1,],
                  gr2_2016_refl[1,],
                  gr3_2016_refl[1,],
                  he1_2016_refl[1,],
                  he2_2016_refl[1,],
                  he3_2016_refl[1,])

ds_20jun <- rbind(gr1_2016_refl[4,],
                 gr2_2016_refl[4,],
                 gr3_2016_refl[4,],
                 he1_2016_refl[4,],
                 he2_2016_refl[4,],
                 he3_2016_refl[4,])


ds_7jul <- rbind(gr1_2016_refl[6,],
                  gr2_2016_refl[6,],
                  gr3_2016_refl[6,],
                  he1_2016_refl[6,],
                  he2_2016_refl[6,],
                  he3_2016_refl[6,])

ds_11aug <- rbind(gr1_2016_refl[11,],
                 gr2_2016_refl[11,],
                 gr3_2016_refl[11,],
                 he1_2016_refl[11,],
                 he2_2016_refl[11,],
                 he3_2016_refl[11,])


ds_6sep <- rbind(gr1_2016_refl[15,],
                 gr2_2016_refl[15,],
                 gr3_2016_refl[15,],
                 he1_2016_refl[15,],
                 he2_2016_refl[15,],
                 he3_2016_refl[15,])


ds_20sep <- rbind(gr1_2016_refl[17,],
                 gr2_2016_refl[17,],
                 gr3_2016_refl[17,],
                 he1_2016_refl[17,],
                 he2_2016_refl[17,],
                 he3_2016_refl[17,])


png(filename = "seasonal_var_6species_date_new.png", width=2500, height=1500, res=300)

par(mfrow=c(2,3), oma=c(1,1,1,1), mar=c(4,4,2,2))

plot(bands, ds_31may[1,dummy], main="31st of May 2016", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:6){
  par(new=T)
  plot(bands, ds_31may[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad2[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

plot(bands, ds_20jun[1,dummy], main="20th of June 2016", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:6){
  par(new=T)
  plot(bands, ds_20jun[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad2[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

plot(bands, ds_7jul[1,dummy], main="7th of July 2016", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:6){
  par(new=T)
  plot(bands, ds_7jul[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad2[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}


plot(bands, ds_11aug[1,dummy], main="11th of August 2016", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:6){
  par(new=T)
  plot(bands, ds_11aug[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad2[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}


plot(bands, ds_6sep[1,dummy], main="6th of September 2016", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:6){
  par(new=T)
  plot(bands, ds_6sep[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad2[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

plot(bands, ds_20sep[1,dummy], main="20th of September 2016", ylim=c(0,0.9), xlim=c(300,2400), col="white", ylab="Reflectance [0-1 = 0-100%]", xlab="wavelength [nm]")
for (i in 1:6){
  par(new=T)
  plot(bands, ds_20sep[i,dummy], lwd="2", ylim=c(0,0.9), type="l", xlim=c(300,2400), col=colgrad2[id[i]], ylab="", xlab="", ann=F, axes=F)
  
}

dev.off()
