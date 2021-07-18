#samp:  community matrix; sites in lines, species in columns
#dis:   dissimilarity matrix
#type:  "both" for results with abundances weighted and non-weighted
#       "abundance" for results with abundances weighted
#       "presence" for results with abundances non-weighted
  
melodic <- function(samp,dis,type="both"){
	if(class(samp)!="matrix"){samp <- as.matrix(samp)}
  	if(class(dis)!="matrix"){dis <- as.matrix(dis)}
  	if(is.null(colnames(samp)) | is.null(colnames(dis)) ){
  		stop("Both samp and dis must have colnames.\n")
  	}
  	N<-dim(samp)[1]
	melodic<-list()
	if (type=="both"){
	    melodic$abundance<-list()
	    melodic$abundance$mpd<-melodic$abundance$rao<-melodic$abundance$simpson<-numeric(N)
	    melodic$presence<-list()
	    melodic$presence$mpd<-melodic$presence$rao<-melodic$presence$simpson<-numeric(N)
  	}
	if (type=="abundance"){ 
    	melodic$abundance<-list()
	    melodic$abundance$mpd<-melodic$abundance$rao<-melodic$abundance$simpson<-numeric(N)
  	}
	if (type=="presence"){ 
    	melodic$presence<-list()
	    melodic$presence$mpd<-melodic$presence$rao<-melodic$presence$simpson<-numeric(N)
  	}
	for (i in 1:N){
    	sppInSample<-names(samp[i,samp[i, ]>0])
	  	melodic$richness[i]<-rowSums(samp>0)[i]
	  	if (length(sppInSample)>1){
      		sample.dis<-dis[sppInSample,sppInSample]
			abund.w<-numeric(length(sppInSample))
			if (type=="both" | type=="abundance"){
  				abund.w <- samp[i , sppInSample] / sum(samp[i , sppInSample])
		        sample.weights <- outer(abund.w , abund.w)
		        melodic$abundance$mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
			    melodic$abundance$rao[i] <- sum(sample.weights * sample.dis)
		    	melodic$abundance$simpson[i] <- sum(2*sample.weights[lower.tri(sample.weights)])
			} 	
    		if (type=="both" | type=="presence"){
		    	abund.nw <- rep(1 , length(sppInSample)) / length(sppInSample)
		        sample.weights <- outer(abund.nw , abund.nw)
		        #melodic$presence$mntd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
		        melodic$presence$mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
		        melodic$presence$rao[i] <- sum(sample.weights * sample.dis)
		        melodic$presence$simpson[i] <- sum(2*sample.weights[lower.tri(sample.weights)])
			}	
		}	else {
			if (type=="both" | type=="abundance"){
		        melodic$abundance$mpd[i] <- NA
		        melodic$abundance$rao[i] <- melodic$abundance$simpson[i] <-0
	    	}
      		if (type=="both" | type=="presence"){
        		melodic$presence$mpd[i] <- NA
		        melodic$presence$rao[i] <- melodic$presence$simpson[i] <-0
			}
    	}
    }  	
	out<-melodic
	return(out)	
}

############
############ EXAMPLES
library(picante)
#data(phylocom)
#distances<-cophenetic(phylocom$phylo)/max(cophenetic(phylocom$phylo))
#test<-mpd(phylocom$sample,distances, abundance.weighted =T) ### Rao values
#test.ab<-melodic(phylocom$sample,distances, type="abundance") ### mpd, rao and simpson values, abundance weighted
#test.pr<-melodic(phylocom$sample,distances, type="presence") ### mpd, rao and simpson values, not abundance weighted
#test.both<-melodic(phylocom$sample,distances, type="both") ### mpd, rao and simpson values, both  weighted and non-weighted



setwd("D:/1_SVH_paper_29_06_21/127_IAVS_special_issue_spectral_diversity/Teja_SVH")
data <- read.csv("canopy_kit_plant_functional_gradient.csv")
names <- data$species[2:21]
traits <- data[2:21,4:15]

traits[,1] <- as.numeric(traits[,1])
traits[,2] <- as.numeric(traits[,2])
traits[,3] <- as.numeric(traits[,3])
traits[,4] <- as.numeric(traits[,4])
traits[,5] <- as.numeric(traits[,5])
traits[,6] <- as.numeric(traits[,6])
traits[,7] <- as.numeric(traits[,7])
traits[,8] <- as.numeric(traits[,8])
traits[,9] <- as.numeric(traits[,9])
traits[,10] <- as.numeric(traits[,10])
traits[,11] <- as.numeric(traits[,11])
traits[,12] <- as.numeric(traits[,12])



setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_17_fd/sim_ras_07_15_fd")
list.files()


### prepare data to calculate

# get FD for 4 species

load("spec4_100cov.RData")
load("spec4_100spec.RData")

data_t4 <- matrix(0, nrow=100, ncol=20)

colnames(data_t4) <- names
rownames(data_t4) <- seq(1,100,1)

for (i4 in 1:100){
  
  cd <- spec4_cov[[i4]]
  data_t4[i4,colnames(data_t4)%in%colnames(cd)] <- as.numeric(cd)
  
}

# get FD for 8 species

load("spec8_100cov.RData")
load("spec8_100spec.RData")

data_t8 <- matrix(0, nrow=100, ncol=20)

colnames(data_t8) <- names
rownames(data_t8) <- seq(1,100,1)

for (i8 in 1:100){
  
  cd <- spec8_cov[[i8]]
  data_t8[i8,colnames(data_t8)%in%colnames(cd)] <- as.numeric(cd)
  
}


# get FD for 12 species

load("spec12_100cov.RData")
load("spec12_100spec.RData")

data_t12 <- matrix(0, nrow=100, ncol=20)

colnames(data_t12) <- names
rownames(data_t12) <- seq(1,100,1)

for (i12 in 1:100){
  
  cd <- spec12_cov[[i12]]
  data_t12[i12,colnames(data_t12)%in%colnames(cd)] <- as.numeric(cd)
  
}

# get FD for 16 species

load("spec16_100cov.RData")
load("spec16_100spec.RData")

data_t16 <- matrix(0, nrow=100, ncol=20)

colnames(data_t16) <- names
rownames(data_t16) <- seq(1,100,1)

for (i16 in 1:100){
  
  cd <- spec16_cov[[i16]]
  data_t16[i16,colnames(data_t16)%in%colnames(cd)] <- as.numeric(cd)
  
}


data_tall <- rbind(data_t4,data_t8,data_t12,data_t16)

## now get functional diversity values for all synhetic images
rownames(traits) <- names
tr_dis <- as.matrix(dist(traits))

#fd_all <-mpd(data_tall, tr_dis, abundance.weighted =T)
fd_all <-melodic(data_tall, tr_dis) 
?mpd


## now load spectal variation data for all synthetic images

######
###### 100 x 100 pixel images
######


load("sv_4spec_100iter.RData")
load("sv_8spec_100iter.RData")
load("sv_12spec_100iter.RData")
load("sv_16spec_100iter.RData")
load("sv_20spec_100iter.RData")

sv_4spec_1 <- unlist(sv_4spec)
sv_8spec_1 <- unlist(sv_8spec)
sv_12spec_1 <- unlist(sv_12spec)
sv_16spec_1 <- unlist(sv_16spec)

sv_all_orig <- c(sv_4spec_1, sv_8spec_1, sv_12spec_1, sv_16spec_1)

load("sv_4spec_100iter_agg2_v2.RData")
load("sv_8spec_100iter_agg2_v.RData")
load("sv_12spec_100iter_agg2_v.RData")
load("sv_16spec_100iter_agg2_v.RData")
load("sv_20spec_100iter_agg2_v.RData")

sv_4spec_1 <- unlist(sv_4spec)
sv_8spec_1 <- unlist(sv_8spec)
sv_12spec_1 <- unlist(sv_12spec)
sv_16spec_1 <- unlist(sv_16spec)

sv_all_agg2 <- c(sv_4spec_1, sv_8spec_1, sv_12spec_1, sv_16spec_1)


load("sv_4spec_agg4_100iter.RData")
load("sv_8spec_agg4_100iter.RData")
load("sv_12spec_agg4_100iter.RData")
load("sv_16spec_agg4_100iter.RData")

sv_4spec_1 <- unlist(sv_4spec)
sv_8spec_1 <- unlist(sv_8spec)
sv_12spec_1 <- unlist(sv_12spec)
sv_16spec_1 <- unlist(sv_16spec)

sv_all_agg4 <- c(sv_4spec_1, sv_8spec_1, sv_12spec_1, sv_16spec_1)


load("sv_4spec_agg10_100iter.RData")
load("sv_8spec_agg10_100iter.RData")
load("sv_12spec_agg10_100iter.RData")
load("sv_16spec_agg10_100iter.RData")
load("sv_20spec_agg10_100iter.RData")

sv_4spec_1 <- unlist(sv_4spec)
sv_8spec_1 <- unlist(sv_8spec)
sv_12spec_1 <- unlist(sv_12spec)
sv_16spec_1 <- unlist(sv_16spec)
sv_20spec_1 <- unlist(sv_20spec)

sv_all_agg10 <- c(sv_4spec_1, sv_8spec_1, sv_12spec_1, sv_16spec_1)


require(RColorBrewer)
colgrad2 = brewer.pal(4, "Set1")
cols <- c(rep(colgrad2[1], 100), rep(colgrad2[2], 100), rep(colgrad2[3], 100), rep(colgrad2[4], 100))

png(filename = "functional_diversity_21_07_17_v1.png", width=2800, height=2400, res=300)

par(mfrow=c(2,2))

plot(sv_all_orig, fd_all$abundance$mpd, xlim=c(0,1.3), col=cols, main="100 × 100 pixels", ylab = "Functional Diversity", xlab = "Spectral Variation", cex.lab=1.3, cex.axis=1.3)
plot(sv_all_agg2, fd_all$abundance$mpd, xlim=c(0,1.3), col=cols, main="50 × 50 pixels", ylab = "Functional Diversity", xlab = "Spectral Variation", cex.lab=1.3, cex.axis=1.3)
plot(sv_all_agg4, fd_all$abundance$mpd, xlim=c(0,1.3), col=cols, main="25 × 25 pixels", ylab = "Functional Diversity", xlab = "Spectral Variation", cex.lab=1.3, cex.axis=1.3)
plot(sv_all_agg10, fd_all$abundance$mpd, xlim=c(0,1.3), col=cols, main="10 × 10 pixels", ylab = "Functional Diversity", xlab = "Spectral Variation", cex.lab=1.3, cex.axis=1.3)

legend("topright", legend=c("4 species", "8 species", "12 species", "16 species"), col = colgrad2[1:4], pch=22, cex=1.3)

dev.off()  
  
  