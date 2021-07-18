require(beanplot)

setwd("D:/1_SVH_paper_29_06_21/5_simulation_graph_new_21_07_13/sim_ras/sim_rasters")
fi <- list.files()
fi

## load spectral variation

## prepare raw data

load("sv_4spec_100iter.RData")
load("sv_8spec_100iter.RData")
load("sv_12spec_100iter.RData")
load("sv_16spec_100iter.RData")
load("sv_20spec_100iter.RData")
load("sv_21spec_100iter.RData")

raw <- data.frame(unlist(sv_4spec),
                  unlist(sv_8spec),
                  unlist(sv_12spec),
                  unlist(sv_16spec),
                  unlist(sv_20spec),
                  unlist(sv_21spec))

rm(sv_4spec, sv_8spec, sv_12spec, sv_16spec, sv_20spec, sv_21spec)

colnames(raw) <- c("4", "8","12", "16","20", "21")


## prepare aggregated data (4 pixels aggregated)

load("sv_4spec_100iter_agg2_v2.RData")
load("sv_8spec_100iter_agg2_v.RData")
load("sv_12spec_100iter_agg2_v.RData")
load("sv_16spec_100iter_agg2_v.RData")
load("sv_20spec_100iter_agg2_v.RData")
load("sv_21spec_100iter_agg2_v.RData")


agg2 <- data.frame(unlist(sv_4spec),
                   unlist(sv_8spec),
                   unlist(sv_12spec),
                   unlist(sv_16spec),
                   unlist(sv_20spec),
                   unlist(sv_21spec))

rm(sv_4spec, sv_8spec, sv_12spec, sv_16spec, sv_20spec, sv_21spec)

colnames(agg2) <- c("4", "8","12", "16","20", "21")


## prepare aggregated data (4 pixels aggregated)

load("sv_4spec_agg4_100iter.RData")
load("sv_8spec_agg4_100iter.RData")
load("sv_12spec_agg4_100iter.RData")
load("sv_16spec_agg4_100iter.RData")
load("sv_20spec_agg4_100iter.RData")
load("sv_21spec_agg4_100iter.RData")


agg4 <- data.frame(unlist(sv_4spec),
                  unlist(sv_8spec),
                  unlist(sv_12spec),
                  unlist(sv_16spec),
                  unlist(sv_20spec),
                  unlist(sv_21spec))

rm(sv_4spec, sv_8spec, sv_12spec, sv_16spec, sv_20spec, sv_21spec)

colnames(agg4) <- c("4", "8","12", "16","20", "21")



## prepare aggregated data (10 pixels aggregated)

load("sv_4spec_agg10_100iter.RData")
load("sv_8spec_agg10_100iter.RData")
load("sv_12spec_agg10_100iter.RData")
load("sv_16spec_agg10_100iter.RData")
load("sv_20spec_agg10_100iter.RData")
load("sv_21spec_agg10_100iter.RData")


agg10 <- data.frame(unlist(sv_4spec),
                   unlist(sv_8spec),
                   unlist(sv_12spec),
                   unlist(sv_16spec),
                   unlist(sv_20spec),
                   unlist(sv_21spec))

rm(sv_4spec, sv_8spec, sv_12spec, sv_16spec, sv_20spec, sv_21spec)

colnames(agg10) <- c("4", "8","12", "16","20", "21")


###############################
## now plot everything
###############################


all <- c(raw, agg2, agg4, agg10)
cols <- c(rep("grey80", 6),rep("grey60", 6),rep("grey40", 6),rep("grey20", 6)) 


### final plot



specnr <- c(rep(4,100),rep(8,100),rep(12,100),rep(16,100),rep(20,100),rep(21,100))
specvar <- c(raw[,1],raw[,2],raw[,3],raw[,4],raw[,5],raw[,6])
plot(specnr, specvar)
cor(specnr, specvar)


specvar2 <- c(agg2[,1],agg2[,2],agg2[,3],agg2[,4],agg2[,5],agg2[,6])
plot(specnr, specvar2)
cor(specnr, specvar2)

specvar3 <- c(agg4[,1],agg4[,2],agg4[,3],agg4[,4],agg4[,5],agg4[,6])
plot(specnr, specvar3)
cor(specnr, specvar3)

specvar4 <- c(agg10[,1],agg10[,2],agg10[,3],agg10[,4],agg10[,5],agg10[,6])
plot(specnr, specvar4)
cor(specnr, specvar4)

plot()

dev.off()


png(filename = "simulations_21_07_13_v1.png", height=1800, width=3400, res=300)
boxplot(all, ylim=c(0, 1.2), col=cols, xlab="Number of Species", ylab="Spec. Var.", cex.axis=1.2, cex.lab=1.4, notch=T)
abline(v=6.5)
abline(v=12.5)
abline(v=18.5)
text(5, 1.1, paste0("cor=", round(cor(specnr, specvar, method="spearman"),2)))
text(11, 1.1, paste0("cor=", round(cor(specnr, specvar2, method="spearman"),2)))
text(17, 1.1, paste0("cor=", round(cor(specnr, specvar3, method="spearman"),2)))
text(23, 1.1, paste0("cor=", round(cor(specnr, specvar4, method="spearman"),2)))

text(0.2, 1.18, "A", cex=2)
text(7.2, 1.18, "B", cex=2)
text(13.2, 1.18, "C", cex=2)
text(19.2, 1.18, "D", cex=2)

dev.off()




#png(filename = "scatter_simulations_21_07_13.png", height=2000, width=2000, res=300)

#par(mfrow=c(2,2))
#par(mar=c(4.5,4.5,2,2))
#set.seed(23)
#plot(specnr+sample(seq(0,1,0.001),100), specvar, main="100 x 100", xlim=c(4,22), ylim=c(0, 1.3), ylab="Spec. Variability", xlab = "Species Count", cex.lab=1.2, cex.axis=1.2)
#text(18, 0.15, paste0("cor=", round(cor(specnr, specvar, method="spearman"),2)))
#plot(specnr+sample(seq(0,1,0.001),100), specvar2, main="50 x 50", xlim=c(4,22), ylim=c(0, 1), ylab="Spec. Variability", xlab = "Species Count", cex.lab=1.2, cex.axis=1.2)
#text(18, 0.85, paste0("cor=", round(cor(specnr, specvar2, method="spearman"),2)))
#plot(specnr+sample(seq(0,1,0.001),100), specvar3, main="25 x 25", xlim=c(4,22), ylim=c(0, 1), ylab="Spec. Variability", xlab = "Species Count", cex.lab=1.2, cex.axis=1.2)
#text(18, 0.85, paste0("cor=", round(cor(specnr, specvar3, method="spearman"),2)))
#plot(specnr+sample(seq(0,1,0.001),100), specvar4, main="10 x 10", xlim=c(4,22), ylim=c(0, 1), ylab="Spec. Variability", xlab = "Species Count", cex.lab=1.2, cex.axis=1.2)
#text(18, 0.85, paste0("cor=", round(cor(specnr, specvar4, method="spearman"),2)))

#dev.off()

