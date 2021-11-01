load("cl_4sp_100cl.RData")
cl_sp_100cl_agg2 <- cl_4sp_100cl

load("cl_sp_100cl_agg4.RData")
cl_sp_100cl_agg4 <- cl_sp_100cl

load("cl_sp_100cl_agg10.RData")
cl_sp_100cl_agg10 <- cl_sp_100cl


png(filename = "kmeans_simulations.png", height=1800, width=3400, res=300)

boxplot(rep(4,100), rep(8,100), rep(12,100),rep(16,100), rep(20,100),
        unlist(cl_sp_100cl_agg2)[1:100], 
        unlist(cl_sp_100cl_agg2)[101:200],
        unlist(cl_sp_100cl_agg2)[201:300],
        unlist(cl_sp_100cl_agg2)[301:400],
        unlist(cl_sp_100cl_agg2)[401:500], 
        unlist(cl_sp_100cl_agg4)[1:100], 
        unlist(cl_sp_100cl_agg4)[101:200],
        unlist(cl_sp_100cl_agg4)[201:300],
        unlist(cl_sp_100cl_agg4)[301:400],
        unlist(cl_sp_100cl_agg4)[401:500], 
        unlist(cl_sp_100cl_agg10)[1:100], 
        unlist(cl_sp_100cl_agg10)[101:200],
        unlist(cl_sp_100cl_agg10)[201:300],
        unlist(cl_sp_100cl_agg10)[301:400],
        unlist(cl_sp_100cl_agg10)[401:500],

        
        
        axes=F, ylab="number of kmeans clusters")
axis(1, at = 1:20, labels=c("4 species", "8 species", "12 species", "16 species", "20 species",
                            "4 species", "8 species", "12 species", "16 species", "20 species",
                            "4 species", "8 species", "12 species", "16 species", "20 species",
                            "4 species", "8 species", "12 species", "16 species", "20 species"), las=2)
axis(2)
box()
abline(v=5.5)
abline(v=10.5)
abline(v=15.5)
text(4, 1.1, paste0("cor=1"))
text(9, 1.1, paste0("cor=", round(cor(c(rep(4,100), rep(8,100), rep(12,100),rep(16,100), rep(20,100)), unlist(cl_sp_100cl_agg2), method="spearman"),2)))
text(14, 1.1, paste0("cor=", round(cor(c(rep(4,100), rep(8,100), rep(12,100),rep(16,100), rep(20,100)), unlist(cl_sp_100cl_agg4), method="spearman"),2)))
text(20, 1.1, paste0("cor=", round(cor(c(rep(4,100), rep(8,100), rep(12,100),rep(16,100), rep(20,100)), unlist(cl_sp_100cl_agg10), method="spearman"),2)))

text(0.1, 98, "A", cex=1.8)
text(6, 98, "B", cex=1.8)
text(11, 98, "C", cex=1.8)
text(16, 98, "D", cex=1.8)

dev.off()
