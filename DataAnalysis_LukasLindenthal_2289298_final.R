# # Script final assigment "Exercise Bruchsal"
# Author: Lukas Lindenthal 
# Matr. nr,: 2289298
# date: 27.03.2022

#set working directory (please adapt)
setwd("C:\\")

# install packages (if needed)
# install.packages("dplyr") #needed for PCA
# install.packages("ape") #needed for PCoA
# install.packages("scales")
# install.packages("ggsn") # Add scale bars and north arrows
# install.packages("maptools") #for map addings as well
# install.packages("prettymapr")

# COMPULSORY: load the packages!!
library(isopam)
library(vegan)
library(ape)
library(dplyr)
library(scales)
library(ggsn)
library(maptools)
library(prettymapr)
library(raster)
library(ggplot2)




# Insert DATA [Bruchsal vegetation smpling] -------------------------------
# load the provided data frames see follwing URL
# (https://github.com/lukaslindenthal/data_analysis_assignment/blob/main/Bruchsal_v2.zip)

load("SPECIES/species_cover_3.RData")
load("SPECIES/species_eco_3.RData")
load("SPECIES/header_3.RData")



# A1:  classification  (isopam::isotab function) --------------------------

ip <- isopam(species_cover, c.fix = FALSE,
                           c.max = "2",
                           c.opt = FALSE,
                           l.max = 1, #for non-hierarchical partitioning
)
ip_class <- isotab(ip) #synoptical table (isopam::isotab - method)




# A3: preparation for Ordination (vegan::decorana function) -------------------------------

ordination <- decorana(species_cover)
ordination <- as.list(ordination)

#check for axis length
ordination_axis <-  c(axis1 = ordination[["evals"]][["DCA1"]],
                      axis2 = ordination[["evals"]][["DCA2"]],
                      axis3 = ordination[["evals"]][["DCA3"]],
                      axis4 = ordination[["evals"]][["DCA4"]]
                      )
ordination_axis


# A4+5: Ordination and plotting of results

source("FUNCTIONS/logmccune.R")  # logarithmic transformation for species cover
log.species_cover <- logmccune(species_cover)

# for species data a NMDS
dist <- vegdist(log.species_cover, method = "bray") #calculate diatance matrix; here Bray-Curtis distance

#function performs NMDS for 1-10 dimensions and plots dim vs. stress
NMDS.scree <- function(log.species_cover) { 
  plot(rep(1, 10), replicate(10, metaMDS(log.species_cover, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "Number of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:5) {
    points(rep(i + 1,10),replicate(10, metaMDS(log.species_cover, autotransform = F, k = i + 1)$stress))
  }
}

NMDS1 <- metaMDS(dist, k = 2,
                 trymax = 100, 
                 trace = FALSE)
plot(NMDS1)
stressplot(NMDS1)


#plot the NMDS results
col.ord <- ip_species_cover$flat  #define color as two classes
# # Define a group variable (first 53 samples belong to group 1, last 10 samples to group 2)
# group <-  c(rep("Group1", 53), rep("Group2", 11))
# # Create a vector of color values with same length as the vector of group values
# colors <-  c(rep("treat1", 53), rep("treat2", 11))

pdf("graphics/A4_ordination/NMDS_speciescover.pdf")

# plot(NMDS1, type = "n",
#      col = col.ord,
#      main = "NMDS ordination [2D]",
#      xlim = c(-1, 1), ylim = c(-1, 1),
# )
ordiplot(NMDS1, type = "n",
         main = "NMDS ordination [2D] for species data",)
ordihull(NMDS1, groups=group, 
         draw="polygon",
         col = "grey90",
         label = FALSE)
# orditorp(NMDS1, display = "sites",
#          col = col.ord,
#          cex = 1
# )
points(NMDS1,
       pch = 20,
       col = col.ord,
       cex = 1,
)
abline(v=0, h=0, lty=3)
# envfit(spe.nmds,species_eco) #add environmental variables as vectors
legend("bottomright",
       inset = c(0.02, 0.02),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 20,
       bty = "n"
)
dev.off()



# A6: Ellenverg indicator Values and Grime values
# calculation according to aggTrait function (see https://github.com/lukaslindenthal/aggTrait.git)
# I do NOT own the code for the aggTRAIT function.


#Ellenberg MOISTURE indicator
Ell_F <- aggTrait(species_cover, species_eco$Ell_F, 
                  grouping = species_eco$Species_Nr
)

log.Ell_F <- aggTrait(log.species_cover, species_eco$Ell_F,
                      grouping = species_eco$Species_Nr
                      )

#Ellenberg SOIL REACTION Zeiger
Ell_R <- aggTrait(species_cover, species_eco$Ell_R, 
                  grouping = species_eco$Species_Nr
)


#Ellenberg NITROGEN indicator
Ell_N <- aggTrait(species_cover, species_eco$Ell_N,
                  grouping = species_eco$Species_Nr
)





#same procedure for the Grime values 
Grm_S <- aggTrait(species_cover, species_eco$Grm_S,
                  grouping = species_eco$Species_Nr
                  )


Grm_C <- aggTrait(species_cover, species_eco$Grm_C,
                  grouping = species_eco$Species_Nr
                  )


Grm_R <- aggTrait(species_cover, species_eco$Grm_R,
                  grouping = species_eco$Species_Nr
                  )





# A7: scaling of the ordination plots with trait values -------------------
# install packages "scale" which includes fkt "rescale". Latest version 0.4.1

###
# now let's plot the result for "Ellenberg-moisture-indicator"
Ell_F.resc <- rescale(Ell_F$means)

log.Ell_F.resc <- rescale(Ell_F$means)

pdf("graphics/A7_EllenbergANDgrime_plots/Ordination_Ellenberg-moisture-indicator.pdf")

plot(NMDS1, type = "n",
     axes = F,
     col = col.ord,
     main = paste("NMDS ordination [2D]", "\n", "Ellenberg indicator for moisture"),
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
)
points(NMDS1, pch = 19,
       col = col.ord,
       cex = Ell_F.resc
)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 19, cex = 1, bty = "n",
)
axis(1, at=seq(-0.5, 0.5, by = 0.5))
axis(2, at=seq(-0.5, 0.5, by = 0.5))

dev.off()



###
# and for "Ellenberg-Nitrogen-indicator" as well

Ell_N.resc <- rescale(Ell_N$means)

pdf("graphics/A7_EllenbergANDgrime_plots/Ordination_Ellenberg-Nitrogen-indicator.pdf")

plot(NMDS1, type = "n",
     axes = F,
     col = col.ord,
     main = paste("NMDS ordination [2D]", "\n", "Ellenberg indicator for nitrogen"),
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
)
points(NMDS1, pch = 19,
       col = col.ord,
       cex = Ell_N.resc
)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 19, cex = 1, bty = "n",
)
axis(1, at=seq(-0.5, 0.5, by = 0.5))
axis(2, at=seq(-0.5, 0.5, by = 0.5))

dev.off()



###
# and for "Ellenberg-soil reaction-indicator" as well

Ell_R.resc <- rescale(Ell_R$means)  

pdf("graphics/A7_EllenbergANDgrime_plots/Ordination_Ellenberg-soilreaction-indicator.pdf")

plot(NMDS1, type = "n",
     axes = F,
     col = col.ord,
     main = paste("NMDS ordination [2D]", "\n", "Ellenberg indicator for soil reaction"),
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
)
points(NMDS1, pch = 19,
       col = col.ord,
       cex = Ell_R.resc
)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 19, cex = 1, bty = "n",
)
axis(1, at=seq(-0.5, 0.5, by = 0.5))
axis(2, at=seq(-0.5, 0.5, by = 0.5))

dev.off()



###
# and for "Grime C indicator" as well

Grm_C.resc <- rescale(Grm_C$means)  

pdf("graphics/A7_EllenbergANDgrime_plots/Ordination_Grime_C_strategie.pdf")

plot(NMDS1, type = "n",
     axes = F,
     col = col.ord,
     main = paste("NMDS ordination [2D]", "\n", "scaled Grime C-strategie"),
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
)
points(NMDS1, pch = 19,
       col = col.ord,
       cex = Grm_C.resc
)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 19, cex = 1, bty = "n",
)
axis(1, at=seq(-0.5, 0.5, by = 0.5))
axis(2, at=seq(-0.5, 0.5, by = 0.5))

dev.off()



###
# and for "Grime s indicator" as well

Grm_S.resc <- rescale(Grm_S$means)  

pdf("graphics/A7_EllenbergANDgrime_plots/Ordination_Grime_S_strategie.pdf")

plot(NMDS1, type = "n",
     axes = F,
     col = col.ord,
     main = paste("NMDS ordination [2D]", "\n", "scaled Grime S-strategie"),
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
)
points(NMDS1, pch = 19,
       col = col.ord,
       cex = Grm_S.resc
)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 19, cex = 1, bty = "n",
)
axis(1, at=seq(-0.5, 0.5, by = 0.5))
axis(2, at=seq(-0.5, 0.5, by = 0.5))

dev.off()



###
# and for "Grime R indicator" as well

Grm_R.resc <- rescale(Grm_R$means)  

pdf("graphics/A7_EllenbergANDgrime_plots/Ordination_Grime_R_strategie.pdf")

plot(NMDS1, type = "n",
     axes = F,
     col = col.ord,
     main = paste("NMDS ordination [2D]", "\n", "scaled Grime R-strategie"),
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
)
points(NMDS1, pch = 19,
       col = col.ord,
       cex = Grm_R.resc
)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend=c("class 1", "class 2"), col = col.ord,
       pch = 19, cex = 1, bty = "n",
)
axis(1, at=seq(-0.5, 0.5, by = 0.5))
axis(2, at=seq(-0.5, 0.5, by = 0.5))

dev.off()




# A8: scaling of the ordination plots with head data ----------------------

###
# now let's plot the result for header: living plants
lp.header <- rescale(header$Lebende_Pflanzen)
pdf("graphics/A8_headdata/Ordination_livingplants.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "living plants scaled"),
     col = col.ord
     )
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = lp.header
       )
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()

###
# and now for header: bare ground
bg.header <- rescale(header$Offener_Untergrund)
pdf("graphics/A8_headdata/Ordination_bareground.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by bare ground"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = bg.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()




###
# and now for header: litter
lit.header <- rescale(header$Streu)
pdf("graphics/A8_headdata/Ordination_litter.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by litter"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = lit.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()


###
# and now for header: deadwood
dw.header <- rescale(header$Totholz_liegend)
pdf("graphics/A8_headdata/Ordination_deadwood.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by litter"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = dw.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()


###
# and now for header: tree cover
covtree.header <- rescale(header$Cover_Trees)
pdf("graphics/A8_headdata/Ordination_treecover.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by tree cover"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = covtree.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()


###
# and now for header: shrubs cover
covshrub.header <- rescale(header$Cover_Shrubs)
pdf("graphics/A8_headdata/Ordination_shrubcover.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by shrub cover"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = covshrub.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()


###
# and now for header: vascular plants
vascplant.header <- rescale(header$Cover_Field_Layer_vascular)
pdf("graphics/A8_headdata/Ordination_vascular.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by vascular plants"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = vascplant.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()


###
# and now for header: NON vascular plants
NONvascplant.header <- rescale(header$Cover_Field_Layer_non_vascular)
pdf("graphics/A8_headdata/Ordination_nonvascular.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by non vascular plants"),
     col = col.ord
)
points(NMDS1,
       col = col.ord,
       pch = 19,
       cex = NONvascplant.header
)
abline(v=0, h=0, lty = 3)
legend("bottomright",
       inset = c(0.03, 0.03),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n",
       cex = 1)

dev.off()


###
# and now for header: exposition
ex.header <- rescale(header$Exposition)
sl.header <- rescale(abs(header$Hangneigung))
sl.header2 <- rescale(header$Hangneigung)

Ex.radians <- ((header$Exposition/360)*2*pi)

pdf("graphics/A8_headdata/Ordination_exposition.pdf")

plot(NMDS1, type = "n",
     main = paste("Ordination [2D]", "\n", "scaled by exposition"),
     xlim = c(-1 , 1), ylim = c(-1, 1)
)
arrows(NMDS1$points[1:64, 1], NMDS1$points[1:64, 2],
       NMDS1$points[1:64, 1] + sin(Ex.radians) * sl.header/2, spe.nmds$points[1:64, 2] + cos(Ex.radians) * sl.header/2,
       col = col.ord, 
       angle = 15,
       length = 0.1
)
points(sl.header, 
       cex = sl.header,
       col = col.ord,
       pch = 19)
abline(v = 0, h = 0, lty = 3)
legend("bottomright",
       inset = c(0.06, 0.04),
       legend = c("class 1", "class 2"),
       # pch = c(NA, NA),
       lwd = 3, col = c(1, 2),
       lty = c(NA, NA),
       bty = "n") #normal legend. Do not plot line
par(font = 5) #different font displays arrows
legend("bottomright",
       inset = c(0.14, 0.04),
       legend = c(NA, NA), pch = c(174, 174),
       lwd = 3, col = c(1, 2),
       lty = c(NA, NA), 
       bty = "n")
par(font = 1) #and we're back to default again;)

dev.off()



# A9: create geographical map of the given data --------------------------

DEM_raster <- raster("SPATIAL/dem.tif")
DEM_aspect <- raster("SPATIAL/aspect.tif")
DEM_hillshade <- raster("SPATIAL/hillshade.tif")
DEM_slope <- raster("SPATIAL/slope.tif")

# DEM_raster <- aggregate(DEM_raster, fact=2) ##check dimensions and if running slow resample to 1/2 resolution
#plot for DEM_raster
df.DEM_raster <- as.data.frame(DEM_raster, xy = TRUE)

pdf("graphics/A9_DEM_headdata/terrain_exposition.pdf")

ggplot()+
  geom_tile(data = df.DEM_raster, aes(x=x, y=y, fill = dem))+ # add the provided data as fill
  scale_fill_gradientn(colors = terrain.colors(10))+
  labs(title = "DEM of study area",
       subtitle = "scaled bay terrain attribute exposition",
       caption = "Data source: course material")+
  theme_bw()
dev.off()


###
#terrrain attribute: slope 
slp <-  terrain(DEM_slope, opt = "slope", unit = "degrees")
df.DEM_slope <- as.data.frame(slp, xy = TRUE)

pdf("graphics/A9_DEM_headdata/terrain_slope.pdf")

ggplot()+
  geom_tile(data = df.DEM_slope, aes(x=x, y=y, fill = slope))+ # same procedure for slope
  scale_fill_gradientn(colors = topo.colors(10))+
  labs(title = "DEM of study area",
       subtitle = "scaled bay terrain attribute slope",
       caption = "Data source: course material")+
  theme_bw()
dev.off()

###
#terrain attribute: aspect
asp <- terrain(DEM_aspect, opt = "aspect", unit = "degrees")
df.DEM_aspect <- as.data.frame(asp, xy = TRUE)

pdf("graphics/A9_DEM_headdata/terrain_aspect.pdf")

ggplot()+
  geom_tile(data = df.DEM_aspect, aes(x=x, y=y, fill = aspect))+# as well as aspect
  scale_fill_gradientn(colors = rainbow(25))+
  labs(title = "DEM of study area",
       subtitle = "scaled bay terrain attribute aspect",
       caption = "Data source: course material")+
  theme_bw()
dev.off()


###
# DEM with gray.colors and living plants
pdf("graphics/A9_DEM_headdata/DEM_ordination_livingplants.pdf")

# plot(DEM_raster,  
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by living plants"))
points(header[, 2:1], pch = 19, col = col.ord, cex = lp.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()


###
# DEM with gray.colors and dead wood
pdf("graphics/A9_DEM_headdata/DEM_ordination_deadwood.pdf")

# plot(DEM_raster,  
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by dead wood"))
points(header[, 2:1], pch = 19, col = col.ord, cex = dw.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()


###
# DEM with gray.colors and litter
pdf("graphics/A9_DEM_headdata/DEM_ordination_litter.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by litter"))
points(header[, 2:1], pch = 19, col = col.ord, cex = lit.header
       )
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()


###
# DEM with gray.colors and bare ground
pdf("graphics/A9_DEM_headdata/DEM_ordination_bareground.pdf")
# plot(DEM_raster,  
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by bare ground"))
points(header[, 2:1], pch = 19, col = col.ord, cex = bg.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()



###
# DEM with gray.colors and tree cover
pdf("graphics/A9_DEM_headdata/DEM_ordination_treecover.pdf")

# plot(DEM_raster,  
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by tree cover"))
points(header[, 2:1], pch = 19, col = col.ord, cex = covtree.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()



###
# DEM with gray.colors and shrub cover
pdf("graphics/A9_DEM_headdata/DEM_ordination_shrubcover.pdf")

# plot(DEM_raster,  
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by shrub cover"))
points(header[, 2:1], pch = 19, col = col.ord, cex = covshrub.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()



## # DEM with gray.colors and non vascular plants
pdf("graphics/A9_DEM_headdata/DEM_ordination_nonvascular.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Classification by geographical area", "\n", "scaled by non vascular plants"))
points(header[, 2:1], pch = 19, col = col.ord, cex = NONvascplant.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black")
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()


###
# A11: PCQM function calculating head data and plotting the results -------

# PCQM-function and map display  ------------------------------------------
#for function see https://github.com/lukaslindenthal/PCQM-function.git
# function provided by the course management

# PCQM function
# Requires a list with matrices, each with 4 rows and columns for species name, distance and DBH in cm (at 130 cm)

pcqm(dat)

load("PCQM/pcqm_2021_v2.RData")
load("PCQM/pcqm_2021_xy_v2.RData")

z1 <- data.frame(Species = c("Carpinus", "Fagus", "Quercus", "Acer"), Distance = c(8, 7, 2, 12), DBH = c(13.5, 23, 74, 23.3))
z2 <- data.frame(Species = c("Fagus", "Fagus", "Fagus", "Fagus"), Distance = c(4, 3, 8, 2), DBH = c(4.5, 27, 6, 4))
z3 <- data.frame(Species = c("Fagus", "Fagus", "Carpinus", "Fagus"), Distance = c(3, 5, 6, 6), DBH = c(14, 15, 10, 10))
z4 <- data.frame(Species = c("Acer", "Fagus", "Fagus", "Fagus"), Distance = c(1, 5, 5, 9), DBH = c(6, 66, 53, 48))
z5 <- data.frame(Species = c("Acer", "Acer", "Fraxinus", "Fagus"), Distance = c(3, 5, 3, 13), DBH = c(12, 9, 4, 60))

dat <- list(z1, z2, z3, z4, z5)

pcqm <- function(dat) {
  result <- list()
  # Nearest Point-to-tree-distances
  # The sum of the nearest point-to-tree distances in the quarters surveyed divided by the number of quarters
  r <- sum(unlist(lapply(dat, function(i) sum(i$Distance)))) / (length(dat) * 4)
  result[["Mean_distance"]] <- r
  
  # Absolute density
  d <- 1 / (r^2)
  result[["Absolute_density"]] <- d
  
  # Absolute density of each species
  species <- mapply('[', dat, TRUE, 1)
  dk <- as.matrix(table(species) / (length(dat) * 4) * d * 10000)
  dkn <- as.numeric(dk)
  names(dkn) <- rownames(dk)
  result[["Absolute_densities_of_species"]] <- dkn
  
  # Relative density of a species
  rd <- dk / 10000 / d * 100
  rdn <- as.numeric(rd)
  names(rdn) <- rownames(rd)
  result[["Relative_densities_of_species"]] <- rdn
  
  # Basal area of each tree
  dbh <- mapply('[', dat, TRUE, 3)
  # c2/4Ï
  ba <- as.numeric(dbh) ^ 2 / 4 * pi
  mba <- aggregate(ba, by = list(as.factor(species)), FUN = mean)
  mban <- as.numeric(mba[, 2])
  names(mban) <- mba[, 1]  
  result[["Mean_basal_area_of_species"]] <- mban
  
  # Basal_area_per_ha
  tba <- mba[, 2] * dk / 10000
  colnames(tba) <- "Basal_area_per_ha"
  tban <- as.numeric(tba)
  names(tban) <- rownames(tba)
  result[["Basal_area_per_ha_and_species"]] <- tban
  
  # Basal_area_per_ha
  sba <- sum(tba)
  result[["Total_basal_area_per_ha"]] <- sba
  
  # Relative Cover (Relative Dominance) of a Species
  rc <- tban / sba * 100
  result[["Relative_cover_by_species"]] <- rc
  
  # Absolute Cover of a Species
  io <- nocc <- table(species, rep(1:length(dat), each = 4))
  io[io > 0] <- 1
  # Number of sample points with a species
  nsp <- rowSums(io)
  result[["Number_of_sample_points_with_species"]] <- nsp
  
  # Frequency of species
  fps <- nsp / length(dat) * 100
  result[["Frequency_by_species"]] <- fps
  
  # Relative frequency of species
  rfps <- fps / sum(fps) * 100
  result[["Relative_frequency_per_species"]] <- rfps
  
  # Theimportance valueof a species is defined as the sum of the three relative measures
  # Relative density+Relative cover+Relative frequency
  ai <- rdn + rc + rfps
  ri <- ai / sum(ai) * 100
  result[["Relative_importance_of_species"]] <- ri
  
  result
}  

pcqm(dat)

pcqm_dat <- pcqm(dat)


### plotting of the pcqm data 
## # DEM with gray.colors, pcqm and living plants
pdf("graphics/A11_PCQM//PCQM_ordination_livingplants.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by living plants")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = lp.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()



###
# DEM with gray.colors, pcqm and bare ground
pdf("graphics/A11_PCQM//PCQM_ordination_bareground.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by bare ground")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = bg.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()



###
# DEM with gray.colors, pcqm and litter
pdf("graphics/A11_PCQM//PCQM_ordination_litter.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by litter")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = lit.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()




###
# DEM with gray.colors, pcqm and dead wood
pdf("graphics/A11_PCQM//PCQM_ordination_deadwood.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by dead wood")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = dw.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()




###
# DEM with gray.colors, pcqm and tree cover
pdf("graphics/A11_PCQM//PCQM_ordination_treecover.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by tree cover")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = covtree.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()




###
# DEM with gray.colors, pcqm and shrub cover
pdf("graphics/A11_PCQM//PCQM_ordination_shrubcover.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by shrub cover")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = covshrub.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()




###
# DEM with gray.colors, pcqm and vascular plants
pdf("graphics/A11_PCQM//PCQM_ordination_vascular.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by vascular plants")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = vascplant.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()




###
# A10: MaxEnt for selected species  ---------------------------------------

load("SPECIES/species_cover_3.RData")

# for species of class 1 we need to convert the data into a .csv format

#for species Galium odoratum
Gal_odoratum <- data.frame(species_cover[, c("Galium odoratum hl")])

Gal_odoratum$species <- Gal_odoratum$Galium.odoratum.hl
MaxEnt_Gal_odoratum <- data.frame(Gal_odoratum$species, header$x , header$y)
MaxEnt_Gal_odoratum[MaxEnt_Gal_odoratum == 0] <- NA
MaxEnt_Gal_odoratum <- na.omit(MaxEnt_Gal_odoratum)
MaxEnt_Gal_odoratum$Gal_odoratum.species <- "Galium odoratum"
write.csv(MaxEnt_Gal_odoratum,"C:/Dokumente und Einstellungen/anon/Documents/R/RdataAnalysisWS21/MaxEnt in R/Galiumodoratum/Galiumodoratum.csv", row.names = FALSE)


###
#now for species of class 2

# Quercus petraea
Querc_pet <- data.frame(species_cover[, c("Quercus petraea sl",
                                          "Quercus petraea t")])

Querc_pet$species <- Querc_pet$Quercus.petraea.sl + Querc_pet$Quercus.petraea.t
MaxEnt_Querc_pet <- data.frame(Querc_pet$species, header$x , header$y)
MaxEnt_Querc_pet[MaxEnt_Querc_pet == 0] <- NA
MaxEnt_Querc_pet <- na.omit(MaxEnt_Querc_pet)
MaxEnt_Querc_pet$Querc_pet.species <- "Querc_petraea"
write.csv(MaxEnt_Querc_pet,"C:/Dokumente und Einstellungen/anon/Documents/R/RdataAnalysisWS21/MaxEnt in R/Quercuspetraea/Quercuspetraea.csv", row.names = FALSE)







###
# DEM with gray.colors, pcqm and non vascular plants
pdf("graphics/A11_PCQM//PCQM_ordination_nonvascular.pdf")
# plot(DEM_raster,
#      col = gray.colors(20), legend = FALSE,
#      # ylim = c(49.095, 49.110), xlim = c(8.565, 8.590)
# )
plot(DEM_hillshade, col = gray.colors(10), legend = FALSE,
     main = paste("Ordination [2D]", "\n", "PCQM scaled by non vascular plants")
)
points(header[, 2:1], pch = 19, col = col.ord, cex = NONvascplant.header)
legend("bottomright",
       inset = c(0.01, 0.01),
       legend = c("class1", "class2"),
       pch = 19,
       col = col.ord,
       bty = "n"
)
addnortharrow("topright", 
              padin = c(0.55, 0.15), scale = 0.5,
              lwd = 0.5, border = "black", cols = c("white", "black"),
              text.col = "black"
)
addscalebar(#model = "WGS84", #auto assuming lat/lon (epsg4326 = Gauß-Krüger)
  # dist = 10,
  # st.size = 1,
  # dist_unit = "km",
  # transform = TRUE
)
dev.off()






# BYPASS: create tables for LaTex -----------------------------------------
# install.packages("stargazer")
library(stargazer)
library(tables)

#synoptical table (isopam::isotab - method)
ip_table <- isotab(ip)
x <- as.data.frame(ip_table)
xy <- x[1:8, 1:5]

stargazer(xy, title = "output isopam classsificaton and diagnostic species",
          style = "aer",
          decimal.mark = "."
          )


#table for ordination; axis length and EV
stargazer(decorana(species_cover), title = "output isopam classsificaton and diagnostic species",
          style = "default",
          decimal.mark = "."
)




# ADD ON: map of study area -----------------------------------------------


# DEM_raster <- aggregate(DEM_raster, fact=2) ##check dimensions and if running slow resample to 1/2 resolution
#plot for DEM_raster
df.DEM_raster <- as.data.frame(DEM_raster, xy = TRUE)
ggplot()+
  geom_tile(data = df.DEM_raster, aes(x=x, y=y, fill = dem))+
  scale_fill_gradientn(colors = terrain.colors(10))+
  theme_bw()



#terrrain attribute: slope 
slp <-  terrain(DEM_slope, opt = "slope", unit = "degrees")
df.DEM_slope <- as.data.frame(slp, xy = TRUE)

ggplot()+
  geom_tile(data = df.DEM_slope, aes(x=x, y=y, fill = slope))+
  scale_fill_gradientn(colors = topo.colors(10))+
  labs(title = "DEM of study area",
       subtitle = "scaled bay terrain attribute slope",
       caption = "Data source: course material")
  theme_bw()

#create GeoTiff file for GIS
write(slp, "Bruchsal_terrainslope.tif")


#terrain attribute: aspect
asp <- terrain(DEM_aspect, opt = "aspect", unit = "degrees")
df.DEM_aspect <- as.data.frame(asp, xy = TRUE)

ggplot()+
  geom_tile(data = df.DEM_aspect, aes(x=x, y=y, fill = aspect))+
  scale_fill_gradientn(colors = rainbow(10))+
  theme_bw()
