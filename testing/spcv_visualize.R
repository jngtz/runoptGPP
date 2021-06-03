library(ggplot2)
library(sf)
library(patchwork)
library(rgeos)
library(raster)
library(rgdal)

# Spatial Data for Plots #############

setwd("/home/jason/Data/Chile/")
rivers <- st_read("clip_rios_principales.shp")
basins <- st_read("sub_catchments.shp")
boundary <- st_read("study_Area_bnd.shp")
water <- st_read("landuse_map_water.shp")


slide_poly_vec <- readOGR(".", "dflow_polygons_v1_reposition_sample_100")
slide_poly_vec$objectid <- 1:100

setwd("saga_data")
hillshade <- raster("hs_filtered_filled.sdat")
hillshade <- aggregate(hillshade, fact = 8, fun = mean)
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
hillshade_df <- hillshade_df[!is.na(hillshade_df$hs_filtered_filled),]
#https://gist.github.com/dirkseidensticker/ce98c6adfe16d5e4590e95c587ea0432
#https://philmikejones.me/tutorials/2015-09-03-dissolve-polygons-in-r/

carea <- raster("carea_alos_12_5m_filled.sdat")

carea_km2 <- carea * 1e-06

#friction_mu_grd <- 0.25*carea_km2^-0.21 # min runout likely
#friction_mu_grd <- 0.19*carea_km2^-0.24 # most likely
friction_mu_grd <- 0.13*carea_km2^-0.25 # max runout

friction_mu_grd[friction_mu_grd > 0.3] <- 0.3
friction_mu_grd[friction_mu_grd < 0.045] <- 0.045

setwd("/home/jason/Scratch/GPP_RW_Paper")

slides <- st_read("dflow_polygons_v1_reposition_sample_100_CV.shp")


ds <- gCentroid(as(slides, "Spatial"), byid=TRUE, id = slides$objectid_2)


d("/home/jason/Scratch/Figures")


# Spatial Folds ######################
setwd("/home/jason/Scratch/GPP_RW_Paper")

slides <- st_read("dflow_polygons_v1_reposition_sample_100_CV.shp")

slides$fold1 <- "Training data"
slides$fold2 <- "Training data"
slides$fold3 <- "Training data"
slides$fold4 <- "Training data"
slides$fold5 <- "Training data"

slides$fold1[slides$spcv_label == 1] <- "Test data"
slides$fold2[slides$spcv_label == 2] <- "Test data"
slides$fold3[slides$spcv_label == 3] <- "Test data"
slides$fold4[slides$spcv_label == 4] <- "Test data"
slides$fold5[slides$spcv_label == 5] <- "Test data"

slides$fold1 <- as.factor(slides$fold1)
slides$fold2 <- as.factor(slides$fold2)
slides$fold3 <- as.factor(slides$fold3)
slides$fold4 <- as.factor(slides$fold4)
slides$fold5 <- as.factor(slides$fold5)


m.f1 <-

  ggplot() +
  geom_sf(data = boundary, alpha = 0.8, colour = NA, fill = "#F8F9F9") +
  geom_sf(data = rivers, colour = "#85C1E9", lwd = 0.3) +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9", lwd = 0.2) +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA, lwd = .2) +
  geom_sf(data = slides, lwd = 1, aes(colour = fold1, fill = fold1),  show.legend = FALSE) +
  scale_fill_manual(values = c("#7DCEA0", "#283747")) +
  scale_colour_manual(values = c("#7DCEA0", "#283747")) +
  xlab("Fold 1") +
  ylab("") +
  coord_sf(datum = NA) +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 10),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10, hjust = 0))


m.f2 <-

  ggplot() +
  geom_sf(data = boundary, alpha = 0.8, colour = NA, fill = "#F8F9F9") +
  geom_sf(data = rivers, colour = "#85C1E9", lwd = 0.3) +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9", lwd = 0.2) +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA, lwd = .2) +
  geom_sf(data = slides, lwd = 1, aes(colour = fold2, fill = fold2),  show.legend = FALSE) +
  scale_fill_manual(values = c("#7DCEA0", "#283747")) +
  scale_colour_manual(values = c("#7DCEA0", "#283747")) +
  xlab("Fold 2") +
  ylab("") +
  coord_sf(datum = NA) +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 10),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10, hjust = 0))


m.f3 <-

  ggplot() +
  geom_sf(data = boundary, alpha = 0.8, colour = NA, fill = "#F8F9F9") +
  geom_sf(data = rivers, colour = "#85C1E9", lwd = 0.3) +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9", lwd = 0.2) +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA, lwd = .2) +
  geom_sf(data = slides, lwd = 1, aes(colour = fold3, fill = fold3),  show.legend = FALSE) +
  scale_fill_manual(values = c("#7DCEA0", "#283747")) +
  scale_colour_manual(values = c("#7DCEA0", "#283747")) +
  xlab("Fold 3") +
  ylab("") +
  coord_sf(datum = NA) +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 10),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 10, hjust = 0))



m.f4 <-

  ggplot() +
  geom_sf(data = boundary, alpha = 0.8, colour = NA, fill = "#F8F9F9") +
  geom_sf(data = rivers, colour = "#85C1E9", lwd = 0.3) +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9", lwd = 0.2) +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA, lwd = .2) +
  geom_sf(data = slides, lwd = 1, aes(colour = fold4, fill = fold4),  show.legend = FALSE) +
  scale_fill_manual(values = c("#7DCEA0", "#283747")) +
  scale_colour_manual(values = c("#7DCEA0", "#283747")) +
  #scale_fill_manual(name = "Debris flow\npolygons", values = c("#00B050", "#323232")) +
  #scale_colour_manual(name = "Debris flow\npolygons",values = c("#00B050", "#323232")) +


  xlab("Fold 4") +
  ylab("") +
  coord_sf(datum = NA) +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 9),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 9, hjust = 0),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.5, "cm"))


m.f5 <-

  ggplot() +
  geom_sf(data = boundary, alpha = 0.8, colour = NA, fill = "#F8F9F9") +
  geom_sf(data = rivers, colour = "#85C1E9", lwd = 0.3) +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9", lwd = 0.2) +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA, lwd = .2) +
  geom_sf(data = slides, lwd = 1, aes(colour = fold5, fill = fold5),  show.legend = FALSE) +

  scale_fill_manual(values = c("#7DCEA0", "#283747")) +
  scale_colour_manual(values = c("#7DCEA0", "#283747")) +

  #scale_fill_manual(name = "Debris flow\npolygons", values = c("#e3636c", "#323232")) +
  #scale_colour_manual(name = "Debris flow\npolygons",values = c("#e3636c", "#323232")) +
  xlab("Fold 5") +
  ylab("") +
  coord_sf(datum = NA) +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 9),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 9, hjust = 0),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"))





m.f1 + m.f2 + m.f3 + m.f4 + m.f5 + plot_layout(ncol = 5)

setwd("/home/jason/Scratch/Figures")

ggsave("spatial_cv_folds.png", dpi = 300, width = 7, height = 2.5, units = "in")



# Map PCM Opt Individual slides ######################

setwd("/home/jason/Scratch/GPP_PCM_Paper")

(load("gridsearch_pcm_settings.Rd"))

polyid.vec <- 1:pcm_settings$n_train
pcmmu.vec <- pcm_settings$vec_pcmmu
pcmmd.vec <- pcm_settings$vec_pcmmd

opt_list <- list() # list of optimal parameters per landslides

for(i in 1:pcm_settings$n_train){


  res_nm <- paste("result_relerr_length_", i, ".Rd", sep="")
  load(res_nm) #res
  res <- relerr_length_result

  # Get optimal paramters per landslide

  wh_opt <- which(res==min(res), arr.ind = TRUE)
  mu_opt <- pcmmu.vec[wh_opt[,1]]
  md_opt <- pcmmd.vec[wh_opt[,2]]

  df_opt <- data.frame(mu = mu_opt, md = md_opt, relerr = min(res))

  opt_list[[i]] <- df_opt
  # calc median for these using apply

}


#In case two optimal params found, take first
sel_opt <- function(x){
  if(nrow(x) > 1){
    opt_comb <- x[1,]
  } else {
    opt_comb <- x
  }

  return(opt_comb)
}


opt_sel_list <- lapply(opt_list, sel_opt)

df_sol <- data.frame(slide_id = 1:pcm_settings$n_train)
df_sol <- cbind(df_sol, do.call("rbind", opt_sel_list))

# summarize optimal parameters across landslides

df_pcm <- data.frame(mu = df_sol$mu, md = df_sol$md)

# Number of optimal solutions/slide

freq_pairs <- table(df_pcm)
freq_pairs != 0

mu_nm <- as.numeric(rownames(freq_pairs))
md_nm <- as.numeric(colnames(freq_pairs))

# Get array index of pairs
pair_ind <- which(freq_pairs !=0, arr.ind = TRUE)
pairs <- data.frame(mu = mu_nm[pair_ind[,1]],
                    md = md_nm[pair_ind[,2]],
                    freq = freq_pairs[pair_ind])
pairs$rel_freq <- pairs$freq/nrow(df_sol) * 100

# find relative error
pair_id <- paste(pairs$mu, pairs$md)
df_sol$pair_id <- paste(df_sol$mu, df_sol$md)

for(i in 1:length(pair_id)){

  relerr_i <- df_sol$relerr[which(df_sol$pair_id == pair_id[i])]
  pairs$rel_err[i] <- median(relerr_i)
  pairs$iqr_relerr[i] <- IQR(relerr_i, na.rm = TRUE)

}



mybreaks <- c(1, 2, 3)

ggplot(pairs, aes(x=md, y=mu)) +
  geom_point(alpha=0.7, aes(colour = rel_err, size = rel_freq)) +


  scale_size(name="Relative\nfrequency (%)", range = c(2, 6), breaks = mybreaks) +

  scale_colour_gradient(high = "#1B4F72", low = "#85C1E9", name = "Median relative\nerror") +

  xlab("Mass-to-drag ratio (m)") +
  ylab("Sliding friction coefficient") +


  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 9),
        axis.text = element_text(size = 8)) +
  guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2))

setwd("/home/jason/Scratch/Figures")
ggsave("individual_scatter_opt_pcm.png", dpi = 300, width = 5.5, height = 3.25, units = "in")


# Number of solutions #####################################
df_sol$n_solutions <- sapply(opt_list, nrow)


h.n_sol <-
  ggplot(df_sol, aes(n_solutions)) +
  geom_histogram(colour = "black", binwidth = 0.1) +
  xlab("No. of distance model solutions") +
  ylab("Count") +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_blank())

h.n_sol

setwd("/home/jason/Scratch/Figures")

ggsave("hist_num_solutions_pcm.png", dpi = 300, width = 3.5, height = 2, units = "in")


# Map of no. optimal solutions


df_sol$x <- coordinates(ds)[,1]
df_sol$y <- coordinates(ds)[,2]

mybreaks <- c(1, 10, 25, 50, 100, 150)

m.n_sol <-
  ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +
  geom_sf(data = boundary, alpha = 0.8, fill = "white") +
  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
#  geom_point(data = df_sol, size = 2.5, alpha=0.8, aes(x = x, y = y, colour = n_solutions)) +
  geom_point(data = df_sol, alpha=0.8, aes(x = x, y = y, size = n_solutions, colour = n_solutions)) +
  geom_point(data = df_sol, shape = 1, fill = NA, alpha = 0.9, colour="#7f8c8d", aes(x = x, y = y, size = n_solutions)) +

  scale_size_continuous(name="No. of\nsolutions", range=c(2,5), breaks=mybreaks) +
  scale_alpha_continuous(name="No. of\nsolutions", range=c(0.1, .9), breaks=mybreaks) +
  scale_colour_gradient(name = "No. of\nsolutions", low = "#fef9e7", high = "#512e5f", breaks = mybreaks) +

  guides( colour = guide_legend()) +
  xlab("") +
  ylab("") +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 10), axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),legend.position = "right")




m.pcm_sol <-
  ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +
  geom_sf(data = boundary, alpha = 0.8, fill = "white") +
  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
  geom_point(data = df_sol, alpha=0.8, aes(x = x, y = y, size = md, colour = mu)) +
  scale_colour_gradient(low = "#fef9e7", high = "#117864",
                        name = "Sliding friction\ncoefficient") +
  geom_point(data = df_sol, shape = 1, fill = NA, alpha = 0.9, colour="#7f8c8d", aes(x = x, y = y, size = md)) +
  scale_size(name="Mass-to-drage\nratio (m)", range = c(1, 3), breaks = c(20, 75, 150 )) +
  xlab("") +
  ylab("") +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),legend.position = "right")


m.n_sol + m.pcm_sol

ggsave("map_nsol_bottom_fnt10.png", dpi = 300, width = 7.5, height = 6.7, units = "in")

# MAP Individual opt error ###################
m.relerr <-
  ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +
  geom_sf(data = boundary, alpha = 0.8, fill = "white") +
  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
  geom_point(data = df_sol, size = 2, alpha=0.8, aes(x = x, y = y, colour = relerr)) +
  scale_colour_gradient(low = "#D0ECE7", high = "#0B5345",
                        name = "Runout distance\nrelative error") +
  geom_point(data = df_sol, shape = 1, size = 2, alpha = 0.9, colour="#7f8c8d", aes(x = x, y = y)) +
  xlab("") +
  ylab("") +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 8),
        axis.text = element_text(size = 7), legend.position = "right")


m.pcm_sol + m.relerr
ggsave("map_indv_relerr_sol_bottom_fnt9.png", dpi = 300, width = 7.5, height = 4, units = "in")

# Possible outliers in data or model cannot run #########################
boxplot(df_sol$relerr)
df_sol[df_sol$relerr >= 0.3,]


# Summary of spatial variation in mu from sample ###################

med_mu <- raster::extract(friction_mu_grd, slide_poly_vec, fun = median)
iqr_mu <- raster::extract(friction_mu_grd, slide_poly_vec, fun = IQR)

length(med_mu[med_mu > 0.21])
hist(med_mu) # most between 0.21 and 0.3
hist(iqr_mu) # not much variation in mu through the slide...

#save(df_sol, file = "indv_opt_result.Rd")

# Correlation to PCM paramters #########################################
setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("indv_opt_result.Rd"))

spdf_sol <- df_sol
coordinates(spdf_sol) <- ~ x + y

library(raster)
setwd("/home/jason/Data/Chile/saga_data")
dem <- raster("dem_alos_12_5m_fill.sdat")
slope <- raster("slope_alos_12_5m_filled.sdat")
carea <- raster("carea_alos_12_5m_filled.sdat")

setwd("/home/jason/Data/Chile/")
landc <- raster("landuse_sub.tif")

layers <- stack(dem, slope, carea, landc)
names(layers) <- c('elev', 'slope', 'carea', 'landuse')
extr <- raster::extract(layers, spdf_sol, df = TRUE)

d <- cbind(spdf_sol, extr)

d$landc <- NA
d$landc[d$landuse == 2] <- "No Vegetation"
d$landc[d$landuse == 3] <- NA #"UrbInd"
d$landc[d$landuse == 4] <- "Vegetation"
d$landc[d$landuse == 5] <- NA #"Watr"
d$landc[d$landuse == 6] <- "Wetland"
d$landc[d$landuse == 7] <- "No Vegetation"#"SnwGl"
d$landc[d$landuse == 8] <- "Vegetation"
d$landc[d$landuse == 9] <- "Agriculture"
d$landc[d$landuse == 10] <- NA #No veg?

d$landc <- as.factor(d$landc)

cor(d$elev, d$mu, method = "spearman")
cor(d$slope, d$mu, method = "spearman")
cor(log10(d$carea), d$mu, method = "spearman")

cor(d$elev, d$md, method = "spearman")
cor(d$slope, d$md, method = "spearman")
cor(log10(d$carea), d$md, method = "spearman")


###

setwd("/home/jason/Scratch/Figures")
#png(filename="bxplt_indv_opt_landcover.png", res = 300, width = 7.5, height = 3,
 #   units = "in", pointsize = 11)

par(family = "Arial", mfrow = c(1,2))

#, mar = c(4, 3, 0.5, 0.5),
 #   mgp = c(2, 0.75, 0))

boxplot(d$mu ~ d$landc, varwidth = TRUE, xlab = "Land cover",
        ylab = "Sliding friction coefficent", cex.axis = 1, cex.lab = 1)

boxplot(d$md ~ d$landc, varwidth = TRUE, xlab = "Land cover",
        ylab = "Mass-to-drag ratio (m)", cex.axis = 1, cex.lab = 1)

#dev.off()

library(ggplot2)
library(patchwork)

b.mu <- ggplot(data = d@data, aes(x = landc, y = mu)) +
  geom_boxplot(fill='gray', color="black", varwidth = TRUE) +
  xlab("Land cover") +
  ylab("Sliding friction coefficient") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 10),
        axis.text = element_text(size = 10), legend.position = "bottom")

b.md <- ggplot(data = d@data, aes(x = landc, y = md)) +
  geom_boxplot(fill='gray', color="black", varwidth = TRUE) +
  xlab("Land cover") +
  ylab("Mass-to-drag ratio (m)") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 10),
        axis.text = element_text(size = 10), legend.position = "bottom")

patch <- b.mu + b.md
patch

ggsave("boxplot_ind_optPCM_landc_veg.png", dpi = 300, width = 7.5, height = 3, units = "in")
