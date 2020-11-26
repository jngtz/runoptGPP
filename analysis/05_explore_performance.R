# LOAD PACKAGES ################################################################

library(runoptGPP)
library(rgdal)
library(raster)
library(rgeos)
library(sf)
library(ggplot2)
library(patchwork) # for multiple plots with ggplot

library(Rsagacmd)

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# LOAD DATA ####################################################################

setwd("/home/jason/Data/Chile/")
# elevation model
dem <- raster("dem_alos_12_5m _no sinks.tif")

# slide start/source points
slide_point_vec <- readOGR(".", "dflow_points_v1_reposition")

# actual/mapped debris flow polygons
slide_poly_vec <- readOGR(".", "dflow_polygons_v1_reposition_sample_100")
slide_poly_vec$objectid <- 1:100

crs(slide_point_vec) <- crs(slide_poly_vec)

rivers <- st_read("clip_rios_principales.shp")
basins <- st_read("sub_catchments.shp")
boundary <- st_read("study_Area_bnd.shp")
water <- st_read("landuse_map_water.shp")

setwd("saga_data")
hillshade <- raster("hs_filtered_filled.sdat")
hillshade <- aggregate(hillshade, fact = 6, fun = mean)
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
hillshade_df <- hillshade_df[!is.na(hillshade_df$hs_filtered_filled),]

setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("rw_gridsearch_multi.Rd"))

setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("pcm_gridsearch_multi.Rd"))

# Get grid search space
pcm_md_vec <- as.numeric(colnames(pcm_gridsearch_multi[[1]][[1]]))
pcm_mu_vec <- as.numeric(rownames(pcm_gridsearch_multi[[1]][[1]]))

rwslp_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[1]])
rwexp_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[2]])
rwper_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[3]])

# PERFORM SPATIAL CROSS VALIDATION #############################################

#RW SPCV
setwd("/home/jason/Scratch/GPP_RW_Paper")
rw_spcv <- rwSPCV(x = rw_gridsearch_multi, slide_plys = slide_poly_vec,
                  n_folds = 5, repetitions = 1000)

save("rw_spcv", file = "rw_1000rep_5fld_spcv.Rd" )

#PCM SPCV
setwd("/home/jason/Scratch/GPP_PCM_Paper")
pcm_spcv <- pcmSPCV(pcm_gridsearch_multi, slide_plys = slide_poly_vec,
                    n_folds = 5, repetitions = 1000, from_save = FALSE)

save(pcm_spcv, file = "pcm_1000rep_5fld_spcv.Rd")

# Plot spatial cross validation results ########################################

setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("rw_1000rep_5fld_spcv.Rd"))
freq_rw <- rwPoolSPCV(rw_spcv, plot_freq = TRUE)

setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("pcm_1000rep_5fld_spcv.Rd"))
freq_pcm <- pcmPoolSPCV(pcm_spcv, plot_freq = TRUE)

# RW
p.rw <- ggplot(freq_rw, aes(x=per, y=exp)) +
  geom_point(alpha=0.7, aes(colour = median_auroc, size = rel_freq)) +
  scale_size(range = c(0.5, 3), name="Relative\nfrequency (%)",
             breaks = c(5, 15, 60)) +
  scale_colour_gradient(low = "#1B4F72", high = "#85C1E9",
                        name = "Median\nAUROC") +
  scale_x_continuous(expression(paste("Persistence factor")),
                     limits = c(min(rwper_vec), max = max(rwper_vec))) +
  scale_y_continuous(expression(paste("Exponent of divergence")),
                     limits = c(min(rwexp_vec), max = max(rwexp_vec)+.1)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7), legend.position = "right", legend.box="vertical",
         plot.title = element_text(family = "Arial", size = 7),
        legend.key.size = unit(.75, 'lines')) +
  guides(
    colour = guide_colourbar(order = 1, barheight = 4, barwidth = 0.75),
    size = guide_legend(order = 2))

# PCM
p.pcm <- ggplot(freq_pcm, aes(x=md, y=mu)) +
  geom_point(alpha=0.7, aes(colour = rel_err, size = rel_freq)) +
  scale_size(range = c(0.5, 3), name="Relative\nfrequency (%)",
             breaks = c(5, 15, 60)) +
  scale_colour_gradient(high = "#1B4F72", low = "#85C1E9",
                        name = "Median\nrelative error") +
  scale_x_continuous(expression(paste("Mass-to-drag ratio (m)")),
                     limits = c(min(pcm_md_vec), max = max(pcm_md_vec))) +
  scale_y_continuous(expression(paste("Sliding friction coefficient")),
                     limits = c(min(pcm_mu_vec), max = max(pcm_mu_vec))) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7), legend.position = "right", legend.box="vertical",
        plot.title = element_text(family = "Arial", size = 7),
        legend.key.size = unit(.75, 'lines')) +
  guides(
    colour = guide_colourbar(order = 1, barheight = 4, barwidth = 0.75),
    size = guide_legend(order = 2))


setwd("/home/jason/Scratch/Figures")

#p.rwpcm <- p.rw + p.pcm + plot_layout(ncol = 1) + plot_annotation(tag_levels = '(a)')
p.rw + p.pcm

ggsave("rwpcm_spcv_scatter.png", dpi = 300, width = 7.5, height = 2.5, units = "in")


# GET RW AND PCM PERFORMANCE FOR INDV EVENTS ###################################

# PCM

pcm_opt <- pcmGetOpt(pcm_gridsearch_multi)

pcm_md_vec <- as.numeric(colnames(pcm_gridsearch_multi[[1]][[1]]))
pcm_mu_vec <- as.numeric(rownames(pcm_gridsearch_multi[[1]][[1]]))

# Index performance for opt param set
ind_md <- which(pcm_md_vec == pcm_opt$pcm_md)
ind_mu <- which(pcm_mu_vec == pcm_opt$pcm_mu)

relerr_single <- rep(NA, length(pcm_gridsearch_multi))
for(i in 1:length(pcm_gridsearch_multi)){
  relerr_single[i] <- pcm_gridsearch_multi[[i]]$relerr[ind_mu, ind_md]
}

# RW
rw_opt <- rwGetOpt(rw_gridsearch_multi)

rwslp_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[1]])
rwexp_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[2]])
rwper_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[3]])

ind_slp <- which(rwslp_vec == rw_opt$rw_slp_opt)
ind_exp <- which(rwexp_vec == rw_opt$rw_exp_opt)
ind_per <- which(rwper_vec == rw_opt$rw_per_opt)

#roc[row, col, matrix]
#roc[rwslp, rwexp, rwper]

auroc_single <- rep(NA, length(rw_gridsearch_multi))
for(i in 1:length(rw_gridsearch_multi)){
  auroc_single[i] <- rw_gridsearch_multi[[i]][ind_slp, ind_exp, ind_per]
}



# MAP PERFORMANCES #############################################################

opt_auroc <- auroc_single
opt_relerr <- relerr_single


dflow_pnt <- gCentroid(slide_poly_vec, byid=TRUE, id = slide_poly_vec$objectid)

dflow_data <- data.frame(objectid = 1:100, opt_auroc = opt_auroc,
                         opt_relerr = opt_relerr)
dflow_pnt <- SpatialPointsDataFrame(dflow_pnt, dflow_data)

#https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
dflow_sf <- st_as_sf(dflow_pnt)
dflow_sf$x <- coordinates(dflow_pnt)[,1]
dflow_sf$y <- coordinates(dflow_pnt)[,2]

dflow_df <- data.frame(dflow_sf)
dflow_df$opt_auroc <- auroc_single
dflow_df$opt_relerr <- relerr_single


m.relerr <-  ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +
  geom_sf(data = boundary, alpha = 0.8, fill = "white") +
  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
  geom_point(data = dflow_sf, size = 2, alpha=0.8, aes(x = x, y = y, colour = dflow_df$opt_relerr)) +
  scale_colour_gradient(low = "#D0ECE7", high = "#0B5345",
                        name = "Runout distance\nrelative error") +
  geom_point(data = dflow_sf, shape = 1, size = 2, alpha = 0.9, colour="#7f8c8d", aes(x = x, y = y)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_text(size = 7), legend.position = "right")

m.auroc <-  ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +
  geom_sf(data = boundary, alpha = 0.8, fill = "white") +
  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
  geom_point(data = dflow_sf, size = 2, alpha=0.8, aes(x = x, y = y, colour = dflow_df$opt_auroc)) +
  scale_colour_gradient(high = "#e6b0aa", low = "#7b241c",
                        name = "Runout path\nAUROC") +
  geom_point(data = dflow_sf, shape = 1, size = 2, alpha = 0.9, colour="#7f8c8d", aes(x = x, y = y)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_text(size = 7), legend.position = "right")

m.perf <- m.auroc + m.relerr
m.perf

setwd("/home/jason/Scratch/Figures")
ggsave("map_relerr_auroc_bottom_fnt7.png", dpi = 300, width = 7.5, height = 6.7, units = "in")


# CALCULATE DISTANCE ERROR #####################################################

setwd("/home/jason/Scratch/GPP_PCM_Paper")

geom_gpp <- data.frame(slide_id = 1:100, est_lng = NA, obs_lng = NA, est_wd = NA, obs_wd = NA)

for(i in 1:100){
  print(i)
  gpp <- pcmPerformance(dem = dem, slide_plys = slide_poly_vec, release_pnts = slide_point_vec,
                          slide_id = i, rw_slp = 40, rw_ex = 3, rw_per = 1.9, pcm_mu = 0.11, pcm_md = 40,
                          buffer_ext = 5000, buffer_source = 50, gpp_iter = 1000, rescale_threshold = 0.5, plot_eval = TRUE, return_features = TRUE)

  geom_obs <- runoutGeom(gpp$actual.poly, elev = gpp$dem)

  # length of estimated
  pred_thres <- rasterCdf(gpp$gpp.parea)
  pred_thres[pred_thres >= 0.5] = 1
  pred_thres[pred_thres < 0.5] = NA
  sp.pred <- rasterToPolygons(pred_thres, n = 4, dissolve = TRUE, na.rm=TRUE)
  geom_est <- runoutGeom(sp.pred, elev = gpp$dem)

  geom_gpp$est_lng[i] <- geom_est$length
  geom_gpp$obs_lng[i] <- geom_obs$length
  geom_gpp$est_wd[i] <- geom_est$width
  geom_gpp$obs_wd[i] <- geom_obs$width

}

geom_gpp$lng_err <- geom_gpp$est_lng - geom_gpp$obs_lng
geom_gpp$wd_err <- geom_gpp$est_wd - geom_gpp$obs_wd


#save(geom_gpp, file = "obs_est_opt_geom.Rd")
(load("obs_est_opt_geom.Rd"))


# EXPLORE ERROR PATTERNS #######################################################

# Correlation between predicted and observed runout distances
cor(geom_gpp$obs_lng, geom_gpp$est_lng, method = "spearman")

# Plot patterns
setwd("/home/jason/Scratch/Figures")

png(filename="hist_error_plots.png", res = 300, width = 7.5, height = 6,
    units = "in", pointsize = 11)

par(family = "Arial", mfrow = c(2,2), mar = c(4, 3, 0.5, 0.5),
    mgp = c(2, 0.75, 0))

hist(dflow_df$opt_auroc, xlab = "Runout path AUROC", main = "", breaks = 20)

hist(dflow_df$opt_relerr, xlab = "Runout distance relative error", main = "", breaks = 20)

hist(geom_gpp$lng_err, xlab = "Runout distance error (m)", breaks = 40,
     cex.axis = 1, cex.lab = 1, col = "lightgrey", main = "" )

plot(geom_gpp$est_lng, geom_gpp$obs_lng,
    xlab = "Estimated runout distance (m)",
    ylab = "Observed runout distance (m)", pch = 20)

dev.off()

# ERROR AND TERRAIN ATTRIBUTES #################################################

geom_gpp$rel_err <- (abs(geom_gpp$est_lng - geom_gpp$obs_lng)) / geom_gpp$obs_lng

slide_poly_vec$slide_id <- 1:100

sel_over_pnt  <- sp::over(slide_point_vec, slide_poly_vec)
smp_pnt  <- slide_point_vec[!is.na(sel_over_pnt$slide_id),]

library(sf)
smp_pnt <- st_as_sf(smp_pnt)
smp_ply <- st_as_sf(slide_poly_vec)

# Join point slides to poly slides
smp_pnt <- st_join(smp_pnt, smp_ply)
smp_pnt [ order(smp_pnt$slide_id , decreasing = FALSE ),]


smp_pnt$lng_err <- geom_gpp$lng_err[match(smp_pnt$slide_id, geom_gpp$slide_id)]
smp_pnt$rel_err <- geom_gpp$rel_err[match(smp_pnt$slide_id, geom_gpp$slide_id)]

setwd("/home/jason/Data/Chile/saga_data")
dem <- raster("dem_alos_12_5m_fill.sdat")
slope <- raster("slope_alos_12_5m_filled.sdat")
carea <- raster("carea_alos_12_5m_filled.sdat")

layers <- stack(dem, slope, carea)
names(layers) <- c('elev', 'slope', 'carea')
extr <- extract(layers,  as_Spatial(smp_pnt), df = TRUE)

smp_pnt$elev <- extr$elev
smp_pnt$slope <- extr$slope
smp_pnt$logcarea <- log10(extr$carea)
smp_pnt$carea <- extr$carea

cor(smp_pnt$elev, smp_pnt$rel_err, method = "spearman")
cor(smp_pnt$slope, smp_pnt$rel_err, method = "spearman")
cor(smp_pnt$logcarea, smp_pnt$rel_err, method = "spearman")

cor(geom_gpp$obs_lng, geom_gpp$rel_err, method = "spearman")
cor(geom_gpp$est_lng, geom_gpp$rel_err, method = "spearman")


# INDIVIDUAL MODEL PERFORMANCES ################################################

# RW

indv_rw <- rwGetOpt_single(rw_gridsearch_multi[[1]])

for(i in 2:length(rw_gridsearch_multi)){
  print(i)
  print(rwGetOpt_single(rw_gridsearch_multi[[i]]))
  indv_rw <- rbind(indv_rw, rwGetOpt_single(rw_gridsearch_multi[[i]]))
}

indv_rw

opt_sets <- data.frame(slp = indv_rw$rw_slp_opt, per = indv_rw$rw_per_opt, exp = indv_rw$rw_exp_opt)
freq_sets <- table(opt_sets)

slp_nm <- as.numeric(attributes(freq_sets)$dimnames$slp)
per_nm <- as.numeric(attributes(freq_sets)$dimnames$per)
exp_nm <- as.numeric(attributes(freq_sets)$dimnames$exp)

set_ind <- which(freq_sets !=0, arr.ind = TRUE)
sets <- data.frame(slp = slp_nm[set_ind[,1]],
                   per = per_nm[set_ind[,2]],
                   exp = exp_nm[set_ind[,3]],
                   freq = freq_sets[set_ind])


sets$rel_freq <- sets$freq/nrow(indv_rw) * 100

#Find median AUROC values
set_id <- paste(sets$slp, sets$per, sets$exp)
indv_rw$set_id <- paste(indv_rw$rw_slp_opt, indv_rw$rw_per_opt, indv_rw$rw_exp_opt)

for(i in 1:length(set_id)){
  auroc_i <- indv_rw$rw_auroc[which(indv_rw$set_id == set_id[i])]
  sets$median_auroc[i] <- median(auroc_i, na.rm = TRUE)
  sets$mean_auroc[i] <- mean(auroc_i, na.rm = TRUE)
}


# Get values for PCM

# PCM

indv_pcm <- pcmGetOpt_single(pcm_gridsearch_multi[[1]])

for(i in 2:length(pcm_gridsearch_multi)){
  print(i)
  indv_pcm <- rbind(indv_pcm, pcmGetOpt_single(pcm_gridsearch_multi[[i]]))
}

df_pcm <- data.frame(mu = indv_pcm$pcm_mu, md = indv_pcm$pcm_md)

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
pairs$rel_freq <- pairs$freq/nrow(indv_pcm) * 100


# find relative error
pair_id <- paste(pairs$mu, pairs$md)
indv_pcm$pair_id <- paste(indv_pcm$pcm_mu, indv_pcm$pcm_md)

for(i in 1:length(pair_id)){

  relerr_i <- indv_pcm$relerr[which(indv_pcm$pair_id == pair_id[i])]
  pairs$rel_err[i] <- median(relerr_i)
  pairs$iqr_relerr[i] <- IQR(relerr_i, na.rm = TRUE)

}


# Plot results

mybreaks <- c(1, 2, 3)

gg.rw.slope <- ggplot(sets, aes(x=per, y=exp)) +
  # Can improve by making different colors for different slope thresholds...

  geom_point(alpha=0.7, aes(col = slp, size = rel_freq)) +

  scale_size(breaks = mybreaks, range = c(0.5, 3)) +

  scale_colour_viridis_c(name = "Slope\nthreshold (°)", option = "cividis", direction = -1) +

  labs(size = "Relative\nfrequency (%)", col = "Slope\nthreshold") +

  #ggtitle("(a) Random walk") +

  guides(colour = guide_coloursteps()) +

  scale_x_continuous(expression(paste("Persistence factor")),
                              limits = c(min(rwper_vec), max = max(rwper_vec))) +
  scale_y_continuous(expression(paste("Exponent of divergence")),
                               limits = c(min(rwexp_vec), max = max(rwexp_vec)+.1)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7), legend.position = "bottom", legend.box="vertical",
        legend.margin=margin(), plot.title = element_text(family = "Arial", size = 7)) +
  guides(
    colour = guide_colourbar(order = 1, barwidth = 5, barheight = 0.6 ,label.position = "bottom"),
    size = guide_legend(order = 2, label.position = "bottom"))

gg.rw.slope


# Now of performance
mybreaks <- c(1, 2, 3)

gg.rw.auroc <- ggplot(sets, aes(x=per, y=exp)) +
  # Can improve by making different colors for different slope thresholds...

  geom_point(alpha=0.7, aes(col = median_auroc, size = rel_freq)) +

  scale_size(breaks = mybreaks, range = c(0.5, 3)) +

  scale_colour_gradient(low = "#1B4F72", high = "#85C1E9", name = "Median\nAUROC") +

  labs(size = "Relative\nfrequency (%)", col = "Slope\nthreshold") +

  #ggtitle("(b) Random walk") +

  scale_x_continuous(expression(paste("Persistence factor")),
                     limits = c(min(rwper_vec), max = max(rwper_vec))) +
  scale_y_continuous(expression(paste("Exponent of divergence")),
                     limits = c(min(rwexp_vec), max = max(rwexp_vec)+.1)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7), legend.position = "bottom", legend.box="vertical",
        legend.margin=margin(), plot.title = element_text(family = "Arial", size = 7)) +
  guides(
    colour = guide_colourbar(order = 1, barwidth = 5, barheight = 0.6 ,label.position = "bottom"),
    size = guide_legend(order = 2, label.position = "bottom"))

gg.rw.auroc


#or
#library(plotly)
#plot_ly(x=sets$per, y=sets$exp, z=sets$slp, type="scatter3d", mode="markers",
#        color=sets$mean_auroc, size =sets$rel_freq)


mybreaks <- c(1, 2, 3)

gg.ind.pcm <- ggplot(pairs, aes(x=md, y=mu)) +
  geom_point(alpha=0.7, aes(colour = rel_err, size = rel_freq)) +


  scale_size(name="Relative\nfrequency (%)", range = c(0.5, 3), breaks = mybreaks) +

  scale_colour_gradient(high = "#1B4F72", low = "#85C1E9", name = "Median\nrelative error") +

  #ggtitle("(c) PCM") +
  xlab("Mass-to-drag ratio (m)") +
  ylab("Sliding friction coefficient") +


  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7), legend.position = "bottom", legend.box="vertical",
        legend.margin=margin(), plot.title = element_text(family = "Arial", size = 7)) +
  guides(
    colour = guide_colourbar(order = 1, barwidth = 5, barheight = 0.6 ,label.position = "bottom"),
    size = guide_legend(order = 2, label.position = "bottom"))
gg.ind.pcm

#setwd("/home/jason/Scratch/Figures")
#ggsave("indv_scatter_opt_pcm.png", dpi = 300, width = 5.5, height = 3.25, units = "in")


gg.rw.slope + gg.rw.auroc + gg.ind.pcm
ggsave("indv_scatter_rw_pcm_perf.png", dpi = 300, width = 7.5, height = 3.25, units = "in")

# MAP OF INDV PERFORMANCE ######################################################

dflow_df$md <- indv_pcm$pcm_md
dflow_df$mu <- indv_pcm$pcm_mu

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
  geom_point(data = dflow_df, alpha=0.8, aes(x = x, y = y, size = md, colour = mu)) +
  scale_colour_gradient(low = "#fef9e7", high = "#117864",
                        name = "Sliding friction\ncoefficient") +
  geom_point(data = dflow_df, shape = 1, fill = NA, alpha = 0.9, colour="#7f8c8d", aes(x = x, y = y, size = md)) +
  scale_size(name="Mass-to-drage\nratio (m)", range = c(1, 3), breaks = c(20, 75, 150 )) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),legend.position = "right")
m.pcm_sol

setwd("/home/jason/Scratch/Figures")
ggsave("map_inv_opt_param.png", dpi = 300, width = 3.5, height = 6.7, units = "in")


# MAP OF RUNOUT ################################################################

setwd("/home/jason/Scratch/GPP_Subcatch_Paper")
parea_cdf <- raster("cdf_gpp_parea_all_cutoff_0.7.tif")
parea_cdf <- aggregate(parea_cdf, fact = 6, fun = mean)

# Prep for ggplot
parea_df <- as.data.frame(parea_cdf, xy = TRUE)
parea_df <- parea_df[!is.na(parea_df[,3]),]


# For close up
setwd("/home/jason/Data/Chile/")


# Load shapefile of sub catchment polygons
sub_catchments <- readOGR("sub_catchments.shp")
sub_catchment <- sub_catchments[sub_catchments$ID == 56,]

# Crop hillshade
hs_sub <- crop(raster("saga_data/hs_filtered_filled.sdat"), sub_catchment)
hs_sub <- mask(hs_sub, sub_catchment)

# Load and crop source area prediction raster
parea_sub <- crop(raster("/home/jason/Scratch/GPP_Subcatch_Paper/cdf_gpp_parea_all_cutoff_0.7.tif"),
                  sub_catchment)
parea_sub <- mask(parea_sub, sub_catchment)

# Aggregate for quicker mapping
parea_sub <- aggregate(parea_sub, fact = 2, fun = mean)
hs_sub <- aggregate(hs_sub, fact = 2, fun = mean)

# Prep for ggplot
parea_sub_df <- as.data.frame(parea_sub, xy = TRUE)
parea_sub_df <- parea_sub_df[!is.na(parea_sub_df[,3]),]

hs_sub_df <- as.data.frame(hs_sub, xy = TRUE)
hs_sub_df <- hs_sub_df[!is.na(hs_sub_df[,3]),]

# Clip river feature
bnd_sub <- st_as_sf(sub_catchment)
river_sub <- st_intersection(bnd_sub_sf, rivers)


# Load and crop source area prediction raster
source_pred <- raster("source_pred_gam.tif")
source_sub <- crop(source_pred, sub_catchment)
source_sub <- mask(source_sub, sub_catchment)

source_pred <- aggregate(source_pred, fact = 6, fun = mean)
source_sub <- aggregate(source_sub, fact = 2, fun = mean)

# Prep for ggplot
source_df <- as.data.frame(source_pred, xy = TRUE)
source_df <- source_df[!is.na(source_df[,3]),]

source_sub_df <- as.data.frame(source_sub, xy = TRUE)
source_sub_df <- source_sub_df[!is.na(source_sub_df[,3]),]

# Crop



# Source area prediction map ###################################################

#

library(ggnewscale)
library(ggspatial)

map.source <- ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +

  geom_sf(data = boundary, alpha = 0.4, fill = "white") +

  new_scale("fill") +
  geom_raster(data=source_df, aes(x=x, y=y, fill = source_pred_gam) ) +
  scale_fill_viridis_c(name = "Source area\nsusceptibility\n(probability)",  alpha = 0.6, direction = -1) +

  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  #geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +

  geom_sf(data = bnd_sub, colour = "#EC7063", fill = NA) +

  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_text(size = 7), legend.position = "right")

#map.source


# For sub catchment


map.source_sub <- ggplot() +
  geom_sf() +
  geom_raster(data=hs_sub_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +

  geom_sf(data = bnd_sub, alpha = 0.4, fill = "white") +

  new_scale("fill") +
  geom_raster(data=source_sub_df, aes(x=x, y=y, fill = source_pred_gam),
              show.legend = FALSE) +
  scale_fill_viridis_c(name = "Source area\nsusceptibility\n(probability)",  alpha = 0.6, direction = -1) +

  geom_sf(data = river_sub, colour = "#85C1E9") +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  #geom_sf(data = bnd_sub, colour = "#7f8c8d", fill = NA) +

  annotation_scale(aes(style = "ticks", location = "br"), text_family = "Arial", text_cex = 0.6,
                   bar_cols = "black", line_width = 0.7) +

  xlab("") +
  ylab("") +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right")

#map.source_sub

map.source + map.source_sub

setwd("/home/jason/Scratch/Figures")
ggsave("regional_source_pred_map.png", dpi = 300, width = 7.5, height = 6, units = "in")




# Process area map #############################################################

#

library(ggnewscale)
library(ggspatial)

map.parea <- ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +

  geom_sf(data = boundary, alpha = 0.4, fill = "white") +

  new_scale("fill") +
  geom_raster(data=parea_df, aes(x=x, y=y, fill = cdf_gpp_parea_all_cutoff_0.7) ) +
  scale_fill_viridis_c(name = "Runout frequency\n(quantiles)",  alpha = 0.6, direction = -1) +

  geom_sf(data = rivers, colour = "#85C1E9") +
  geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  #geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +

  geom_sf(data = bnd_sub, colour = "#EC7063", fill = NA) +

  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_text(size = 7), legend.position = "right")

#map.parea


# For sub catchment

map.parea_sub <- ggplot() +
  geom_sf() +
  geom_raster(data=hs_sub_df, aes(x=x, y=y, fill = hs_filtered_filled),
              show.legend = FALSE) +
  scale_fill_gradient(high = "black", low = "white", na.value = "#FFFFFF") +

  geom_sf(data = bnd_sub, alpha = 0.4, fill = "white") +

  new_scale("fill") +
  geom_raster(data=parea_sub_df, aes(x=x, y=y, fill = cdf_gpp_parea_all_cutoff_0.7),
              show.legend = FALSE) +
  scale_fill_viridis_c(name = "Runout frequency\n(quantiles)",  alpha = 0.6, direction = -1) +

  geom_sf(data = river_sub, colour = "#85C1E9") +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = bnd_sub, colour = "#7f8c8d", fill = NA) +

  xlab("") +
  ylab("") +
  theme_void() +
  annotation_scale(aes(style = "ticks", location = "br"), text_family = "Arial", text_cex = 0.6,
                  bar_cols = "black", line_width = 0.7) +


  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right")


#map.parea_sub

map.parea + map.parea_sub

setwd("/home/jason/Scratch/Figures")
ggsave("regional_runout_impact_map.png", dpi = 300, width = 7.5, height = 6, units = "in")

