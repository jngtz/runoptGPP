
# Libraries ---------------------------------------------------
library(runoptGPP)
library(rgdal)
library(raster)
library(rgeos)
library(sf)

library(Rsagacmd)
saga <- saga_gis()


#myenv <- rsaga.env(path="D:\\Software\\saga-6.1.0_x64", parallel = TRUE)

# FUNCTIONS ######################################################################


# DATA ############################################################################

# Set workspace
setwd("/home/jason/Data/Chile/")

# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m.tif")

# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons and assign object ID based on row number
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
runout_polygons$objectid <- 1:length(runout_polygons)

# Select a debris flow and source point for this example
runout_polygon <- runout_polygons[77,]

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]

setwd("saga_data")
hillshade <- raster("hs_filtered_filled.sdat")

setwd("/home/jason/Scratch/Figures")


# PCM Optimization #################################################################################

pretty_PCM_example <- function(dem, slide_ply, slide_src, src_thresh=0.5, MD, MU, hillshade){

  require(broom)
  require(tibble)
  require(dplyr)


  gpp <- pcmPerformance(dem, slide_plys = slide_ply, slide_src = slide_src,
                        rw_slp = 40, rw_ex = 3, rw_per = 1.5,
                        pcm_mu = MU, pcm_md = MD,
                        gpp_iter = 1000, buffer_ext = 600, buffer_source = 50,
                        plot_eval = TRUE, return_features = TRUE, saga_lib = saga)

  # Get bounding boxes for obs. and simulated runout
  pred_thres <- rasterCdf(gpp$gpp.parea)
  pred_thres[pred_thres >= src_thresh] = 1
  pred_thres[pred_thres < src_thresh] = NA
  pred_poly <- raster::rasterToPolygons(pred_thres, n = 4, dissolve = TRUE, na.rm=TRUE)

  bb_est <- minBBoxSpatialPolygons(pred_poly)
  bb_obs <- minBBoxSpatialPolygons(slide_ply)

  obs_df <- tidy(gpp$actual.poly, region = "objectid" )
  bbEst_df <- tidy(bb_est, region = "value" )
  bbObs_df <- tidy(bb_obs, region = "value" )

  obs_df <- obs_df %>% add_column(label = "Obs. debris flow runout")
  bbObs_df <- bbObs_df %>% add_column(label = "Obs. min. area bounding box")
  bbEst_df <- bbEst_df %>% add_column(label = "Est. min. area bounding box")

  obs_df <- obs_df %>% add_column(fill_ply = 1)
  bbObs_df <- bbObs_df %>% add_column(fill_ply = 2)
  bbEst_df <- bbEst_df %>% add_column(fill_ply = 3)

  layers <- bind_rows(list(obs_df, bbObs_df, bbEst_df))

  parea_df <- as.data.frame(rasterCdf(gpp$gpp.parea), xy = TRUE)
  names(parea_df)[3] <- "value"
  parea_df <- parea_df[!is.na(parea_df$value),]

  hs <- crop(hillshade, gpp$dem)

  hs_df <-  as.data.frame(hs, xy = TRUE)

  return(list(
    parea_df = parea_df,
    hs_df = hs_df,
    layers = layers,
    mu = MU,
    md = MD,
    roc = gpp$roc,
    length.relerr = gpp$length.relerr
  ))
}


#https://egallic.fr/en/scale-bar-and-north-arrow-on-a-ggplot2-map/

# Run GPP for Opt Mu example


gpp <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                      rw_slp = 40, rw_ex = 3, rw_per = 1.5,
                      pcm_mu = 0.11, pcm_md = 40,
                      gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                      plot_eval = TRUE, return_features = TRUE, saga_lib = saga)

#save(gpp, file = "gpp_example_opt_plot.Rd")

gppA <- pretty_PCM_example(dem, runout_polygon, source_point,
                           src_thresh=0.5, MD = 40, MU = 0.03, hillshade)
#save(gppA, file = "gppA_example_opt_plot.Rd")

gppB <- pretty_PCM_example(dem, runout_polygon, source_point,
                           src_thresh=0.5, MD = 40, MU = 0.1, hillshade)
#save(gppB, file = "gppB_example_opt_plot.Rd")

gppC <- pretty_PCM_example(dem, runout_polygon, source_point,
                           src_thresh=0.5, MD = 40, MU = 0.5, hillshade)
#save(gppC, file = "gppC_example_opt_plot.Rd")


# GGPLOTS  ##################################################################
library(ggplot2)
library(patchwork)

#(load("gpp_example_opt_plot.Rd"))
#(load("gppA_example_opt_plot.Rd"))
#(load("gppB_example_opt_plot.Rd"))
#(load("gppC_example_opt_plot.Rd"))


# define plot limits

minx <- min(c(
  minxA <- min( c(min(gppA$layers$long), min(gppA$parea_df$x)) ),
  minxB <- min( c(min(gppB$layers$long), min(gppB$parea_df$x)) ),
  minxC <- min( c(min(gppC$layers$long), min(gppC$parea_df$x)) )
))

maxx <- min(c(
  maxxA <- max( c(max(gppA$layers$long), max(gppA$parea_df$x)) ),
  maxxB <- max( c(max(gppB$layers$long), max(gppB$parea_df$x)) ),
  maxxC <- max( c(max(gppC$layers$long), max(gppC$parea_df$x)) )
))

miny <- min(c(
  minyA <- min( c(min(gppA$layers$lat), min(gppA$parea_df$y)) ),
  minyB <- min( c(min(gppB$layers$lat), min(gppB$parea_df$y)) ),
  minyC <- min( c(min(gppC$layers$lat), min(gppC$parea_df$y)) )
))

maxy <- max(c(
  maxyA <- max( c(max(gppA$layers$lat), max(gppA$parea_df$y)) ),
  maxyB <- max( c(max(gppB$layers$lat), max(gppB$parea_df$y)) ),
  maxyC <- max( c(max(gppC$layers$lat), max(gppC$parea_df$y)) )
))


maxx <- maxx + 100

src_pnt <- coordinates(gpp$source.pnt)
src_pnt <- as.data.frame(src_pnt)
names(src_pnt) <- c("x", "y")

# mu 0.05

# reproject extent
( e <- raster::extent(minx, maxx, miny, maxy) )
e <- as(e, "SpatialPolygons")
sp::proj4string(e) <- crs(gpp$dem)

e.geo <- sp::spTransform(e, CRS("+proj=longlat +datum=WGS84 +no_defs
                             +ellps=WGS84 +towgs84=0,0,0"))


geo.bbox <- c(
  left = extent(e.geo)[1],
  bottom = extent(e.geo)[3],
  right = extent(e.geo)[2],
  top = extent(e.geo)[4]
)

src_pnt_wgs84 <- gpp$source.pnt
src_pnt_wgs84 <- spTransform(src_pnt_wgs84, CRS("+proj=longlat +datum=WGS84"))
src_pnt_wgs84_crd <- coordinates(src_pnt_wgs84 )

lat_pnt <- as.numeric(src_pnt_wgs84_crd[,2])
lon_pnt <- as.numeric(src_pnt_wgs84_crd[,1])

# A

library(ggnewscale)
library(ggspatial)

m.muA <- ggplot() +

  geom_raster(data = gppA$parea_df, aes(x = x, y = y, fill = value), show.legend = FALSE) +
  geom_polygon(data = gppA$layers, alpha=0.2, fill = NA,  aes(x = long, y = lat, col = label, linetype = label), show.legend = FALSE) +
  geom_point(data = src_pnt, fill = "black", aes(x=x, y=y)) +


  #scale_fill_viridis_c(name = "Runout frequency\n(Quantile)", direction = -1) +


  scale_fill_gradient(low = "#F8F9F9", high = "#2874A6", name = "Runout frequency\n(Quantile)") +


  scale_color_manual(name = "", values = c("black", "black", "red")) +
  scale_linetype_manual(name = "", values = c("solid", "dotdash", "dashed")) +

  ggtitle(paste("Sliding friction coefficient:", gppA$mu,
                "\nRelative error:", round(gppA$length.relerr, digits = 2),
                " AUROC:", round(gppA$roc, digits = 2))) +
  xlab("Easting") +
  ylab("Northing") +

  xlim(minx, maxx) +
  ylim(miny, maxy) +

  scale_x_continuous(breaks=seq(397400, 398800, 200), lim = c(397500, 398800)) +

  coord_fixed()+
  theme_void() +
  #annotation_scale(aes(style = "ticks", location = "br"), text_family = "Arial", text_cex = 0.6,
  #                 bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right",
        plot.title = element_text(family = "Arial", size = 8))


# B

m.muB <- ggplot() +
  #geom_sf() +
  #geom_raster(data = hs_df, aes(x = x, y=y, alpha = hs_filtered_filled )) +
  geom_raster(data = gppB$parea_df, aes(x = x, y = y, fill = value), show.legend = FALSE) +
  geom_polygon(data = gppB$layers, alpha=0.2, fill = NA,  aes(x = long, y = lat, col = label, linetype = label), show.legend = FALSE) +
  geom_point(data = src_pnt, fill = "black", aes(x=x, y=y)) +


  #scale_fill_viridis_c(name = "Runout frequency\n(Quantile)", direction = -1) +


  scale_fill_gradient(low = "#F8F9F9", high = "#2874A6",
                      name = "Runout frequency\n(Quantile)") +

  scale_color_manual(name = "", values = c("black", "black", "red")) +
  scale_linetype_manual(name = "", values = c("solid", "dotdash", "dashed")) +

  ggtitle(paste("Sliding friction coefficient:", gppB$mu,
                "\nRelative error:", round(gppB$length.relerr, digits = 2),
                " AUROC:", round(gppB$roc, digits = 2))) +
  xlab("Easting") +
  ylab("Northing") +

  xlim(minx, maxx) +
  ylim(miny, maxy) +

  scale_x_continuous(breaks=seq(397400, 398800, 200), lim = c(397500, 398800)) +

  coord_fixed()+
  theme_void() +
  #annotation_scale(aes(style = "ticks", location = "br"), text_family = "Arial", text_cex = 0.6,
  #                 bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right",
        plot.title = element_text(family = "Arial", size = 8))



# C

m.muC <- ggplot() +
  #geom_sf() +
  geom_raster(data = gppC$parea_df, aes(x = x, y = y, fill = value)) +
  geom_polygon(data = gppC$layers, alpha=0.2, fill = NA,  aes(x = long, y = lat, col = label, linetype = label), show.legend = TRUE) +
  geom_point(data = src_pnt, fill = "black", aes(x=x, y=y)) +

  #scale_fill_viridis_c(name = "Runout frequency\n(Quantile)", direction = -1) +

  scale_fill_gradient(low = "#F8F9F9", high = "#2874A6",
                      name = "Runout frequency\n(quantile)") +

  scale_color_manual(name = "", values = c("black", "red", "black"),
                     labels = c("Est. min. area\nbounding box",
                                "Obs. min. area\nbounding box",
                                "Obs. debris flow\nrunout")) +

  scale_linetype_manual(name = "", values = c("solid", "dashed", "dotdash"),
                        labels = c("Est. min. area\nbounding box",
                                   "Obs. min. area\nbounding box",
                                   "Obs. debris flow\nrunout")) +

  ggtitle(paste("Sliding friction coefficient:", gppC$mu,
                "\nRelative error:", round(gppC$length.relerr, digits = 2),
                " AUROC:", round(gppC$roc, digits = 2))) +
  xlab("Easting") +
  ylab("Northing") +

  xlim(minx, maxx) +
  ylim(miny, maxy) +

  #guides(fill=guide_legend(
  #  keywidth=0.5,
  #  keyheight=0.5,
  #  default.unit="cm")
  #) +

  scale_x_continuous(breaks=seq(397400, 398800, 200), lim = c(397500, 398800)) +

  coord_fixed()+

  xlab("") +
  ylab("") +
  theme_void() +
  annotation_scale(aes(style = "ticks", location = "br"), text_family = "Arial", text_cex = 0.6,
                   bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=7),
        plot.title = element_text(family = "Arial", size = 8)) +

  guides(fill = guide_colourbar(order = 1, barheight = 0.75, barwidth = 5))



# Plot all together with patchwork


#with patch work, but cannot maintain fixed scales (coord_fixed, coord_sf)
m.muA + m.muB + m.muC + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = 'bottom')

#plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_sep = '',tag_suffix = ')')

ggsave("map_3level_pcm_mu_opt_blue.png", dpi = 300, width = 7.5, height = 2.5, units = "in")

