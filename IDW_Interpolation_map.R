library(dplyr)
library(sp)
library(gstat)
library(raster)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)


# load f2 matrix
f2_mean <- read.table("f2_matrix.txt")

# remove outgroup
mat2 <- mat2[rownames(mat2) != 'Mbuti',
                colnames(mat2) != 'Mbuti']

# read metadata
meta <- read.table("metameta_for_interpolation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# choose time interval (in this case 14-11 ka BP)
df2 <- meta %>% filter(Date < 14000 & Date > 11000)

common_samples <- intersect(df2$Sample_name, rownames(mat2))
mat2 <- as.matrix(mat2)
mat2 <- mat2[common_samples, common_samples]

# compute MDS
mds <- cmdscale(mat2, k = 2)
mds_df <- as.data.frame(mds)
colnames(mds_df) <- c("PC1", "PC2")
mds_df$Sample_name <- rownames(mds_df)

df3 <- merge(df2, mds_df, by = "Sample_name")

coordinates(df3) <- ~Long + Lat
proj4string(df3) <- CRS("+proj=longlat +datum=WGS84")

# grid for interpolation
xlim <- c(-13.826067, 60.308187)
ylim <- c(35.782773, 64.910148)

grd <- expand.grid(
  Long = seq(xlim[1], xlim[2], by = 0.1),
  Lat  = seq(ylim[1], ylim[2], by = 0.1)
)

coordinates(grd) <- ~Long + Lat
gridded(grd) <- TRUE
proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")

# IDW interpolation of the first two MDS components
idw_PC1 <- idw(PC1 ~ 1, locations = df3, newdata = grd, idp = 2)
idw_PC2 <- idw(PC2 ~ 1, locations = df3, newdata = grd, idp = 2)

# mask sea using land polygons
world <- ne_countries(scale = "medium", returnclass = "sf")
world_sp <- as(world, "Spatial")

r_PC1_masked <- mask(raster(idw_PC1), world_sp)
r_PC2_masked <- mask(raster(idw_PC2), world_sp)

# convert to dataframes for ggplot
df_PC1_masked <- as.data.frame(r_PC1_masked, xy = TRUE)
colnames(df_PC1_masked) <- c("Long", "Lat", "PC1")

df_PC2_masked <- as.data.frame(r_PC2_masked, xy = TRUE)
colnames(df_PC2_masked) <- c("Long", "Lat", "PC2")

# plot dimension 1
p1 <- ggplot() +
  geom_raster(data = df_PC1_masked, aes(x = Long, y = Lat, fill = PC1), alpha = 0.8) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90") +
  new_scale_fill() +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.5) +
  geom_point(
    data = as.data.frame(df3),
    aes(x = Long, y = Lat, shape = group, fill = group),
    color = "black",
    stroke = 0.1,
    size = 3
  ) +
  scale_shape_manual(values = c(
    "East_11-8_kaBP" = 21,
    "East_14-11_kaBP" = 21,
    "East_17-14_kaBP" = 21,
    "South_11-8_kaBP" = 25,
    "South_14-11_kaBP" = 25,
    "South_17-14_kaBP" = 25,
    "West-North_11-8_kaBP" = 24,
    "West-North_14-11_kaBP" = 24,
    "West-North_17-14_kaBP" = 24
  )) +
  scale_fill_manual(values = c(
    "East_11-8_kaBP"  = "olivedrab1",
    "East_14-11_kaBP" = "olivedrab3",
    "East_17-14_kaBP" = "olivedrab4",
    "South_11-8_kaBP"  = "darkgoldenrod1",
    "South_14-11_kaBP" = "darkgoldenrod3",
    "South_17-14_kaBP" = "darkgoldenrod4",
    "West-North_11-8_kaBP"  = "#A6CEE3",
    "West-North_14-11_kaBP" = "#1F78B4",
    "West-North_17-14_kaBP" = "#08306B"
  )) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "IDW Interpolation – Dim1 (14-11 kaBP)",
       fill = "Group / Dim1", shape = "Group") +
  theme_minimal() +
  theme(panel.grid = element_blank())

p1
