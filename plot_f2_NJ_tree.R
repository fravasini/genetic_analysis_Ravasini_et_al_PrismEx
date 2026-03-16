library(ggtree)
library(ape)
library(phangorn)
library(ggnewscale)
library(dplyr)
library(ggplot2)

# load tree
tree <- read.tree("f2_tree.nwk")

# root tree
tree <- root(tree, outgroup = "Mbuti", resolve.root = TRUE)

# load metadata
meta <- read.table("metadata_for_tree_figure.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta$Sample_name <- as.character(meta$Sample_name)
meta$Sample_name1 <- meta$Sample_name

# find "clade marker" nodes
cultures <- sort(unique(meta$culture))
node_ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
threshold <- 1

node_info <- lapply(node_ids, function(nd) {
  desc_idx <- Descendants(tree, nd, type = "tips")[[1]]
  tips <- tree$tip.label[desc_idx]
  cats <- meta$culture[match(tips, meta$Sample_name)]
  n_total <- length(cats)                          # denominatore = TUTTI i tip
  cats <- cats[!is.na(cats)]                       # togli NA solo per trovare la dominante
  if (length(cats) == 0) return(data.frame(node = nd, dominant = NA_character_, prop = 0))
  tab <- table(factor(cats, levels = cultures))
  max_cat <- names(which.max(tab))
  max_prop <- max(tab) / n_total                   # proporzione su tutti i tip, NA inclusi
  data.frame(node = nd, dominant = max_cat, prop = max_prop, stringsAsFactors = FALSE)
})
node_info_df <- do.call(rbind, node_info)

# keep only highest node for clade
eligible <- node_info_df %>% filter(prop >= threshold, !is.na(dominant))

to_keep <- sapply(eligible$node, function(nd) {
  ancestors <- Ancestors(tree, nd, type = "all")
  !any(ancestors %in% setdiff(eligible$node, nd))
})

clade_markers <- eligible[to_keep, ]

# colors for cultures (ATU 2.1 in original publication)
culture_cols <- c(
  "Asturian"              = "#E41A1C",
  "Azilian"               = "#377EB8",
  "Beuronian"             = "#4DAF4A",
  "Butovo"                = "#984EA3",
  "Castelnovian"          = "#FF7F00",
  "Epigravettian"         = "#A65628",
  "Federmesser"           = "#F781BF",
  "Fosna"                 = "#E6AB02",
  "Iron_Gates_Mesolithic" = "#66C2A5",
  "Komornica"             = "#FC8D62",
  "Kongemose"             = "#8DA0CB",
  "Kukrek"                = "#E78AC3",
  "Kunda"                 = "#A6D854",
  "Magdalenian"           = "#FFD92F",
  "Maglemose"             = "#1B9E77",
  "Muge"                  = "#D95F02",
  "na"                    = "#999999",
  "Onega"                 = "#7570B3",
  "Rhine-Meuse-Scheld"    = "#E7298A",
  "Sauveterrian"          = "#66A61E",
  "Tardenoisian"          = "#CAB2D6",
  "Veretye"               = "#BEBADA"
)




# color and shape for groups (called sectors/time slices in original publication)
group_cols <- c(
  "East_11-8_kaBP"        = "olivedrab1",
  "East_14-11_kaBP"       = "olivedrab3",
  "East_17-14_kaBP"       = "olivedrab4",
  "South_11-8_kaBP"       = "darkgoldenrod1",
  "South_14-11_kaBP"      = "darkgoldenrod3",
  "South_17-14_kaBP"      = "darkgoldenrod4",
  "North-West_11-8_kaBP"  = "#A6CEE3",
  "North-West_14-11_kaBP" = "#1F78B4",
  "North-West_17-14_kaBP" = "#08306B"
)

group_shapes <- c(
  "East_11-8_kaBP"        = 21,
  "East_14-11_kaBP"       = 21,
  "East_17-14_kaBP"       = 21,
  "South_11-8_kaBP"       = 25,
  "South_14-11_kaBP"      = 25,
  "South_17-14_kaBP"      = 25,
  "North-West_11-8_kaBP"  = 24,
  "North-West_14-11_kaBP" = 24,
  "North-West_17-14_kaBP" = 24
)


# base tree plot
p <- ggtree(tree) %<+% meta +
  geom_tree() +
  theme_tree2() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 260, 5, 5))

tree_plot_data <- p$data

# build 'collapsed' clades
triangle_polys <- do.call(rbind, lapply(seq_len(nrow(clade_markers)), function(i) {
  nd <- clade_markers$node[i]
  cult <- clade_markers$dominant[i]
  
  mrca_row <- tree_plot_data[tree_plot_data$node == nd, ]
  
  desc_tips <- Descendants(tree, nd, type = "tips")[[1]]
  tip_rows <- tree_plot_data[tree_plot_data$node %in% desc_tips, ]
  
  data.frame(
    x = c(mrca_row$x, max(tip_rows$x), max(tip_rows$x)),
    y = c(mrca_row$y, min(tip_rows$y), max(tip_rows$y)),
    dominant = cult,
    id = i
  )
}))

# final plot
p_final <- p +
  # collapsed clade on full branches
  geom_polygon(data = triangle_polys,
               aes(x = x, y = y, group = id, fill = dominant),
               colour = "black", linewidth = 0.5, alpha = 1) +
  scale_fill_manual(values = culture_cols, name = "Culture") +
  
  new_scale_fill() +
  
  geom_tippoint(aes(shape = group, fill = group),
                size = 5, stroke = 0.8, color = "black") +
  scale_fill_manual(values = group_cols, name = "Group") +
  scale_shape_manual(values = group_shapes, name = "Group") +
  
  geom_tiplab(aes(label = Sample_name1),
              align = FALSE, offset = 0.005,
              linetype = "dotted", size = 5) +
  
  theme(legend.position = "right")

p_final


# IMPORTANT! the final plot was manually modified to remove the remaining branches under the collapsed triangles
