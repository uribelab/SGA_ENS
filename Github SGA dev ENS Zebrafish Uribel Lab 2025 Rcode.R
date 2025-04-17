

# Setting up the working directory (once downloaded expression and text file in any computer change the directory path).
setwd("C:/Users/NIBOR/Desktop/1-12-24 phox2bb harmony 4dfp3fish vs 7dpf1fish")

# Loading Giotto  and setting up the directory where spatial and expression files are saved (inside the working directory create a folder named gobject_data/ ).
library(Giotto)
rm(list = ls())
data_directory = paste0(getwd(),'/gobject_data/')

# Creating the expression and location objects for the larvae based on their location.
fourdpf_F1_exprs = paste0(data_directory, "1-13-25_4dpf_F1_R1-R4_expression.txt")
fourdpf_F1_locs = paste0(data_directory, "1-13-25_4dpf_F1_R1-R4_coordinates.txt")

fourdpf_F3_exprs = paste0(data_directory, "1-10-25_4dpf_F3_R1-R4_expression.txt")
fourdpf_F3_locs = paste0(data_directory, "1-10-25_4dpf_F3_R1-R4_coordinates.txt")

fourdpf_F7_exprs = paste0(data_directory, "1-10-25_4dpf_F7_R1-R4_expression.txt")
fourdpf_F7_locs = paste0(data_directory, "1-10-25_4dpf_F7_R1-R4_coordinates.txt")

fourdpf_F9_exprs = paste0(data_directory, "1-17-25_4dpf_F9_R1-R4_expression.txt")
fourdpf_F9_locs = paste0(data_directory, "1-17-25_4dpf_F9_R1-R4_coordinates.txt")

sevendpf_F6_exprs = paste0(data_directory, "1-10-25_7dpf_F6_R1-R4_expression.txt")
sevendpf_F6_locs = paste0(data_directory, "1-10-25_7dpf_F6_R1-R4_coordinates.txt")

sevendpf_F8_exprs = paste0(data_directory, "1-16-25_7dpf_F8_R1-R4_expression.txt")
sevendpf_F8_locs = paste0(data_directory, "1-16-25_7dpf_F8_R1-R4_coordinates.txt")

# Creating the Giotto Objects that combine expression and spatial location for the larvae.
fourdpf_F1_object = createGiottoObject(expression = fourdpf_F1_exprs,
                                       spatial_locs = fourdpf_F1_locs)

fourdpf_F3_object = createGiottoObject(expression = fourdpf_F3_exprs,
                                       spatial_locs = fourdpf_F3_locs)

fourdpf_F7_object = createGiottoObject(expression = fourdpf_F7_exprs,
                                       spatial_locs = fourdpf_F7_locs)

fourdpf_F9_object = createGiottoObject(expression = fourdpf_F9_exprs,
                                       spatial_locs = fourdpf_F9_locs)

sevendpf_F6_object = createGiottoObject(expression = sevendpf_F6_exprs,
                                        spatial_locs = sevendpf_F6_locs)

sevendpf_F8_object = createGiottoObject(expression = sevendpf_F8_exprs,
                                        spatial_locs = sevendpf_F8_locs)

# Join Giotto Objects by integrate data from different larvae.
combined4fish4dpf_vs_2fish7dpf_fishobject <- joinGiottoObjects(gobject_list = list(fourdpf_F1_object, fourdpf_F3_object, fourdpf_F7_object,fourdpf_F9_object, sevendpf_F6_object , sevendpf_F8_object),
                                                               gobject_names = c("fourdpf_F1", "fourdpf_F3", "fourdpf_F7", "fourdpf_F9", "sevendpf_F6_object","sevendpf_F8_object"),
                                                               join_method = "shift", x_padding = 100)

# Apply filtering of the integrated object.
combined4fish4dpf_vs_2fish7dpf_fishobject <- filterGiotto(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
                                                          expression_threshold = 10,
                                                          feat_det_in_min_cells = 10,
                                                          min_det_feats_per_cell = 10,
                                                          expression_values = c('raw'),
                                                          verbose = T)
# Normalize Giotto integrated object.
combined4fish4dpf_vs_2fish7dpf_fishobject <- normalizeGiotto(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject)

#run Principal Component Analysis (PC) of the integrated, normalized and filtered object. Also generate a scree plot.
combined4fish4dpf_vs_2fish7dpf_fishobject <- runPCA(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
                                                    expression_values = "scaled", 
                                                    scale_unit = FALSE, 
                                                    center = FALSE)
screePlot(combined4fish4dpf_vs_2fish7dpf_fishobject,
          ncp = 30
)

# Load and run Harmony R package algorithm designed for data integration and batch effect correction.
library(harmony)
combined4fish4dpf_vs_2fish7dpf_fishobject <- runGiottoHarmony(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                                                              vars_use = "list_ID",
                                                              dimensions_to_use = 1:11,
                                                              do_pca = FALSE,
                                                              sigma = 0.1,
                                                              theta = 2,
                                                              lambda = 1,
                                                              nclust = NULL)

# Add a nearest neighbor network.
combined4fish4dpf_vs_2fish7dpf_fishobject <- createNearestNetwork(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, dimensions_to_use = 1:11, k = 15)

#Use Leiden Cluster algorithm to partition the cells into distinct groups (clusters) based on the similarity network
combined4fish4dpf_vs_2fish7dpf_fishobject <- doLeidenClusterIgraph(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
                                                                   name = 'leiden_clus',
                                                                   spat_unit = 'cell',
                                                                   feat_type = 'rna',
                                                                   resolution = 1,
                                                                   #the lower the number for the resolution the less number of clusters. 1 has been used
                                                                   n_iterations = 1000)


# Run UMAP and plot UMAP based on leiden cluster and larvae coloring (Figure 1A and B).
combined4fish4dpf_vs_2fish7dpf_fishobject <- runUMAP(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                                                     dimensions_to_use = 1:11,
)

plotUMAP(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
         cell_color = "leiden_clus", 
         show_NN_network = FALSE,
         show_legend = FALSE,
         point_size = 5)

plotUMAP(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
         cell_color = "list_ID", 
         show_NN_network = FALSE,
         show_center_label = FALSE,
         show_legend = FALSE,
         point_size = 5)

# Generate spatial plots in 2D coloring the cells by their leiden clusters and larvae type (Figure 1D-I).
spatPlot2D(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, 
           cell_color = "leiden_clus",
           cow_n_col = 1,
           point_size = 3)
spatPlot2D(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, 
           group_by = "list_ID",
           cell_color = "list_ID",
           cow_n_col = 3,
           point_size = 6)

# Show Heat maps based on the leiden clustering and larvae type (Figure 2A and B).
showClusterHeatmap(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, 
                   expression_values = "normalized", 
                   cluster_column = "leiden_clus"
)
showClusterHeatmap(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, 
                   expression_values = "normalized", 
                   cluster_column = "list_ID"
)

# Plot UMAP coloring based on early (green) or late (violet) larvae (Figure 2C).
plotUMAP(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
         cell_color = "list_ID", 
         color_as_factor = TRUE,
         show_NN_network = FALSE,
         show_center_label = FALSE,
         cell_color_code = c ("darkseagreen","darkseagreen2", "darkseagreen3", "darkseagreen4", "violetred","violetred2"),
         show_legend = TRUE,
         point_size = 7
)

# Identify marker genes for cell clusters that show significantly different expression levels in one group compared to others.
scran_markers <- findMarkers_one_vs_all(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "leiden_clus",
                                        min_feats = 5)
topgenes_scran <- scran_markers[, head(.SD, 1), by = "cluster"]$feats

# Heatmap of the marker genes for the different leiden clusters (Figure 2D).
plotMetaDataHeatmap(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                    selected_feats = unique(topgenes_scran),
                    metadata_cols = "leiden_clus",
                    x_text_size = 10, y_text_size = 10)

# Heatmap of the marker genes for the different larvae (Figure 2E).
plotMetaDataHeatmap(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                    selected_feats = unique(topgenes_scran),
                    metadata_cols = "list_ID",
                    x_text_size = 10, y_text_size = 10)


# UMAPs plots of the different marker genes coloring cells by expression (Figure 2F).
dimFeatPlot2D(combined4fish4dpf_vs_2fish7dpf_fishobject, 
              expression_values = "scaled",
              feats = sort(unique(topgenes_scran)),
              cow_n_col = 2, 
              point_size = 3,
              save_param = list(base_width = 20, base_height = 20))

# Create spatial networks based on kNN method.
combined4fish4dpf_vs_2fish7dpf_fishobject <- createSpatialNetwork(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
                                                                  method = "kNN", 
                                                                  k = 6,
                                                                  maximum_distance_knn = 400,
                                                                  name = "spatial_network")

# Spatial plot coloring elavl3 gene expression of different larvae with the network lines in magenta (Figure 2G' and H').
spatFeatPlot2D(combined4fish4dpf_vs_2fish7dpf_fishobject,
               feats = "4-R elavl3 CH11", 
               expression_values = "scaled",
               #also I used normalized
               show_network= TRUE,
               spatial_network_name = "spatial_network",
               network_color = "magenta1",
               point_size = 5,
               cow_n_col = 2,
               group_by = "list_ID",
)

# Cell proximity enrichment analysis and its plot based on the leiden clustering (Figure 2K).
cell_proximities <- cellProximityEnrichment(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject,
                                            cluster_column = "leiden_clus",
                                            spatial_network_name = "spatial_network",
                                            number_of_simulations = 1000)
cellProximityBarplot(gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, 
                     CPscore = cell_proximities, 
                     min_orig_ints = 25, 
                     min_sim_ints = 25,
)

# Cell proximity network visualization based on the Cell proximity enrichment analysis (Figure 2I).
cellProximityNetwork (gobject = combined4fish4dpf_vs_2fish7dpf_fishobject, 
                      CPscore = cell_proximities, 
                      remove_self_edges = FALSE, 
                      only_show_enrichment_edges = TRUE,
                      self_loop_strength = 0.3
)

# Rank the genes on the spatial dataset depending on whether they exhibit a spatial pattern location or not.
ranktest <- binSpect(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                     calc_hub = TRUE,
                     spatial_network_name = "spatial_network",
                     implementation = c("simple"),
                     hub_min_int = 5,
                     bin_method = "rank")

# Do a heatmap of the top 6 spatial genes that each leiden cluster has (Figure 4B).
plotMetaDataHeatmap(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                    selected_feats = ranktest$feats[1:6],
                    metadata_cols = "leiden_clus",
                    #gradient_style = c("sequential"),
                    x_text_size = 10, y_text_size = 10)

# Do a heatmap of the top 6 spatial genes that each larvae has (Figure 4D).
plotMetaDataHeatmap(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                    selected_feats = ranktest$feats[1:6],
                    metadata_cols = "list_ID",
                    #gradient_style = c("sequential"),
                    x_text_size = 10, y_text_size = 10)

# Spatial co-expression modules. Cluster based on spatially correlated features.
ext_spatial_genes <- ranktest[1:12,]$feats

spat_cor_netw_DT <- detectSpatialCorFeats ( combined4fish4dpf_vs_2fish7dpf_fishobject,
                                            method = "network",
                                            spatial_network_name = "spatial_network",
                                            subset_feats = ext_spatial_genes)

# Cluster spatial genes by defining 5 groups and do heat map (Figure 4A).
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 5)


heatmSpatialCorFeats(combined4fish4dpf_vs_2fish7dpf_fishobject,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     show_cluster_annot = TRUE,
                     show_column_names = TRUE,
                     show_row_names = TRUE,
                     heatmap_legend_param = list(title = NULL))

# Create an spatial and UMAP plots for the top 6 spatial genes colored by gene expression (Figure 4E).
spatDimFeatPlot2D(combined4fish4dpf_vs_2fish7dpf_fishobject, 
                 expression_values = "scaled",
                 feats = ranktest$feats[1:6], 
                 plot_alignment = "horizontal",
                 #gradient_style = c("sequential"),
                 spat_point_size = 3,
                 dim_point_size = 3,
                 cow_n_col = 1,
)

# Generate an spatial plot of etv1, ret and hoxb5b colored by gene expression (Figure 4F).
spatFeatPlot2D(combined4fish4dpf_vs_2fish7dpf_fishobject,
               feats = c(
                 "2-R etv1 CH5",
                 "1-R ret  CH3",
                 "3-R hoxb5b CH8"),
               expression_values = "scaled",
               show_network= TRUE,
               spatial_network_name = "spatial_network",
               network_color = "black",
               point_size = 3,
               cow_n_col = 1,
)



