#ANPC Lipidomics PCA quality control visualisation

# FUNC_individual_multivariate_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lipids_pca <- function(FUNC_individual_multivariate_data, FUNC_colour_by, FUNC_plot_label, FUNC_scaling){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  #browser()
  
  lipid <- FUNC_individual_multivariate_data %>% select(contains("(")) %>% colnames()
  
  #create data matrix for PCA
  pca_x <- FUNC_individual_multivariate_data %>%  select(all_of(lipid)) %>% as.matrix()+1 
  pca_x <- log(pca_x)
  title_text <- "individual lipid species"
  pca_x[pca_x == 0] <- NA #remove all 0 values
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
  
  #create PCA model
  pca_model <- pca(pca_x, pc=2, scale = paste(FUNC_scaling), center = TRUE)
  
  # extract score values for plotting in plot_ly
  PC1 <- as.numeric(as.matrix(pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_model@t[,2]))
  
  # extract loadings values for plotting in plot_ly
  plotly_loadings_data <- pca_model@p %>% as_tibble(rownames = "lipid") %>% rename(PC1 = V1, PC2 = V2)

  #produce plot_ly PCA scores plot
  
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  pca_colour <- FUNC_individual_multivariate_data %>% select(all_of(FUNC_colour_by)) %>% as.matrix()
  pca_colour[is.na(pca_colour)] <- "none"
  pca_plot_label <- FUNC_individual_multivariate_data %>% 
    select(all_of(FUNC_plot_label)) %>% 
    as.matrix()
  
  # create plot values
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$pca_colour <- c(pca_colour)
  plot_Val$pca_plot_label <- c(pca_plot_label)
  
  plot_colors <- RColorBrewer::brewer.pal(name = "Set3",
                                          n = length(unique(pca_colour)))
  
  x_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("PC1 (", round(pca_model@Parameters$R2[1]*100,1), " %)", sep = "")
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("PC2 (", round(pca_model@Parameters$R2[2]*100,1), " %)", sep = "")
  )
  
 plotly_pca <- plot_ly(type = "scatter", 
                       mode = "markers", 
                       data = plot_Val, 
                       x = ~PC1, 
                       y = ~PC2, 
                       text = ~pca_plot_label, 
                       color = ~pca_colour, 
                       colors = c(plot_colors[1:length(unique(pca_colour))]), 
                        marker = list(size = 7, 
                                      #color = '#1E90FF', 
                                      opacity = 0.5,
                                      line = list(
                                        color = '#000000',
                                        width = 1)
                        )) %>% 
    layout(title = paste(" Plotly PCA - ", title_text, sep = ""),
           xaxis = x_axis_settings_scores,
           yaxis = y_axis_settings_scores)
  
 
 
 # create loadings plot
  x_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("")
  )
  
  y_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("")
  )
  
  plotly_loadings <- plot_ly(type = "scatter", 
                             mode = "markers", 
                             data = plotly_loadings_data, 
                             x = ~PC1, 
                             y = ~PC2, 
                             text = ~lipid, 
                             marker = list(size = 7, color = '#808080', opacity = 0.5,
                                           line = list(color = '#000000', width = 1)
                             )) %>% 
    layout(title = paste(" Plotly PCA - ", title_text, sep = ""),
           xaxis = x_axis_settings_loading,
           yaxis = y_axis_settings_loading
    )
  
  combined_plotly <- subplot(plotly_pca, plotly_loadings, 
                             margin = c(0.05, 0.05, 0.01, 0.01),
                             titleX = TRUE,
                             titleY = TRUE
  ) %>% layout(showlegend = TRUE, title =  "")
  

  combined_plotly

}