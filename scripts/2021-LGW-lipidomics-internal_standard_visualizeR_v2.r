#visualise the RSD of the ltrs

ltr_rsd <- apply(as_tibble(colnames(ratio_data %>% select(contains("(")))), 1, function(FUNC_LTR_RSD){
  #browser()
  func_data <- ratio_data %>% filter(grepl(paste0(qc_type), sampleID)) %>% select(all_of(FUNC_LTR_RSD))
  (sd(func_data$value)*100)/mean(func_data$value)
}) %>% as_tibble() %>% add_column(colnames(ratio_data %>% select(contains("("))), .before = 1)

colnames(ltr_rsd) <- c("lipid", "RSD")

temp_answer <- "blank"

if(workflow_choice == "default"){
  temp_answer <- "yes"
  }

while(temp_answer != "yes" & temp_answer != "no"){
  temp_answer <- dlgInput(paste(nrow(ltr_rsd)-length(which(ltr_rsd$RSD < 30)), " lipid targets had an RSD in ", paste0(qc_type), " of > 30%. Do you wnt to remove them?", sep=""), "yes/no")$res
  }


if(temp_answer == "yes"){
  dlg_message(paste(nrow(ltr_rsd)-length(which(ltr_rsd$RSD < 30)), " lipid targets had an RSD in ", paste0(qc_type), " of > 30% so were removed from the dataset", sep=""))
  lipid_keep_list <- ltr_rsd %>% filter(RSD < 30)
}

if(temp_answer == "no"){
  lipid_keep_list <- ltr_rsd
}

#dlg_message(paste(nrow(ltr_rsd)-length(which(ltr_rsd$RSD < 30)), " lipid targets had an RSD in ", paste0(qc_type), " of > 30% so were removed from the dataset", sep=""))
dlg_message(paste("number of feature ratios with with an RSD in ", paste0(qc_type), " of <30% =", length(which(ltr_rsd$RSD < 30))), type = 'ok')
dlg_message(paste("number of feature ratios with with an RSD in ", paste0(qc_type), " of <20% =", length(which(ltr_rsd$RSD < 20))), type = 'ok')
dlg_message(paste("number of feature ratios with with an RSD in ", paste0(qc_type), " of <15% =", length(which(ltr_rsd$RSD < 15))), type = 'ok')
dlg_message(paste("number of feature ratios with with an RSD in ", paste0(qc_type), " of <10% =", length(which(ltr_rsd$RSD < 10))), type = 'ok')

#lipid_keep_list <- ltr_rsd %>% filter(RSD < 30)


rsd_filtered_data <- ratio_data %>% select(sampleID, plateID, all_of(lipid_keep_list$lipid))
#rsd_filtered_data[is.na(rsd_filtered_data)] <- 0

#non_filtered_dataset <- ratio_data_2 %>% select(-plateID)
#non_filtered_dataset[is.na(non_filtered_dataset)] <- 0

# visualisation of normalized data
# first - produce a plot of all normalized features to see if there are any overall trends in the data

total_summed_ratio <- apply(rsd_filtered_data %>% select(sampleID), 1, function(summedTIC){
  #browser()
  temp_data <- rsd_filtered_data %>% filter(sampleID == summedTIC) %>% select(-sampleID, -plateID) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(rsd_filtered_data$sampleID, rsd_filtered_data$plateID, .before = 1) %>% 
  rename(summed_TIC = value, sampleID = "rsd_filtered_data$sampleID", plateID = "rsd_filtered_data$plateID")
total_summed_ratio$sample_idx <- c(1:nrow(total_summed_ratio))

#sd(total_summed_ratio$summed_TIC*100)/mean(total_summed_ratio$summed_TIC)

total_summed_ratio$sample <- "sample"
total_summed_ratio$sample[grep(paste0(qc_type), total_summed_ratio$sampleID)] <- paste0(qc_type)

# create a plotly plot to visualize
total_summed_ratio$log_summed_TIC <- log(total_summed_ratio$summed_TIC+1)

total_summed_ratio_samples <- total_summed_ratio %>% filter(!grepl(paste0(qc_type), sampleID))
total_summed_ratio_LTR <- total_summed_ratio %>% filter(grepl(paste0(qc_type), sampleID))

# create a plate list ID
plate_number <- unique(plateID) 
plateIDx <- lapply(unique(plateID), function(FUNC_plateID){
  #browser()
  grep(FUNC_plateID, total_summed_ratio$plateID)[1]}) %>% unlist()
  

# create a layout list of extra lines to add
p_plate_list <- lapply(plateIDx[2:length(plateIDx)], function(FUNC_P_PLATE_LIST){
  list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, 
       y0=log(min(total_summed_ratio_samples$summed_TIC)-(min(total_summed_ratio_samples$summed_TIC)/100*50)), 
       y1=log(max(total_summed_ratio_samples$summed_TIC)+(max(total_summed_ratio_samples$summed_TIC)/100*25)),
       line=list(dash='dot', width=2, color = '#808080'))
})

#only add plate lines if multiple plates exist
# if(is.na(plateIDx)){
#   p_plot_lines <- NULL
# }

if(length(plateIDx) == 1){
  p_plot_lines <- NULL
}

if(length(plateIDx) > 1){
  p_plot_lines <- c(p_plate_list)
} 

#create a list of axis settings for plot_ly
x_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = FALSE,
  range = c(0, max(total_summed_ratio_samples$sample_idx)),
  title = "Sample index"
)

y_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = TRUE,
  range = c(log(min(total_summed_ratio_samples$summed_TIC)-(min(total_summed_ratio_samples$summed_TIC)/100*50)), 
            log(max(total_summed_ratio_samples$summed_TIC)+(max(total_summed_ratio_samples$summed_TIC)/100*25))),
  title = "Summed lipid target response ratio (Log)"
)

  p <- plot_ly(
    type = "scatter", mode = "markers", data = total_summed_ratio_samples, x = ~sample_idx, y = ~log_summed_TIC, text = ~sampleID, color = ~sample, colors = c('#1E90FF', '#FF0000'), 
    marker = list(size = 7, color = '#1E90FF', opacity = 0.5,
                  line = list(color = '#000000',width = 1))
  ) %>% 
    add_trace(type = "scatter", data = total_summed_ratio_LTR, x = ~sample_idx, y = ~log_summed_TIC, text = ~sampleID, color = ~sample, 
              marker = list(size = 8, color = '#FF0000')
    ) %>%
    layout(xaxis = x_axis_settings,
           yaxis = y_axis_settings
    ) %>%
    layout(shapes=p_plot_lines) 

normalized_check_p <- p

saveWidget(normalized_check_p, file = paste(project_dir_html, "/", project_name, "_", user_name, "_normalized_check_plot.html", sep=""))# save plotly widget
#browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_normalized_check_plot.html", sep="")) #open plotly widget in internet browser

#dlg_message("Check plot for summed all normilized features. Press OK to continue", type = 'ok')

# second - produce a plot of normalized features, summed by class to see if there are any overall trends in the lipid class data

#dlg_message("Now we are going to look at the data summed by lipid class", type = 'ok')

rsd_filtered_class_data <- create_lipid_class_data_summed(rsd_filtered_data)

lipid_class_list <- rsd_filtered_data %>% select(contains("(")) %>% colnames() 
lipid_class_list <- sub("\\(.*", "", lipid_class_list) %>% unique()
lipid_class_list <- lipid_class_list[!grepl("sampleID", lipid_class_list)] %>% as_tibble()

#add ltr TRUE/FALSE column

rsd_filtered_class_data$is_qc <- "sample"
rsd_filtered_class_data$is_qc[grep(paste0(qc_type), rsd_filtered_class_data$sampleID)] <- paste0(qc_type)

plotlist <- apply(lipid_class_list %>% select(value), 1, function(lipidClass){
  #browser()
  plot_data <- rsd_filtered_class_data %>% select("sampleID", "is_qc", all_of(lipidClass)) %>% rename(ms_response = value) 
  plate_id <- str_extract(plot_data$sampleID, "PLIP.*")
  plate_id <- substr(plate_id, 0,15)
  plot_data$sample_index <- paste(plate_id, sub(".*\\_", "", plot_data$sampleID), sep="_")
  plot_data$sample_idx <- 1:nrow(plot_data)
  
  plot_data_ltr <- plot_data %>% filter(is_qc == paste0(qc_type))
  plot_data <- plot_data %>% filter(is_qc == "sample")
  
  # plateIDx <- lapply(unique(plateID), function(FUNC_plateID){
  #   #browser()
  #   grep(FUNC_plateID, plot_data$plateID)[1]}) %>% unlist()
  
  # create a plate list ID
  plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
  plateIDx <- lapply(unique(plateID), function(FUNC_plateID){
    #browser()
    grep(FUNC_plateID, total_summed_ratio$plateID)[1]}) %>% unlist()
  
  # create a layout list of extra lines to add
  p_plate_list <- lapply(plateIDx[2:length(plateIDx)], function(FUNC_P_PLATE_LIST){
    list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, 
         y0=(log(min(plot_data$ms_response)+1)-(log(min(plot_data$ms_response)+1)/100*50)), 
         y1=(log(max(plot_data$ms_response)+1)+(log(max(plot_data$ms_response)+1)/100*25)),
         line=list(dash='dot', width=2, color = '#808080'))
  })
  
  #only add plate lines if multiple plates exist
  # if(is.na(plateIDx)){
  #   p_plot_lines <- NULL
  # }
  
  if(length(plateIDx) == 1){
    p_plot_lines <- NULL
  }
  
  if(length(plateIDx) > 1){
    p_plot_lines <- c(p_plate_list)
  } 
  
  xmax <- max(plot_data$sample_idx)
    
  #create a list of axis settings for plot_ly
  x_axis_settings <- list(
    zeroline = FALSE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = FALSE,
    range = c(0, xmax+1),
    title=""
    #title = paste(lipidClass)
  )
  
  ymin <- log(min(plot_data$ms_response)+1)-(log(min(plot_data$ms_response)+1)/100*50)
  ymax <- log(max(plot_data$ms_response)+1)+(log(max(plot_data$ms_response)+1)/100*25)
  
  y_axis_settings <- list(
    zeroline = FALSE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    range = c(ymin, (ymax*1.1)),
    title = "Summed lipid target/internal standard ratio (Log)"
  )
  
  #browser()
  p <- plot_ly(
    type = "scatter", mode = "markers",  colors = c('#1E90FF', '#FF0000'), data = plot_data, x = ~sample_idx, y = ~log(ms_response+1), text = ~sampleID, color = ~is_qc, 
    marker = list(size = 5, color = '#1E90FF', opacity = 0.5,
                  line = list(color = '#000000',width = 1)),
    showlegend = FALSE
  ) %>% 
    add_trace(type = "scatter", data = plot_data_ltr, x = ~sample_idx, y = ~log(ms_response+1), text = ~sampleID, color = ~is_qc, 
              marker = list(size = 6, color = '#FF0000'),
              showlegend = FALSE
    ) %>%
    layout(xaxis = x_axis_settings,
           yaxis = y_axis_settings
    ) %>%
    add_annotations(x = xmax/2, 
                    y = ymax*1.1,
                    text = paste(lipidClass), showarrow = F) %>%
    layout(shapes=p_plot_lines)
    
  p
})


normalized_check_class_p <- subplot(plotlist, nrows = 4, titleX = TRUE, margin = c(0.015,0.015, 0.05,0.05))

saveWidget(normalized_check_class_p, file = paste(project_dir_html, "/", project_name, "_", user_name, "_normalized_check_class_plot.html", sep=""))# save plotly widget
#browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_normalized_check_class_plot.html", sep="")) #open plot_ly widget in internet browser

#dlg_message("Check plot for summed lipid class normilized features. Press OK to continue", type = 'ok')

