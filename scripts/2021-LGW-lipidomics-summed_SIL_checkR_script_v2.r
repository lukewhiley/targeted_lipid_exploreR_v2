###################################
###### summed SIL checkR #########
##################################

#this script fits into the lipidomics lipid exploreR. 
# the script will look for samples that have been got preparation errors
# the script sums the total intensity from SIL internal standard
# Samples are removed if the summed intensity is either > or < x from the median summed SIL internal standard signal of the dataset. 
# users define the SD cuttoff in custom mode (x)
# default is 50%
# Sample SIL intensity > x indicates high concentration of SIL, so likely as a result of excess SIL IS added
# Sample SIL intensity < x  below mean indicates too little volume added of SIL internal standard added
# users can select if failed samples are removed in custom mode
# automatically removed in default mode

dlg_message("Internal standard check. This next step will assess the internal standards accross all of the samples. If internal standards have been incorrectly added the summed signal intensity will be too low/high.", type = 'ok')

total_summed_sil <- apply(lipid_exploreR_data[["individual_lipid_data_unprocessed"]] %>% select(sampleID), 1, function(summedSIL){
  #browser()
  temp_data <- lipid_exploreR_data[["individual_lipid_data_unprocessed"]] %>% 
    filter(sampleID == summedSIL) %>% 
    select(-sampleID) %>% 
    select(contains("SIL")) %>% 
    rowSums(na.rm = TRUE)
}) %>% 
  c() %>% 
  as_tibble() %>%  
  add_column(lipid_exploreR_data[["individual_lipid_data_unprocessed"]]$sampleID, .before = 1) %>% 
  rename(SIL_TIC = value, sampleID = `lipid_exploreR_data[["individual_lipid_data_unprocessed"]]$sampleID`)

# arrange samples into correct run order, remove blanks and conditioning runs
total_summed_sil <- new_project_run_order %>% 
  left_join(total_summed_sil, by = "sampleID") %>% 
  arrange(injection_order) %>%
  filter(!grepl("blank", sampleID)) %>%
  filter(!grepl("COND", sampleID)) %>%
  filter(!grepl("conditioning", sampleID))

total_summed_sil$sample_idx <- c(1:nrow(total_summed_sil))
total_summed_sil$LOG_SIL_TIC <- log(total_summed_sil$SIL_TIC)

# while loop here
sil_check_status <- "change"
while(sil_check_status == "change"){

# #flag samples with SIL x standard deviations below mean
  temp_answer <- "blank"
  
  if(workflow_choice == "default"){
    temp_answer <- 75
  }

 while(is.na(as.numeric(temp_answer))){
   temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x % from the median", "e.g. recommended default x = 50")$res
 }
 
  median_sil_tic <- median(total_summed_sil$SIL_TIC)

  #sil_cut_off_lower <- median_sil_tic - (median_sil_tic*as.numeric(temp_answer)/100)
  #sil_cut_off_upper <- median_sil_tic + (median_sil_tic*as.numeric(temp_answer)/100)
  
  sil_cut_off_lower <- boxplot.stats(x = log(total_summed_sil$SIL_TIC), coef = 3)$stats[1]
  sil_cut_off_upper <- boxplot.stats(x = log(total_summed_sil$SIL_TIC), coef = 3)$stats[5]

#create lists of which samples have failed the SIL internal standard check
  #for both upper and lower filter
  
sil_qc_fail <- total_summed_sil$sampleID[which(log(total_summed_sil$SIL_TIC) < sil_cut_off_lower | log(total_summed_sil$SIL_TIC) > sil_cut_off_upper)] %>% as_tibble %>% rename(sampleID = value)
  #just lower filter
  #sil_qc_fail <- total_summed_sil$sampleID[which(log(total_summed_sil$SIL_TIC) < sil_cut_off_lower)] %>% as_tibble %>% rename(sampleID = value)
sil_qc_fail$fail_point <- "sil"
sil_qc_fail_ltr <- sil_qc_fail %>% filter(grepl(paste0(qc_type), sampleID))
sil_qc_fail_samples <- sil_qc_fail %>% filter(!grepl(paste0(qc_type), sampleID))

# visualise for reports
total_summed_sil$outlier <- "pass_qc"
total_summed_sil$outlier[total_summed_sil$sampleID %in% sil_qc_fail$sampleID] <- "outlier"

total_summed_sil_outlier <- total_summed_sil %>% filter(grepl("outlier", outlier))
total_summed_sil_pass <- total_summed_sil %>% filter(grepl("pass_qc", outlier))

# create a plate list ID
plate_number <- unique(plateID) %>% substr(14,14) %>% unique()
plate_number <- seq(1, length(plate_number))

plateIDx <- lapply(unique(total_summed_sil$plateID), function(FUNC_plateID){
  #browser()
  grep(FUNC_plateID, total_summed_sil$plateID)[1]}) %>% unlist()

#y_limit_lower <- min(log(total_summed_sil$SIL_TIC))-(min(log(total_summed_sil$SIL_TIC))/100*25)
#y_limit_upper <- max(log(total_summed_sil$SIL_TIC))+(max(log(total_summed_sil$SIL_TIC))/100*25)

# #set y axis limits
# if(sil_cut_off_lower < min(total_summed_sil$SIL_TIC)){
#   #y_limit_lower <- log(sil_cut_off_lower-(sil_cut_off_lower/100*25))
#   y_limit_lower <- min(total_summed_sil$LOG_SIL_TIC)-(total_summed_sil$LOG_SIL_TIC/100*25)
# }
# if(sil_cut_off_lower > min(total_summed_sil$SIL_TIC)){
#   #y_limit_lower <- log(min(total_summed_sil$SIL_TIC)-(min(total_summed_sil$SIL_TIC)/100*25))
#   y_limit_lower <- min(total_summed_sil$SIL_TIC)-(min(total_summed_sil$SIL_TIC)/100*25)
# }
# if(sil_cut_off_upper > max(total_summed_sil$SIL_TIC)){
#   #y_limit_upper <- log(max(total_summed_sil$SIL_TIC)+(max(total_summed_sil$SIL_TIC)/100*25))
#   y_limit_upper <- max(log(total_summed_sil$SIL_TIC))+(max(log(total_summed_sil$SIL_TIC))/100*25)
# }
# if(sil_cut_off_upper < max(total_summed_sil$SIL_TIC)){
#   #y_limit_upper <- log(max(total_summed_sil$SIL_TIC)+(max(total_summed_sil$SIL_TIC)/100*25))
#   y_limit_upper <- max(log(total_summed_sil$SIL_TIC))+(max(log(total_summed_sil$SIL_TIC))/100*25)
# }

y_limit_lower <- total_summed_sil$LOG_SIL_TIC %>% min() * 0.8
y_limit_upper <- total_summed_sil$LOG_SIL_TIC %>% max() * 1.1

# create a layout list of extra lines to add
p_threshold_lines <- list(list(type='line', x0= min(total_summed_sil$sample_idx), x1= (max(total_summed_sil$sample_idx)+10), y0=sil_cut_off_lower, y1=sil_cut_off_lower,
                          line=list(dash='dot', width=3, color = '#FF0000')),
                     list(type='line', x0= min(total_summed_sil$sample_idx), x1= (max(total_summed_sil$sample_idx)+10), y0=sil_cut_off_upper, y1=sil_cut_off_upper,
                          line=list(dash='dot', width=3, color = '#FF0000')),
                     list(type='line', x0= min(total_summed_sil$sample_idx), x1= (max(total_summed_sil$sample_idx)+10), y0=log(median_sil_tic), y1=log(median_sil_tic),
                          line=list(dash='dot', width=3, color = '#000000'))
)
p_plate_list <- lapply(plateIDx[2:length(plateIDx)], function(FUNC_P_PLATE_LIST){
   list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, y0=y_limit_lower-log(median_sil_tic)*.05, y1=y_limit_upper+log(median_sil_tic)*.05,
            line=list(dash='dot', width=2, color = '#808080'))
})

if(length(plateIDx) == 1){
p_plot_lines <- p_threshold_lines
}

if(length(plateIDx) > 1){
  p_plot_lines <- c(p_threshold_lines, p_plate_list)
}

#create a list of axis settings for plot_ly
x_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = FALSE,
  range = c(0, max(total_summed_sil$sample_idx)+10),
  title = "Sample index"
)

y_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = TRUE,
  title = "Lipid total ion count (Log)",
  range = c(y_limit_lower, 
            y_limit_upper)
)

p <- plot_ly(
  type = "scatter", mode = "markers", data = total_summed_sil_pass, x = ~sample_idx, y = ~LOG_SIL_TIC, text = ~sampleID, color = ~outlier, colors = c('#1E90FF', '#FF0000'), 
  marker = list(size = 7, color = '#1E90FF', opacity = 0.5,
                line = list(color = '#000000',width = 1))
 ) %>% 
  add_trace(type = "scatter", data = total_summed_sil_outlier, x = ~sample_idx, y = ~LOG_SIL_TIC, text = ~sampleID, color = ~outlier, 
            marker = list(size = 8, color = '#FF0000')
            ) %>%
  layout(xaxis = x_axis_settings,
         yaxis = y_axis_settings
         ) %>%
  layout(shapes=p_plot_lines)

#create html widget and display it in the users internet browser
sil_check_p <- p

saveWidget(sil_check_p, file = paste(project_dir_html, "/", project_name, "_", user_name, "_SIL_check_plot.html", sep=""))# save plotly widget
#browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_SIL_check_plot.html", sep="")) #open plotly widget in internet browser

#sil_qc_fail - ask the user if they wish to continue or change the threshold
sil_check_status <- "blank"

if(workflow_choice == "default"){
  sil_check_status <- "continue"
}

while(sil_check_status != "continue" & sil_check_status != "change"){
  sil_check_status <- dlgInput(paste(nrow(sil_qc_fail), "samples FAILED the SIL QC check.  continue or change the exclusion threshold?"), "continue/change")$res
}
}

#sil_qc_fail - ask the user if they wish to remove all/none/samples/LTR which failed the QC check
temp_answer <- "blank"

if(workflow_choice == "default"){
  temp_answer <- "all"
  dlg_message(paste0(nrow(sil_qc_fail), "  samples FAILED the SIL QC check  ", nrow(sil_qc_fail_ltr),"  were ", paste0(qc_type), ".  These have been removed from the dataset."), 
              type = 'ok')
}

while(temp_answer != "all" & temp_answer != "none" & temp_answer != "samples" & temp_answer != "LTR"){
  temp_answer <- dlgInput(paste("of the ", nrow(sil_qc_fail), "FAILED samples.  ",  nrow(sil_qc_fail_ltr), " were ", paste0(qc_type), ".  Do you want to remove failed samples?"), paste0("all/none/samples/", qc_type))$res
}

if(temp_answer == "all"){lipid_exploreR_data[["individual_lipid_data_sil_filtered"]] <- lipid_exploreR_data[["individual_lipid_data_unprocessed"]] %>% filter(!sampleID %in% sil_qc_fail$sampleID)}
if(temp_answer == "samples"){lipid_exploreR_data[["individual_lipid_data_sil_filtered"]] <- lipid_exploreR_data[["individual_lipid_data_unprocessed"]] %>% filter(!sampleID %in% sil_qc_fail_samples$sampleID)}
if(temp_answer == paste0(qc_type)){lipid_exploreR_data[["individual_lipid_data_sil_filtered"]] <- lipid_exploreR_data[["individual_lipid_data_unprocessed"]] %>% filter(!sampleID %in% sil_qc_fail_ltr$sampleID)}
if(temp_answer == "none"){lipid_exploreR_data[["individual_lipid_data_sil_filtered"]] <- lipid_exploreR_data[["individual_lipid_data_unprocessed"]]}

#filter out blanks and conditioning runs to take dataset forwards
lipid_exploreR_data[["individual_lipid_data_sil_filtered"]] <- lipid_exploreR_data[["individual_lipid_data_sil_filtered"]] %>% 
  filter(!grepl("blank", sampleID)) %>%
  filter(!grepl("COND", sampleID)) %>%
  filter(!grepl("conditioning", sampleID))

