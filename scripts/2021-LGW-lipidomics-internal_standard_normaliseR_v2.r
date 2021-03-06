################################################
###### Internal standard normalization #########
#################################################

## data preparation from previous script

data_for_sil_normalisation <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered

#replace NAs/missing values with a low value 
#first convert all 0 values to NA
data_for_sil_normalisation[data_for_sil_normalisation==0] <- NA

#then impute all NAs with the minimum value in the spreadsheet (estimate of a noise peak)

data_for_sil_normalisation[is.na(data_for_sil_normalisation)] <- data_for_sil_normalisation %>% 
  select(contains("(")) %>%
  as.matrix() %>%
  min(na.rm = TRUE)

# this section normalizes each SRM lipid target with the appropriate internal standard
# Requires a template guide with internal standard transition for each target lipid SRM transition

dlg_message("Time for normalization using the internal standards :-)", type = 'ok')

dlg_message("First we need to create response ratios by dividing the lipid target peak area by the peak area for the appropiate internal standard", type = 'ok')

dlg_message("Step 1 - please prepare a reference csv file listing each lipid target and the assigned internal standard. Column headings required are 'Precursor Name' containing the lipid target  and 'Note' containing the IS name", type = 'ok')

#import transition report 3
# filtered_data <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered  %>% 
#   filter(!grepl("conditioning", sampleID))

dlg_message("Please select this template file now.", type = 'ok')




temp_answer <- "change"
while(temp_answer == "change"){
sil_target_list <- read_csv(file = file.choose(.)) %>% clean_names

lipid_data_for_normalisation <- data_for_sil_normalisation %>% 
  select(sampleID, plateID, contains("(")) %>% 
  select(!contains("SIL"))

sil_data_for_normalisation <- data_for_sil_normalisation %>% 
  select(sampleID, plateID, contains("SIL")) #%>% filter(grepl("LTR", sampleID))

# this section checks each of the SIL IS used in the target list template in the LTRs. It evaluates if:
##  a: is the internal standard present in the LTR samples? Some batches of IS do not contain every IS availible. This alos prevents user error if the IS batch has not been made correctly.
##  b: the RSD of the internal standard signal intensity in LTRs. If the SIL IS is suitable for use in the dataset its signal should be stable in the LTRs. (raises a warning for >30%)

dlg_message("Checks to see if all internal standards are present in the SIL internal standard mix", type = 'ok')



#create a list of SIL internal standards used
sil_list <- sil_target_list %>% filter(grepl("SIL", note)) %>% select(note) %>% unique() 

# check for missing SIL values indicating SIL has not been correctly added
missing_sil_list <- NULL
for(idx_sil in 1:length(sil_list$note)){
  temp_sil_data <- sil_data_for_normalisation %>% select(sampleID, plateID, sil_list$note[idx_sil])
  colnames(temp_sil_data) <- c("sampleID", "plateID", "area")
  temp_sil_data$lipidID <- sil_list$note[idx_sil]
  #if(length(which(temp_sil_data[,3]==0))>0){
  missing_sil_list <- rbind(missing_sil_list,
                              temp_sil_data[which(temp_sil_data[,3]==0),])
  missing_sil_list <- rbind(missing_sil_list,
                              temp_sil_data[which(is.na(temp_sil_data[,3])),])
#  }
}

# create a warning list of missing SIL compounds
sil_list_warning <- missing_sil_list %>% select(sampleID, plateID, lipidID)
sil_list_warning$reason_for_flag <- "missing value in sample - check skyline"

#creates a summed sil value for all samples per SIL compound
#looks for SIL compounds that are well below the interquartile range for the SILS
# again indicates the SIL hasn't been added correctly

# sil_sum <- lapply(sil_list$note, function(FUNC_SIL){
#   #browser()
#   temp_func_data_sum <- sil_data_for_normalisation %>% 
#     select(all_of(FUNC_SIL)) %>% 
#     as.matrix() %>% 
#     sum() %>% 
#     log()
# }) %>% 
#   c() %>% 
#   unlist() %>% 
#   as_tibble() %>% 
#   rename(SIL_SUM = value) %>% 
#   add_column(sil_list, .before = 1) %>% 
#   arrange(SIL_SUM)
# 
# sil_sum_q1 <- quantile(sil_sum$SIL_SUM, 0.25, na.rm = TRUE) %>% as.numeric()
# inter_quantile_range <- as.numeric(quantile(sil_sum$SIL_SUM, 0.75, na.rm = TRUE)) - as.numeric(quantile(sil_sum$SIL_SUM, 0.25, na.rm = TRUE))
# sil_sum_lower_threshold <- sil_sum_q1 - inter_quantile_range
# 
# #create a list of IS that fail the test
# sil_list_warning_2 <- sil_sum$note[which(sil_sum$SIL_SUM < sil_sum_lower_threshold)]


# looks at SIL in LTR only for signal (area) %RSD across the LTRs
sil_data_check_ltr <- sil_data_for_normalisation %>% 
  filter(grepl(paste0(qc_type), sampleID))

sil_rsd <- lapply(sil_list$note, function(FUNC_SIL){
  #browser()
  temp_func_data <- sil_data_check_ltr %>% select(all_of(FUNC_SIL))
  temp_func_data_mean <- temp_func_data %>% as.matrix() %>% mean()
  temp_func_data_sd <- temp_func_data %>% as.matrix() %>% sd()
  temp_func_data_rsd <- (temp_func_data_sd/temp_func_data_mean)*100
  temp_func_data_rsd
}) %>% 
  c() %>% 
  unlist() %>% 
  as_tibble() %>% 
  rename(SIL_RSD = value) %>% 
  add_column(sil_list, .before = 1)

sil_list_warning_3 <- sil_rsd %>% filter(SIL_RSD > 30) %>% rename(lipidID = note) %>% select(lipidID)
sil_list_warning_3$plateID <- NA 
sil_list_warning_3$reason_for_flag <- "RSD > 30% across dataset"
sil_list_warning_3$sampleID <- NA
sil_list_warning_3 <- sil_list_warning_3 %>% select(sampleID, plateID, lipidID, plateID, lipidID, reason_for_flag)

sil_list_warning <- sapply(sil_list_warning, as.character) %>% as_tibble()
sil_list_warning_3 <- sapply(sil_list_warning_3, as.character) %>% as_tibble()

sil_list_warning <- bind_rows(sil_list_warning,
                              sil_list_warning_3)

sil_list_warning_html <- htmlTable(sil_list_warning)

htmltools::save_html(sil_list_warning_html, file = paste(project_dir_html, "/", project_name, "_", user_name, "_SIL_internal_standard_fail_list.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_SIL_internal_standard_fail_list.html", sep="")) #open plotly widget in internet browser


dlg_message(paste("SIL WARNING: CHECK BROWSER FOR TABLE OF FAILED SIL SRTANDARDS AND CHECK SKYLINE IF NECESSARY"), 
            type = 'ok')

temp_answer <- "blank"
if(workflow_choice == "default"){
  temp_answer <- "continue"
}

#add in a check in case the user enters the incorrect entry. It must be "continue" or "change" to continue
while(temp_answer != "continue" & temp_answer!="change"){
  temp_answer <- dlgInput("Do you wish to continue or use different internal standards?", "continue/change")$res
}

# if change has been selected the user now has an opportunity to change the internal standard import file template
if(temp_answer == "change"){
  dlg_message("OK - please edit and select a new internal standard reference file now")
}
}

ratio_concentration_choice <- "blank"
if(workflow_choice == "default"){
  ratio_concentration_choice <- "concentration"
}
while(ratio_concentration_choice != "ratio" & ratio_concentration_choice != "concentration"){
  ratio_concentration_choice <- dlgInput("Do you also want to estimate the concentrations or just continue to use the response ratio?", "ratio/concentration")$res
}


if(ratio_concentration_choice == "concentration"){

# now multiply by the IS concentration to create a concentration factor
dlg_message("We also need to calculate lipid concentrations from the internal standard. Please select the concentration template csv - NOTE: use the correct lot for your analysis", type = 'ok')

#read in SIL-internal standard concentration file. aim to put on github to centralise
sil_concentrations <- read_csv(file = file.choose(.)) %>% clean_names

sil_batch <- "blank"
while(is.na(as.numeric(sil_batch))){
  sil_batch <- dlgInput("What batch did you use?", "101/102/103")$res
  sil_batch <- as.numeric(sil_batch)
}
}


# this apply function creates:
# 1. a response ratio by dividing the signal area for each target lipid by the peak area from the appropriate SIL IS metabolite. 
# 2. a final estimated concentration using the pre-defined internal standard as a single point calibration

ratio_data <- lapply(lipid_data_for_normalisation %>% select(contains("(")) %>% names(), function(FUNC_IS_RATIO){
  #browser()
  # step 1 - create a ratio between lipid target and the appropriate internal standard as pre-defined in the reference file
  func_data <- data_for_sil_normalisation %>%
    select(sampleID, all_of(FUNC_IS_RATIO))
  
  sil_to_use <- sil_target_list$note[which(sil_target_list$precursor_name==FUNC_IS_RATIO)]
  
  func_data_sil <- sil_data_for_normalisation %>% select(sampleID, all_of(sil_to_use))
  normalised_data <- func_data %>% left_join(func_data_sil, by = "sampleID")
  normalised_data$ratio <- (normalised_data %>% select(all_of(FUNC_IS_RATIO)) %>% as.matrix()/normalised_data %>% select(all_of(sil_to_use)) %>% as.matrix()) %>% 
    as.numeric()
  
  normalised_data <- normalised_data %>% select(sampleID, ratio) 
  names(normalised_data) <- c("sampleID", FUNC_IS_RATIO)
  
  concentration_data <- normalised_data
 
  if(ratio_concentration_choice == "ratio"){
    normalised_data
  }
  
  #step 2 - select concentration factor from csv template and multiply normalised data by concentration factor
  if(ratio_concentration_choice == "concentration"){
  func_concentration_factor <- sil_concentrations %>% 
    filter(sil_name == sil_to_use) %>% 
    select(concentration_factor) %>% 
    as.numeric()
  
  concentration_data[,2] <- concentration_data[,2]*func_concentration_factor
  concentration_data
  }
  
}) %>% 
  reduce(left_join, by = "sampleID")

ratio_data <- data_for_sil_normalisation %>%
  select(sampleID, plateID) %>% 
  left_join(ratio_data, by = "sampleID")

lipid_exploreR_data[["individual_lipid_data_sil_tic_intensity_filtered_ratio"]] <- ratio_data

