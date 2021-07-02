
#correlation check between  

corr_out = NULL

for(idx in unique(sil_target_list$precursor_name)){
  loop_data <- master_lipid_data %>%
    filter(lipid_target == idx) %>%
    rename(sampleID = replicate) %>%
    select(sampleID, area)
  
  is_used <- sil_target_list$note[which(sil_target_list$precursor_name == idx)]
  
  loop_is <- master_lipid_data %>%
    filter(lipid_target == is_used) %>%
    rename(sampleID = replicate, is_area = area) %>%
    select(sampleID, is_area)
  
  checkpoint <- which(names(old_data) == idx)
  
  if(length(checkpoint) > 0){
  
  processed_data <- final_corrected_data %>% 
    select(sampleID, idx)
  
  test_data <- processed_data %>%
    left_join(loop_data, by = "sampleID") %>%
    left_join(loop_is, by = "sampleID")
  test_data$normalised <- test_data$area/test_data$is_area
  
  corr_result <- cor(test_data[,2], test_data$normalised)
  
  corr_out <- rbind(corr_out, 
                    corr_result)
  
  }
}


