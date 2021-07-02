# ######### read in data etc
# 
# # # load required packages not required if running as part of the notebook
# 
# package_list <- c("plyr", "tidyverse", "janitor", "gridExtra", "ggpubr", "readxl", "cowplot", "scales", "stats", "devtools", "metabom8", "shiny", "plotly", "svDialogs", "DataEditR", "htmlwidgets", "httr", "htmlTable", 'statTarget')
# loaded_packages <- lapply(package_list, require, character.only = TRUE)
# rm(loaded_packages, package_list)
# 
# # load custom functions from github
# lipidomics_class_sum_function <- GET(url = "https://raw.githubusercontent.com/lukewhiley/targeted_lipid_exploreR/main/ANPC_lipidomics_tools/functions/2021-LGW-lipidomics-class_sumR_function.r") %>% content(as = "text")
# eval(parse(text = lipidomics_class_sum_function), envir = .GlobalEnv)
# rm(lipidomics_class_sum_function)
# 
# dlg_message("Welcome to lipid exploreR! :-)", type = 'ok')
# 

dlg_message("Please select your project folder", type = 'ok'); project_dir <- rstudioapi::selectDirectory() # save project directory root location
setwd(project_dir) # switch the project directory

# create a new directory to store html widgets
if(!dir.exists(paste(project_dir, paste("/",Sys.Date(), "_html_files", sep=""), sep=""))){
  dir.create(paste(project_dir, paste("/",Sys.Date(), "_html_files", sep=""), sep=""))
} 
project_dir_html <- paste(project_dir, paste("/",Sys.Date(), "_html_files", sep=""), sep="")

#user input here for project name and user initials
temp_answer <- "blank"
if(exists("project_name") == TRUE){temp_answer <- dlgInput(paste("the project name is ", project_name, "is this correct?", sep=" "), "yes/no")$res}
while(temp_answer != "yes"){
project_name <- dlgInput("what is the name of the project?", "example_project")$res
temp_answer <- dlgInput(paste("the project name is ", project_name, "is this correct?", sep=" "), "yes/no")$res
}

temp_answer <- "blank"
if(exists("user_name") == TRUE){temp_answer <- dlgInput(paste("the user is ", user_name, "is this correct?", sep=" "), "yes/no")$res}
while(temp_answer != "yes"){
  user_name <- dlgInput("Insert your initials", "example_initials")$res
  temp_answer <- dlgInput(paste("the user is ", user_name, "is this correct?", sep=" "), "yes/no")$res
}

# read in master data
temp_answer <- "blank"
while(temp_answer != "yes" & temp_answer != "no"){
temp_answer <- dlgInput("Do you want to read in new data?", "yes/no")$res
}
if(temp_answer == "yes"){dlg_message("Open metabolite data csv by selecting the CSV generated by the SkylineR data processing script", type = 'ok')
  master_lipid_data <- read_csv(file = file.choose(.)) %>% clean_names %>% rename(lipid_target = peptide, lipid_class = protein)
}

sampleID <- master_lipid_data$replicate %>% unique() # create list of sample IDs
lipid <- master_lipid_data$lipid_target %>% unique() # create list of lipid targets
lipid_class_list <- master_lipid_data %>% select(lipid_class) %>% unique()

individual_lipid_data <- apply(as_tibble(lipid), 1, function(lip){
  #browser()
  sampleID <- master_lipid_data$replicate %>% unique() %>% as_tibble() # create list of sample IDs
  temp_data <- master_lipid_data %>% filter(lipid_target == lip) %>% select(replicate, area)
  colnames(temp_data) <- c("value", lip) 
  temp_data <- left_join(sampleID, temp_data, by = "value") 
  #temp_data <- temp_data %>% select(lip)
}) %>% bind_cols() %>% select(all_of(lipid)) %>% add_column(sampleID, .before = 1)

plateID <- str_extract(individual_lipid_data$sampleID, "PLIP.*")
plateID <- substr(plateID, 0,15)
plateID <- paste(plateID, sub(".*\\_", "", individual_lipid_data$sampleID), sep="_")

individual_lipid_data <- individual_lipid_data %>% add_column(plateID, .before = 2) %>% arrange(plateID)


project_run_order <- individual_lipid_data %>% select(sampleID, plateID)
project_run_order$injection_order <- 1:nrow(project_run_order)
project_run_order_html <- htmlTable(project_run_order)

htmltools::save_html(project_run_order_html, file = paste(project_dir_html, "/", project_name, "_", user_name, "_run_order_check.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_run_order_check.html", sep="")) #open plotly widget in internet browser

temp_answer <- "blank"
temp_answer_2 <- "blank"
while(temp_answer_2 != "yes"){
while(temp_answer != "yes" & temp_answer != "no"){
  temp_answer <- dlgInput("A worklist has just opened in your browser.  Does this match the run order of your analysis?", "yes/no")$res
}

if(temp_answer == "no"){
  dlg_message("OK. Please upload a worklist template csv file now. It will need 3x columns: sampleID, PlateID and injection_order. A template file has been created in your project directory (run_order_template.csv)", type = 'ok')
  temp_tibble <- project_run_order
  temp_tibble$injection_order <- NA
  write_csv(temp_tibble, 
            file = paste(project_dir, "/", Sys.Date(), "_run_order_template.csv", sep=""))
  dlg_message("Select this file now", type = 'ok')
  new_project_run_order <- file.choose(.) %>% read_csv()
  colnames(new_project_run_order) <- c("sampleID", "plateID", "injection_order")
}

individual_lipid_data$run_order <- NA
for(idx_ro in 1:nrow(new_project_run_order)){
  #browser()
  #add run order value from worklist template to individual_lipid_data 
  individual_lipid_data$run_order[grep(new_project_run_order$sampleID[idx_ro], individual_lipid_data$sampleID)] <- new_project_run_order$injection_order[idx_ro]
  #add plate number order value from worklist template to individual_lipid_data 
  individual_lipid_data$plateID[grep(new_project_run_order$sampleID[idx_ro], individual_lipid_data$sampleID)] <- new_project_run_order$plateID[idx_ro]
}

individual_lipid_data <- individual_lipid_data %>% arrange(run_order)

new_project_run_order <- individual_lipid_data %>% select(sampleID, plateID, run_order)
colnames(new_project_run_order) <- c("sampleID", "plateID", "injection_order")
new_project_run_order <- new_project_run_order %>% filter(!is.na(injection_order))
new_project_run_order_html <- htmlTable(new_project_run_order)

htmltools::save_html(new_project_run_order_html, file = paste(project_dir_html, "/", project_name, "_", user_name, "_run_order_check.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_run_order_check.html", sep="")) #open plotly widget in internet browser

temp_answer_2 <- dlgInput("A new worklist order has just opened in your browser.  Does this match the run order of your analysis?", "yes/no")$res
if(temp_answer_2 == "no"){
  temp_answer <- "no"
}
}

individual_lipid_data <- individual_lipid_data %>% add_column(individual_lipid_data$run_order, .before = 3, .name_repair = "minimal") %>% select(-run_order)
colnames(individual_lipid_data)[3] <- "run_order"


individual_lipid_data <- individual_lipid_data %>% filter(!grepl("conditioning", sampleID))
class_lipid_data <- create_lipid_class_data_summed(individual_lipid_data)

new_project_run_order <- new_project_run_order %>% filter(!grepl("conditioning", sampleID))
plateID <- individual_lipid_data$plateID
run_order <- individual_lipid_data$run_order

