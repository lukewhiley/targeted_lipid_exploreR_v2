# produces a PCA of the normalized samples

#label data
rsd_filtered_data$sample_class <- "sample"
rsd_filtered_data$sample_class[grep(paste0(qc_type), rsd_filtered_data$sampleID)] <- paste0(qc_type)

rsd_filtered_class_data$sample_class <- "sample"
rsd_filtered_class_data$sample_class[grep(paste0(qc_type), rsd_filtered_class_data$sampleID)] <- paste0(qc_type)

#run function
pca_check_status <- "change"
while(pca_check_status == "change"){

pca_p <- lipids_pca_ltr(FUNC_individual_multivariate_data = rsd_filtered_data, 
                        FUNC_family_multivariate_data = rsd_filtered_class_data, 
                        FUNC_multivariate_class = "sample_class", 
                        FUNC_plot_label = "sampleID", 
                        FUNC_scaling = "option", 
                        FUNC_qc_label = qc_type)

saveWidget(pca_p[[1]][[1]], file = paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_all_lipids.html", sep=""))# save plotly widget
#browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_all_lipids.html", sep="")) #open plotly widget in internet browser
saveWidget(pca_p[[2]][[1]], file = paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_lipid_class.html", sep=""))# save plotly widget
#browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_lipid_class.html", sep="")) #open plotly widget in internet browser


# pca_check_status <- dlgInput("Check the PCA plots. Are you happy to continue? or do wish to change the scalling type?", "continue/change")$res
# while(pca_check_status != "continue" & pca_check_status != "change"){
#   pca_check_status <- dlgInput("Check the PCA plots. Are you happy to continue? or do wish to change the scalling type?", "continue/change")$res
# }
pca_check_status <- "continue"
}

scale_used <- pca_p[[3]][1]


