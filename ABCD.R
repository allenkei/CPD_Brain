#install.packages("remotes")
#remotes::install_github("jonclayden/RNifti")
library(RNifti)


##### MID
#df1 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/MID/sub-NDARINVB9V0KP7H_ses-baselineYear1Arm1_task-MID_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
#df2 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/MID/sub-NDARINVDX7DPE6L_ses-baselineYear1Arm1_task-MID_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
#df3 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/MID/sub-NDARINVEWDP96RH_ses-baselineYear1Arm1_task-MID_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
#df4 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/MID/sub-NDARINVKF3MM2R1_ses-baselineYear1Arm1_task-MID_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")



##### SST
#df1 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/SST/sub-NDARINVB9V0KP7H_ses-baselineYear1Arm1_task-SST_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
#df2 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/SST/sub-NDARINVDX7DPE6L_ses-baselineYear1Arm1_task-SST_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
#df3 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/SST/sub-NDARINVEWDP96RH_ses-baselineYear1Arm1_task-SST_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
#df4 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/SST/sub-NDARINVKF3MM2R1_ses-baselineYear1Arm1_task-SST_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")


##### nBack
df1 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/nBack/sub-NDARINVB9V0KP7H_ses-baselineYear1Arm1_task-nback_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
df2 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/nBack/sub-NDARINVDX7DPE6L_ses-baselineYear1Arm1_task-nback_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
df3 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/nBack/sub-NDARINVEWDP96RH_ses-baselineYear1Arm1_task-nback_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
df4 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/nBack/sub-NDARINVKF3MM2R1_ses-baselineYear1Arm1_task-nback_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")



##### rest
df1 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/rest/sub-NDARINVB9V0KP7H_ses-baselineYear1Arm1_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
df2 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/rest/sub-NDARINVDX7DPE6L_ses-baselineYear1Arm1_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
df3 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/rest/sub-NDARINVEWDP96RH_ses-baselineYear1Arm1_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")
df4 <- readNifti("/Users/allen/Desktop/Brain Networks/abcd_sample/rest/sub-NDARINVKF3MM2R1_ses-baselineYear1Arm1_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii")




df1 <- df1[,,,1,,] # sub 1
df2 <- df2[,,,1,,] # sub 2
df3 <- df3[,,,1,,] # sub 3
df4 <- df4[,,,1,,] # sub 4

dim(df1)
dim(df2)
dim(df3)
dim(df4)





data_list <- list(df1, df2, df3, df4); #rm(df1, df2, df3, df4)



# Parameters
window_size_sec <- 50  # seconds
stride_sec <- 5        # seconds
TR_sec <- 1            # Assume TR = 1 second (change if needed)


n_timepoints <- min(nrow(df1), nrow(df2), nrow(df3), nrow(df4))
window_size <- window_size_sec / TR_sec  # window size in TRs
stride <- stride_sec / TR_sec            # stride in TRs
n_layers <- length(data_list)
n_regions <- ncol(df1)
n_windows <- floor((n_timepoints - window_size) / stride) + 1


network_array <- array(NA, dim = c(n_windows, n_regions, n_regions, n_layers))


for (l in 1:n_layers) {
  fmri_data <- data_list[[l]]
  
  for (i in 1:n_windows) {
    start_idx <- (i - 1) * stride + 1
    end_idx <- start_idx + window_size - 1
    
    window_data <- fmri_data[start_idx:end_idx, ]
    corr_matrix <- cor(window_data)
    binary_matrix <- ifelse(corr_matrix > 0.5, 1, 0)
    network_array[i, , , l] <- binary_matrix
  }
}


dim(network_array)








library(rTensor)
source("SBS.R")
source("CUSUM.R") # source("utility.R")
#source("competitor.R")


dim(network_array)

A.tensor <- network_array[1:136,,,] # check


num_layer <- 4
num_node <- ncol(data_list[[1]])
num_time <- dim(A.tensor)[1]


# proposed method
hat.rank <- c(10, 10, num_layer)
threshold_ratio <- c(2, 5, 10)
threshold_list <- rev(threshold_ratio * num_node*sqrt(num_layer)*(log(num_time/2))^(3/2))
intervals <- construct_intervals(num_time/2, sqrt(1/2), 4)




A.tensor.even <- A.tensor[seq(2, num_time, by = 2), , , ]
B.tensor.odd  <- A.tensor[seq(1, num_time-1, by = 2), , , ] 

gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = TRUE, intervals, obj.B = B.tensor.odd)
results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_time/2, CUSUM_res = gains, verbose = FALSE,
                               threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)

output <- list()

for (i in 1:(length(results_g)-1)) {
  detected_CP_g <- sort(results_g[[i+1]]$results[, 1])
  detected_CP_gl1 <- refinement1(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
  
  output[[i]] <- list()
  output[[i]]$threshold <- results_g[[i+1]]$threshold
  output[[i]]$thres_ratio <- rev(threshold_ratio)[i]
  output[[i]]$detected_CP <- 2* detected_CP_gl1
  
}

output
#save(output, file = paste0("real_data/results/ABCD_proposed_nBack.RData"))


