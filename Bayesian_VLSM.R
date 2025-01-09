#install.packages("readxl")
library("readxl")
#install.packages("RNifti")
library(RNifti)
#install.packages("neurobase")
library(neurobase)
#install.packages("stats")
library(stats)
#install.packages("oronifti")
library(oro.nifti)
#install.packages(BayesFactor)
library(BayesFactor)


# Loading behavioral and lesion data --------------------------------------

# loading behavioral scores of ABI patients (MASC total scores)
df <- read_excel("SC_results.xlsx")
symptom_scores <- df$MASCtot[1:20] # removing control subjects

# Load the lesion mask (e.g., from a .nii file)
lesion_data <- readNifti("all_lesions_4d.nii.gz")

# VLSM frequentist approach (FDR corrected t-tests) -------------------------------

# Initialize arrays for t-values and p-values
dim_3d <- dim(lesion_data)[1:3]
t_values <- array(NA, dim = dim_3d)
p_values <- array(NA, dim = dim_3d)

# Voxel-wise t-tests
for (x in 1:dim_3d[1]) {
  for (y in 1:dim_3d[2]) {
    for (z in 1:dim_3d[3]) {
      lesion_voxel <- lesion_data[x, y, z, ]  # Voxel values for all patients
      
      # Check if voxel has sufficient data
      if (sum(lesion_voxel) > 2) {  # At least 2 patients with lesions
        test <- lm(symptom_scores ~ lesion_voxel)
        t_values[x, y, z] <- summary(test)$coefficients[2, "t value"]
        p_values[x, y, z] <- summary(test)$coefficients[2, "Pr(>|t|)"]
      }
    }
  }
}

# Flatten p-values for FDR correction
p_values_flat <- as.vector(p_values)
p_values_flat <- p_values_flat[!is.na(p_values_flat)]  # Remove NAs

# Apply FDR correction
adjusted_p <- p.adjust(p_values_flat, method = "fdr")

# Map corrected p-values back into 3D space
p_values_corrected <- array(NA, dim = dim_3d)
p_values_corrected[!is.na(p_values)] <- adjusted_p

# Load the header of the original NIfTI file
template <- readNifti("all_lesions_4d.nii.gz")

# Save t-map
writeNifti(t_values, file = "t_map.nii.gz", template = template)

# Save corrected p-value map
writeNifti(p_values_corrected, file = "p_map_corrected.nii.gz", template = template)

#Save uncorrected p-values
writeNifti(p_values, file = "p_map.nii.gz", template = template)


# VLSM Bayesian approach --------------------------------------------------

# Loading behavioral scores of ABI patients (MASC total scores)
df <- read_excel("SC_results.xlsx")
symptom_scores <- df$MASCtot[1:20]  # Removing control subjects

# Load the lesions of all patients
lesion_data <- readNifti("all_lesions_4d.nii.gz")

# Initialize arrays for Bayes Factors
dim_3d <- dim(lesion_data)[1:3]
logBF <- array(NA, dim = dim_3d)

# Voxel-wise Bayesian t-tests
for (x in 1:dim_3d[1]) {
  for (y in 1:dim_3d[2]) {
    for (z in 1:dim_3d[3]) {
      lesion_voxel <- lesion_data[x, y, z, ]  # Voxel values for all patients
      
      # Check if voxel has sufficient data (at least 2 in each group)
      if (sum(lesion_voxel == 1) >= 2 && sum(lesion_voxel == 0) >= 2) {
        # Split symptom scores based on lesion presence/absence
        group_present <- symptom_scores[lesion_voxel == 1]
        group_absent <- symptom_scores[lesion_voxel == 0]
        
        # Bayesian t-test for group comparison
        bf_result <- ttestBF(x = group_present, y = group_absent)
        bf<- as.numeric(extractBF(bf_result)['bf']) # forceful extraction of bf from weirdly set-up S4 class ^^'
        
        # log Bayes Factor
        logBF[x, y, z] <- as.numeric(log(bf))
      }
    }
  }
}

# Save log Bayes Factor map
template<- lesion_data
writeNifti(logBF, file = "logBF_map.nii.gz", template = template)