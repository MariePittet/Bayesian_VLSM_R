# Setting up the crazy big environment --------------------------------------
#install.packages("readxl")
library(readxl)
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
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("patchwork")
library(patchwork) 
#install.packages("ggnewscale")
library(ggnewscale)
#install.packages("dplyr")
library(dplyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("viridis")
library(viridis)

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

# Visualization in axial view with manual extraction of interesting slices ------------------
logBF_map <- readNIfTI("logBF_map.nii.gz", reorient = FALSE)
mni_template <- readNIfTI("MNI152_T1_1mm.nii.gz", reorient = FALSE)

# Converting to arrays
logBF_array <- as.array(logBF_map)
mni_array   <- as.array(mni_template)

# Creating a data frame with x, y, z for the logBF map
dims <- dim(logBF_array)
voxel_data <- expand.grid(x = 1:dims[1], y = 1:dims[2], z = 1:dims[3])
voxel_data$logBF <- as.vector(logBF_array)

# Creating a data frame for the MNI template
mni_voxel_data <- expand.grid(x = 1:dims[1], y = 1:dims[2], z = 1:dims[3])
mni_voxel_data$intensity <- as.vector(mni_array)

# Defining the slices of interest
selected_slices <- c(40,50,60,102,112)

# Filtering on these slices
voxel_data_filtered <- voxel_data %>%
  filter(z %in% selected_slices) %>%
  drop_na()

mni_voxel_data_filtered <- mni_voxel_data %>%
  filter(z %in% selected_slices)

# Converting z to a factor (so that facet_wrap displays them in the right order)
voxel_data_filtered$z <- factor(voxel_data_filtered$z, levels = selected_slices)
mni_voxel_data_filtered$z <- factor(mni_voxel_data_filtered$z, levels = selected_slices)

# Initializing the plot
p_manual <- ggplot() +
  
  # Plotting the MNI background in grayscale
  geom_raster(data = mni_voxel_data_filtered, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient(name = "MNI Intensity", low = "black", high = "white", guide = "none") +
  
  # Starting a new fill scale for the logBF map so that ggplot doesnt take the MNI intensity
  new_scale_fill() +
  
  # Overlaying the logBF map
  geom_raster(data = voxel_data_filtered, aes(x = x, y = y, fill = logBF), alpha = 0.7) +
  scale_fill_viridis_c(option = "magma", name = expression(logBF), limits = c(0, NA),oob = scales::squish)+
  
  # Facetting by z
  facet_wrap(~z, ncol = 5) +
  
  # Working on the theme
  theme_void(base_size = 14) +
  coord_fixed() +
  theme(legend.position = "right", panel.spacing = unit(0.3, "lines")) +
  labs(title = "Bayesian VLSM")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_manual

# Visualization in 3 planes with manual slice extraction ------------------

#Selecting the interesting slices (I visualized the in ITK-SNAP)
x_slices<- c(60,102,109)
y_slices<- c(176,140,133)
z_slices <- c(50,60,100)

# Preparing slice data for z-plane
voxel_data_z <- voxel_data %>%
  filter(z %in% z_slices) %>%
  drop_na()

mni_voxel_data_z <- mni_voxel_data %>%
  filter(z %in% z_slices)

voxel_data_z$z <- factor(voxel_data_z$z, levels = z_slices)
mni_voxel_data_z$z <- factor(mni_voxel_data_z$z, levels = z_slices)

# Converting z to a factor (so that facet_wrap displays them in the right order)
voxel_data_filtered$z <- factor(voxel_data_filtered$z, levels = selected_slices)
mni_voxel_data_filtered$z <- factor(mni_voxel_data_filtered$z, levels = selected_slices)

# Preparing slice data for x-plane
voxel_data_x <- voxel_data %>%
  filter(x %in% x_slices) %>%
  drop_na()

mni_voxel_data_x <- mni_voxel_data %>%
  filter(x %in% x_slices)

voxel_data_x$x <- factor(voxel_data_x$x, levels = x_slices)
mni_voxel_data_x$x <- factor(mni_voxel_data_x$x, levels = x_slices)

# Preparing slice data for y-plane
voxel_data_y <- voxel_data %>%
  filter(y %in% y_slices) %>%
  drop_na()

mni_voxel_data_y <- mni_voxel_data %>%
  filter(y %in% y_slices)

voxel_data_y$y <- factor(voxel_data_y$y, levels = y_slices)
mni_voxel_data_y$y <- factor(mni_voxel_data_y$y, levels = y_slices)

# Defining colors for the logBF scale. Here you can change the thing if you don't like my colors :(
color_scale <- scale_fill_gradientn(
  colours = magma(256)[50:256],  # Skip the first 25 tones of the magma palette
  limits = c(0.5, 3.15),            # Set limits for the color scale
  name = expression(logBF),      # Set the legend name
  oob = scales::squish          # Handle out-of-bound values
)

# Plotting for z-plane
plot_z <- ggplot() +
  geom_raster(data = mni_voxel_data_z, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient(name = "MNI Intensity", low = "black", high = "white", guide = "none") +
  new_scale_fill() +
  geom_tile(data = voxel_data_z, aes(x = x, y = y, fill = logBF), alpha = 0.7) +
  color_scale +
  facet_wrap(~z, ncol = 3, strip.position = "bottom") +
  theme_void() +
  coord_fixed() +
  labs(title = "Axial view")+
  theme(legend.position = "right", panel.spacing = unit(0.5, "lines")) +
  theme(plot.title = element_text(size = 14, face = "italic", hjust = 0.5))

# Plotting for y-plane
plot_y <- ggplot() +
  geom_tile(data = mni_voxel_data_y, aes(x = x, y = z, fill = intensity)) +
  scale_fill_gradient(name = "MNI Intensity", low = "black", high = "white", guide = "none") +
  new_scale_fill() +
  geom_tile(data = voxel_data_y, aes(x = x, y = z, fill = logBF), alpha = 0.7) +
  color_scale +
  facet_wrap(~y, ncol = 3, strip.position = "bottom") +
 theme_void() +
  coord_fixed() +
  labs(title = "Coronal view") +
  theme(legend.position = "right", panel.spacing = unit(0.3, "lines")) +
  theme(plot.title = element_text(size = 14, face = "italic", hjust = 0.5))
  
# Plotting for x-plane
plot_x <- ggplot() +
  geom_tile(data = mni_voxel_data_x, aes(x = y, y = z, fill = intensity)) +
  scale_fill_gradient(name = "MNI Intensity", low = "black", high = "white", guide = "none") +
  new_scale_fill() +
  geom_tile(data = voxel_data_x, aes(x = y, y = z, fill = logBF), alpha = 0.7) +
  color_scale +
  facet_wrap(~x, ncol = 3,strip.position = "bottom") +
  theme_void() +
  coord_fixed() +
  labs(title = "Sagittal view") +
  theme(legend.position = "right", panel.spacing = unit(0.1, "lines")) +
  theme(plot.title = element_text(size = 14, face = "italic", hjust = 0.5))

# Combining the plots into one with a single legend 
combined_plot <- (plot_x | plot_y | plot_z) +
  plot_layout(nrow = 3) +
  plot_annotation(title = "") &
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    strip.text = element_text(hjust = 0.5, size = 12),
    legend.position = "right",  # Position the legend to the right
    legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10))

combined_plot 
