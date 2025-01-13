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
#install.packages("tidyr")
library(tidyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("viridis")
library(viridis)
#install.packages("xml2")
library(xml2)
#install.packages("gt")
library(gt)

# Loading behavioral and lesion data --------------------------------------

# Loading behavioral scores of ABI patients (MASC total scores)
df <- read_excel("SC_results.xlsx")
symptom_scores <- df$MASCtot[1:20] # removing control subjects

# Load the lesion masks of all participants
lesion_data <- readNifti("all_lesions_4d.nii.gz")

# VLSM frequentist approach (FDR corrected t-tests) -------------------------------

# Initializing arrays for t-values and p-values
dim_3d <- dim(lesion_data)[1:3]
t_values <- array(NA, dim = dim_3d)
p_values <- array(NA, dim = dim_3d)

# Voxel-wise t-tests (in each voxel, comparing MASC scores for lesioned vs unlesioned patients in the precise voxel)
for (x in 1:dim_3d[1]) {
  for (y in 1:dim_3d[2]) {
    for (z in 1:dim_3d[3]) {
      lesion_voxel <- lesion_data[x, y, z, ]  # Voxel values for all patients
      
      # Checking if voxel has sufficient data
      if (sum(lesion_voxel) > 2) {  # At least 2 patients with lesions otherwise it's not relevant
        test <- lm(symptom_scores ~ lesion_voxel)
        t_values[x, y, z] <- summary(test)$coefficients[2, "t value"]
        p_values[x, y, z] <- summary(test)$coefficients[2, "Pr(>|t|)"]
      }
    }
  }
}

# Flattening p-values for FDR correction
p_values_flat <- as.vector(p_values)
p_values_flat <- p_values_flat[!is.na(p_values_flat)]  # Remove NAs

# Applying FDR correction
adjusted_p <- p.adjust(p_values_flat, method = "fdr")

# Mapping corrected p-values back into 3D space
p_values_corrected <- array(NA, dim = dim_3d)
p_values_corrected[!is.na(p_values)] <- adjusted_p

# Loading the header of the original NIfTI file
template <- readNifti("all_lesions_4d.nii.gz")

# Saving t-map
writeNifti(t_values, file = "t_map.nii.gz", template = template)

# Saving corrected p-value map
writeNifti(p_values_corrected, file = "p_map_corrected.nii.gz", template = template)

#Saving uncorrected p-values
writeNifti(p_values, file = "p_map.nii.gz", template = template)


# VLSM Bayesian approach --------------------------------------------------

# Loading behavioral scores of ABI patients (MASC total scores) in case you didn't perform the frequentist analysis, you lazy thing
df <- read_excel("SC_results.xlsx")
symptom_scores <- df$MASCtot[1:20]  # Removing control subjects (this is specific to my dataset which contains data from control subjects

# Loading the lesions of all patients (if you didn't already)
lesion_data <- readNifti("all_lesions_4d.nii.gz")

# Initializing arrays for log Bayes Factors
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

# Defining the slices of interest (CHANGE HERE FOR YOUR OWN ANALYSIS)
selected_slices <- c(40,50,60,102,112)

# Filtering on these slices
voxel_data_filtered <- voxel_data %>%
  filter(z %in% selected_slices) %>%
  mutate(logBF = ifelse(logBF < 0.5, NA, logBF))%>%   # thresholding so that values below 0.5 (no evidence for H0 or H1) do not appear
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

# Selecting interesting slices (I visualized them in ITK-SNAP; CHANGE HERE FOR YOUR OWN ANALYSIS)
x_slices<- c(60,102,109)
y_slices<- c(176,140,133)
z_slices <- c(50,60,100)

# Preparing slice data for z-plane
voxel_data_z <- voxel_data %>%
  filter(z %in% z_slices) %>%
  mutate(logBF = ifelse(logBF < 0.5, NA, logBF))%>%   # thresholding so that values below 0.5 do not appear since logBF < 0.5 do not provide evidence for H1
  drop_na()

mni_voxel_data_z <- mni_voxel_data %>%
  filter(z %in% z_slices)

voxel_data_z$z <- factor(voxel_data_z$z, levels = z_slices)
mni_voxel_data_z$z <- factor(mni_voxel_data_z$z, levels = z_slices)

# Converting z to a factor (so that facet_wrap displays them in the right order, hopefully)
voxel_data_filtered$z <- factor(voxel_data_filtered$z, levels = selected_slices)
mni_voxel_data_filtered$z <- factor(mni_voxel_data_filtered$z, levels = selected_slices)

# Preparing slice data for x-plane
voxel_data_x <- voxel_data %>%
  filter(x %in% x_slices) %>%
  mutate(logBF = ifelse(logBF < 0.5, NA, logBF))%>%   # thresholding so that values below 0.5 do not appear
  drop_na()

mni_voxel_data_x <- mni_voxel_data %>%
  filter(x %in% x_slices)

voxel_data_x$x <- factor(voxel_data_x$x, levels = x_slices)
mni_voxel_data_x$x <- factor(mni_voxel_data_x$x, levels = x_slices)

# Preparing slice data for y-plane
voxel_data_y <- voxel_data %>%
  filter(y %in% y_slices) %>%
  mutate(logBF = ifelse(logBF < 0.5, NA, logBF))%>%   # thresholding so that values below 0.5 do not appear
  drop_na()

mni_voxel_data_y <- mni_voxel_data %>%
  filter(y %in% y_slices)

voxel_data_y$y <- factor(voxel_data_y$y, levels = y_slices)
mni_voxel_data_y$y <- factor(mni_voxel_data_y$y, levels = y_slices)

# Defining colors for the logBF scale. Here you can change the thing if you don't like my colors :(
color_scale <- scale_fill_gradientn(
  colours = magma(256)[50:256],  # Skipping the first 25 tones of the magma palette because they are ugly as hell
  limits = c(0.5, 3.15),            # Setting limits for the color scale (CHANGE HERE IF YOU HAVE CRAZY HIGH logBFs, YOU LUCKY PERSON)
  name = expression(logBF),      # Setting the legend name
  oob = scales::squish          # Handling out-of-bound values
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

# Comparison with atlases - Talairach -------------------------------------------------

# Loading logBF map
logBF_map <-  readNIfTI("logBF_map.nii.gz", reorient = FALSE)

# Loading the Talairach label file
label_data <- readNIfTI("Talairach-labels-1mm.nii.gz", reorient = FALSE)

# Dimensions of the label file
cat("Label file dimensions:", dim(label_data), "\n")
cat("Label range (region indices):", range(label_data, na.rm = TRUE), "\n")

# Data frame with logBF data and atlas region data
logBF_data <- data.frame(
  x = rep(1:dim(label_data)[1], each = dim(label_data)[2] * dim(label_data)[3]),
  y = rep(rep(1:dim(label_data)[2], each = dim(label_data)[3]), times = dim(label_data)[1]),
  z = rep(1:dim(label_data)[3], times = dim(label_data)[1] * dim(label_data)[2]),
  logBF = as.vector(logBF_map),  
  region = as.vector(label_data)  # Region indices from the Talairach atlas
)

# Checking the first few rows to verify that I didn't mess up
head(logBF_data)

# Summarizing logBF by region
region_summary <- logBF_data %>%
  group_by(region) %>%
  summarise(
    max_logBF = max(logBF, na.rm = TRUE),
    num_voxels = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(max_logBF))

# Threshold LogBF values
threshold_value <- 0.5
significant_regions <- region_summary %>%
  filter(max_logBF >= threshold_value)

# Reading the .txt file with labels
labels_data <- read.table("Talairach-labels.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Assigning proper column names
colnames(labels_data) <- c("Region_ID", "Label")

# List of significant region indices
significant_regions1 <- significant_regions$region 

# Merging significant regions with labels
significant_labels <- labels_data[labels_data$Region_ID %in% significant_regions, ]

# Merging labels into the analysis results
analysis_results_with_labels <- merge(significant_regions, labels_data, by.x = "region", by.y = "Region_ID")

# Analysis results with labels
analysis_results_with_labels

# Reordering the data by max_logBF, and by n_voxels
sorted_results <- analysis_results_with_labels %>%
  arrange(desc(max_logBF), desc(num_voxels))

sorted_results

# Comparison with atlases - AAL  -------------------------------------------------

# Loading the AAL atlas and labels
atlas_data <- readNIfTI("AAL_resampled_1mm.nii.gz", reorient = TRUE) # had to resample the AAL atlas in FSL to fit the MNI slices/dimensions
xml_file <- read_xml("AAL.xml") # xml file with labels

# Extracting the labels and indices
labels <- xml_find_all(xml_file, ".//label")
indices <- xml_find_all(labels, ".//index")
names <- xml_find_all(labels, ".//name")

# Converting the extracted data into a dataframe
label_data <- tibble(
  index = as.numeric(xml_text(indices)),
  label = xml_text(names))

# Viewing the dimensions and range of the atlas
cat("Atlas dimensions:", dim(atlas_data), "\n")
cat("Label data dimensions:", dim(label_data), "\n")
cat("Atlas region range:", range(atlas_data, na.rm = TRUE), "\n")

# Creating a data frame with all the data
logBF_data <- data.frame(
  x = rep(1:dim(atlas_data)[1], each = dim(atlas_data)[2] * dim(atlas_data)[3]),
  y = rep(rep(1:dim(atlas_data)[2], each = dim(atlas_data)[3]), times = dim(atlas_data)[1]),
  z = rep(1:dim(atlas_data)[3], times = dim(atlas_data)[1] * dim(atlas_data)[2]),
  logBF = as.vector(logBF_map),   
  region = as.vector(atlas_data)  # regions from AAL atlas
)

# Merging with region labels
logBF_data <- logBF_data %>%
  left_join(label_data, by = c("region" = "index"))

# Summarizing the logBF data by region
region_summary <- logBF_data %>%
  group_by(region, label) %>%
  summarise(
    max_logBF = max(logBF, na.rm = TRUE),
    num_voxels = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(max_logBF))

# Thresholding LogBF values to focus on significant regions
threshold_value <- 0.5
significant_regions <- region_summary %>%
  filter(max_logBF >= threshold_value)

# Reordering the data by max_logBF, and by n_voxels
sorted_results <- significant_regions %>%
  arrange(desc(max_logBF), desc(num_voxels))

sorted_results

# Removing rows where region_name is NA
filtered_results <- sorted_results %>%
  filter(!is.na(label))

filtered_results

# Creating a cute publication-ready table
  gt(data = filtered_results) %>%
  tab_header(
    title = "Regions showing the highest log(BF) values ",
  ) %>%
  cols_label(
    region = "Region ID",
    max_logBF = "Max logBF",
    num_voxels = "Voxel Count",
    label = "Region Name"
  ) %>%
 tab_options(
      table.font.size = px(12),    # Reduce font size because otherwise the thing gets huge. You can adjust
      data_row.padding = px(2)    # Adjust row padding to make it more or less compact
  ) %>%
  fmt_number(
    columns = c(max_logBF),
    decimals = 2
  ) %>%
  fmt_number(
    columns = num_voxels,
    decimals = 0
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  opt_table_outline() %>%
  opt_align_table_header(align = "center")
