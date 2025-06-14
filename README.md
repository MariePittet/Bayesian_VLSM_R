# Bayesian_VLSM_R
Contributor: Marie Pittet, marie.pittet93@gmail.com, Jan 2025

Full R-based Bayesian Voxel-Based Lesion-Symptom Mapping with publication ready plots of selected slices. This is particularly suitable for the exploration of the effects of lesions on behavioral outcomes (symptoms) in small samples. 
A first analysis using the frequentist approach (t-tests) is performed for comparison. The second approach is the bVLSM using bayesian t-tests (with the "BayesFactor" package), with Cauchy priors. BVLSM -also called Bayesian Lesion-Deficit Inference (BLDI), does not really benefit from lesion-volume correction (Sperber et al., 2023), hence it is not applied here. A visualization in the axial plane is provided, as well as a visualization of selected slices in the sagittal, coronal, and axial planes. Additionnally, the significant voxels are attributed to Talairach atlas' regions or Automated anatomical labeling (AAL) regions in a table for easy interpretation. The AAL atlas had to be resampled to fit the MNI1mm template used for normalization.

## Input: 
- A 4d file containing the normalized lesions of all patients
- A vector containing the behavioral scores of interest of all patients. They must appear in the same order as the patients' lesion file so that the lesions and the behavioral scores are correctly attributed to the right patient. Adapt the symptom_scores variable to reflect your behavioral scores of interest.
- For visualization, selected slices are given as examples but can be modified by modifying the numbers in the selected_slices variable for the axial plane, or in the x_slices, y_slices, and y_slices for the visualization in 3 planes. 

## Output:
- A map of log10 Bayesian Factor values (logBF) instead of BF factors. This is because the latter can range between 0 and very large values with non-linear interpretation and is not easy to handle within visualization softwares. The logBF map can be opened with neuroimaging visualization softwares (such as ITK-SNAP for a free option) on top of the brain template used for lesion normalization (e.g. MNI152).
- a publication-ready map of the logBF map on top of the MNI template in one plane (axial), or three planes (sagittal, axial, coronal).
- a table with the attribution of significant voxels to Talairach atlas regions or AAL regions.

 
 
 ## Log(BF) values are usually interpreted as follows (Kass & Raftery, 1995):
- logBF = 0                no evidence for either the null or the alternative hypothesis
- logBF = [0, 0.5[         weak evidence for alternative hypothesis
- logBF = [0.5, 1[         Substantial evidence for alternative hypothesis
- logBF = [1, 1.5[         Strong evidence for alternative hypothesis
- logBF = [1.5, 2[         Very strong evidence for alternative hypothesis
- logBF = [2, inf]         Decisive evidence for alternative hypothesis

 ## Quick note on bVLSM or VLSM on stroke patients
I feel that VLSM conducted on samples composed exclusively of stroke patients is a tricky subject, even with the bayesian approach. This caution stems from a fundamental anatomical issue tied to stroke pathology: most vascular occlusions occur near the Circle of Willis, aka the central arterial hub supplying the brain with blood from its ventral surface (Sperber & Karnath, 2015).
Because of this, lesions often originate in proximal arterial territories and extend to more distal regions. When the function of interest lies in dorsal regions (e.g., parietal cortex), VLSM may falsely implicate upstream areas, not because they are causally related to the behavior, but because they are simply part of the same vascular tree.
This goes beyond the usual concerns about lesion coverage or lesion dependence: even with sufficient coverage, VLSM can produce misleading statistical associations along the entire arterial path leading to the functionally relevant region.
For example, a sensorimotor variable that depends on parietal cortex function may result in significant VLSM hits in more ventral regions like Heschl's gyrus, simply because those regions are commonly affected in MCA strokes and almost always co-occur with damage to the true site.
Think of these proximal regions as the office manager who gets all the credit while the real work is done upstream (by the parietal cortex in my toy example). The boss didn’t do the job, but still shows up in the VLSM results.



_References:
Sperber, C., et al. (2023). Bayesian lesion-deficit inference with Bayes factor mapping: Key advantages, limitations, and a toolbox. NeuroImage. \
Talairach, J., & Tournoux, P. (1988). Co-Planar Stereotaxic Atlas of the Human Brain. Thieme Medical Publishers.\
Tzourio-Mazoyer, N., et al. (2002). Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of the MNI MRI single-subject brain. Neuroimage.\
Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical Association.\
Sperber, C., Karnath, O. (2015). Topography of acute stroke in a sample of 439 right brain damaged patients. NeuroImage Clin._
