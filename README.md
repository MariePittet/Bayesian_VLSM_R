# Bayesian_VLSM_R
Contributor: Marie Pittet, marie.pittet93@gmail.com
Jan 2025

Full R-based Bayesian Voxel-Based Lesion-Symptom Mapping. This is particularly suitable for the exploration of the effects of lesions on behavioral outcomes (symptoms) in small sample sizes. 
A first analysis using the frequentist approach (t-tests) is performed for comparison. The second approach is the bVLSM using bayesian t-tests (with the "BayesFactor" package), with Cauchy priors. BVLSM -also called Bayesian Lesion-Deficit Inference (BLDI), does not really benefit from lesion-volume correction (Sperber et al., 2023), hence it is not applied here.

Input: 
- A 4d file containing the normalized lesions of all patients
- A vector containing the behavioral scores of interest of all patients. They must appear in the same order as the patients' lesion file so that the lesions and the behavioral scores are correctly attributed to the right patient. Adapt the symptom_scores variable to reflect your behavioral scores of interest.

Output:
A map of log10 Bayesian Factor values (logBF) rather thanBF factors, since the latter can range between 0 and very large values with non-linear interpretation. The logBF map can be opened with neuroimaging visualization softwares (such as ITK-SNAP for a free option) on top of the brain template used for lesion normalization (e.g. MNI152). Log(BF) values are usually interpreted as follows:

- logBF = 0                no evidence for either the null or the alternative hypothesis
- logBF = [0, 0.5[         weak evidence for alternative hypothesis
- logBF = [0.5, 1[         Substantial evidence for alternative hypothesis
- logBF = [1, 1.5[         Strong evidence for alternative hypothesis
- logBF = [1.5, 2[         Very strong evidence for alternative hypothesis
- logBF = [2, inf]         Decisive evidence for alternative hypothesis
