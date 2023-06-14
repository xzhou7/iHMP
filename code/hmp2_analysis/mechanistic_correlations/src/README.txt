####################################################
# DESCRIPTION OF MECHANISTIC CORRELATION FUNCTIONS #
####################################################

core_mechanistic_functions

	- All utility functions for mechanistic correlations (see below).

test_mechanistic_correlations

	- Extracts two types of residuals for each data type: (i) from a full model that includes Diagnosis as fixed effect and Disease Active Status (Binary) as nested fixed effect within Diagnosis while adjusting for Immunosuppressants, Antibiotics, sex, consent age, and Bowel Surgery, and (ii) from an empty model that does not include any covariates. Both of these models consider subjects as random effects. Filtering + Normalization/Transformation is applied to each data type separately before running the model. After residuals are generated, they are matched pairwise to have same samples across paired multi’omics and finally appropriately formatted output are generated for downstream HAllA analysis.

test_mechanistic_correlations_baseline

	- Extracts two types of residuals for each data type, after subsetting the data to the baseline samples: (i) from a full model that includes Diagnosis as fixed effect and Disease Active Status (Binary) as nested fixed effect within Diagnosis while adjusting for Immunosuppressants, Antibiotics, sex, consent age, and Bowel Surgery, and (ii) from an empty model that does not include any covariates. Filtering + Normalization/Transformation is applied to each data type separately before running the model. After residuals are generated, they are matched pairwise to have same samples across paired multi’omics and finally appropriately formatted output are generated for downstream HAllA analysis.