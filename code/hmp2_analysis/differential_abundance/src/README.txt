###############################
# DESCRIPTION OF DA FUNCTIONS #
###############################

core_DA_functions
	- All differential abundance (DA) functions (see below).

convert_tables
	- Converts raw output files to a single formatted table.

test_DA
	-Analysis 1: Tests for DA in a full model that includes Diagnosis as fixed effect and Disease Active Status (Binary) as nested fixed effect within Diagnosis. Also adjusts for Immunosuppressants, Antibiotics, sex, consent age, and Bowel Surgery. Subjects are modelled as random effects. Filtering + Normalization/Transformation is applied to each data type separately before running the model.

	Analysis 2: Tests for DA in a full model that includes Diagnosis, Disease Active Score as fixed effects. Also adjusts for Immunosuppressants, Antibiotics, sex, consent age, and Bowel Surgery. Subjects are modelled as random effects. Filtering + Normalization/Transformation is applied to each data type separately before running the model.

	Analysis 3: Repeats Analysis 1 with only baseline samples (hence no random effects rest being the same).
	
	Analysis 4: Tests for DA in a full model that includes PreActive Disease as main effect. NOTE, this DOES NOT include Diagnosis as covariate, and fits separate models to CD and UC patients. Also adjusts for Immunosuppressants, Antibiotics, sex, consent age, and Bowel Surgery. Subjects are modelled as random effects. Filtering + Normalization/Transformation is applied to each data type separately before running the model.


test_DA_CLR
	All the Analyses Above with CLR-transformed data, as applicable.


