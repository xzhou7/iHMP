#TODO: set HAllA input and output path

halla_input_path=~/HMP2/analysis/mechanistic_correlations/halla_input
halla_output_path=~/HMP2/analysis/mechanistic_correlations/halla_output  #/Dropbox\ \(Huttenhower\ Lab\)

# load HAllA on server
hutlab load centos6/halla/latest

declare -a types=("empty" "full")
declare -a set1=("bugs" "ecDNAs" "ecRNAs" "metabolites" "proteinECs" "bugkoRNADNAs")
declare -a set2=("serology" "fecalcal" "HTXrectum" "HTXileum")
#"bugs" "ecDNAs" "ecRNAs" "metabolites" "proteinECs" 

for residuals in "${types[@]}"
do
	echo "$residuals"
	for i in "${set1[@]}"
	do
		echo "$i"
		for j in "${set2[@]}"
		do
			
			if [ -f ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv ];
	        then
	        	   echo "$j"
	        	   ls ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv
	               halla -X ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv \
				   -Y ${halla_input_path}/${i}-${j}_${j}_residuals_${residuals}.tsv \
				   -m spearman --header -o ${halla_output_path}/${i}-${j}_residuals_${residuals}_q05 \
				   -q .05 --missing-char "NA"  --format-feature-names &
				   halla -X ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv \
				   -Y ${halla_input_path}/${i}-${j}_${j}_residuals_${residuals}.tsv \
				   -m spearman --header -o ${halla_output_path}/${i}-${j}_residuals_${residuals}_q05_alla -a AllA \
				   -q .05 --missing-char "NA"  --format-feature-names &
	        fi
	        if [ -f "${halla_input_path}"/${j}-${i}_${j}_residuals_${residuals}.tsv ];
	        then
	              echo "$j"
	              halla -X ${halla_input_path}/${j}-${i}_${j}_residuals_${residuals}.tsv \
				  -Y ${halla_input_path}/${j}-${i}_${i}_residuals_${residuals}.tsv \
				  -m spearman --header -o ${halla_output_path}/${j}-${i}_residuals_${residuals}_q05 \
				  -q .05 --missing-char "NA"  --format-feature-names &
				  halla -X ${halla_input_path}/${j}-${i}_${j}_residuals_${residuals}.tsv \
				  -Y ${halla_input_path}/${j}-${i}_${i}_residuals_${residuals}.tsv \
				  -m spearman --header -o ${halla_output_path}/${j}-${i}_residuals_${residuals}_q05_alla -a AllA \
				  -q .05 --missing-char "NA"  --format-feature-names &
	        fi
	        sleep 2
		done
	done
done
# this part runs HAllA and AllA  for each data set against itself 
# to get clusters within each data set and statistics info for pairs respectively 
declare -a set1=("bugs" "ecDNAs" "ecRNAs" "metabolites" "proteinECs" "bugkoRNADNAs" "serology" "fecalcal" "HTX-rectum" "HTX-ileum")
declare -a set2=("bugs" "ecDNAs" "ecRNAs" "metabolites" "proteinECs" "bugkoRNADNAs" "serology" "fecalcal" "HTX-rectum" "HTX-ileum")
for residuals in "${types[@]}"
do
	echo "$residuals"
	for i in "${set1[@]}"
	do
		echo "$i"
		for j in "${set2[@]}"
		do
			
			if [ -f ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv ];
	        then
	        	   echo "$j"
	        	   ls ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv
	               halla -X ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv \
				  -Y ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv \
				  -m spearman --header -o ${halla_output_path}/${i}-${i}_residuals_${residuals}_q05 \
				  -q .05 --missing-char "NA"  --format-feature-names &
				  halla -X ${halla_input_path}/${i}-${j}_${i}_residuals_${residuals}.tsv \
				  -Y ${halla_input_path}/${i}-${j}_${j}_residuals_${residuals}.tsv \
				  -m spearman --header -o ${halla_output_path}/${i}-${i}_residuals_${residuals}_q05_alla -a AllA \
				  -q .05 --missing-char "NA"  --format-feature-names &
				  break
	        fi
	        if [ -f ${halla_input_path}/${j}-${i}_${j}_residuals_${residuals}.tsv ];
	        then
	              echo "$j"
	              halla -X ${halla_input_path}/${j}-${i}_${j}_residuals_${residuals}.tsv \
				  -Y ${halla_input_path}/${j}-${i}_${j}_residuals_${residuals}.tsv \
				  -m spearman --header -o ${halla_output_path}/${j}-${j}_residuals_${residuals}_q05 \
				  -q .05 --missing-char "NA"  --format-feature-names &
				  halla -X ${halla_input_path}/${j}-${i}_${j}_residuals_${residuals}.tsv \
				  -Y ${halla_input_path}/${j}-${i}_${i}_residuals_${residuals}.tsv \
				  -m spearman --header -o ${halla_output_path}/${j}-${j}_residuals_${residuals}_q05_alla -a AllA \
				  -q .05 --missing-char "NA"  --format-feature-names &
				  break
	        fi
	        sleep 2
		done
	done
done



#  proteinECs-metabolites
## HAllA
#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &
## AllA
#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &


#  proteinECs-metabolites
## HAllA
#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &
## AllA
#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &

#  proteinECs-bugs
## HAllA
#halla -X ${halla_input_path}/proteinECs-bugs_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-bugs_bugs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-bugs_bugs_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &
## AllA
#halla -X ${halla_input_path}/proteinECs-bugs_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-bugs_bugs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-bugs_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &


#  ecDNAs-bugs
## HAllA
#halla -X ${halla_input_path}/ecDNAs-bugs_ecDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecDNAs-bugs_bugs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecDNAs-bugs_bugs_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &
## AllA
#halla -X ${halla_input_path}/ecDNAs-bugs_ecDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecDNAs-bugs_bugs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecDNAs-bugs_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &


#  proteinECs-metabolites_proteinECs
## HAllA
#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &
## AllA
#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs-metabolites_metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &


# bugs-metabolites
## HAllA
#halla -X ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/bugs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/bugs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &
## AllA
#halla -X ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/bugs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/bugs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &


# bugs-koRNAs
## HAllA
#halla -X ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/bugs-koRNAs_residuals_${residuals_type}_q25 -q .25 --missing-char "NA"  --format-feature-names &

## AllA
#halla -X ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/bugs-koRNAs_residuals_${residuals_type}_q25_alla -q .25 --missing-char "NA"  --format-feature-names -a AllA &

# ecDNAs-metabolites
## HAllA
#halla -X ${halla_input_path}/ecDNAs-metabolites_ecDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecDNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecDNAs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &

## AllA
#halla -X ${halla_input_path}/ecDNAs-metabolites_ecDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecDNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecDNAs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &

# ecRNAs-metabolites
## HAllA
#halla -X ${halla_input_path}/ecRNAs-metabolites_ecRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecRNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecRNAs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &

##AllA
#halla -X ${halla_input_path}/ecRNAs-metabolites_ecRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecRNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecRNAs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &

# koDNAs-metabolites
## HAllA
#halla -X ${halla_input_path}/koDNAs-metabolites_koDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koDNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/koDNAs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &

## AllA
#halla -X ${halla_input_path}/koDNAs-metabolites_koDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koDNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/koDNAs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &

# koRNAs-metabolites
## HAllA
#halla -X ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koRNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/koRNAs-metabolites_residuals_${residuals_type}_q05 -q .05 --missing-char "NA"  --format-feature-names &

## AllA
#halla -X ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koRNAs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/koRNAs-metabolites_residuals_${residuals_type}_q05_alla -q .05 --missing-char "NA"  --format-feature-names -a AllA &

# Get stats for within cluster by pairwise testing within datasets (e.g. bugs-bugs)

#halla -X ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/bugs_bugs_alla_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a AllA &

#halla -X ${halla_input_path}/bugs-metabolites_metabolites_residuals_${residuals_type}.tsv -Y ${halla_input_path}/bugs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/metabolites_metabolites_alla_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a AllA &

#halla -X ${halla_input_path}/ecRNAs-metabolites_ecRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecRNAs-metabolites_ecRNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecRNAs_ecRNAs_alla_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a AllA &

#halla -X ${halla_input_path}/ecDNAs-metabolites_ecDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecDNAs-metabolites_ecDNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecDNAs_ecDNAs_alla_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a AllA &

#halla -X ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/koRNAs_koRNAs_alla_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a AllA &

#halla -X ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/proteinECs-metabolites_proteinECs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/proteinECs_proteinECs_alla_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a AllA &

# Get within dataset clusters by blockwise testing within datasets (e.g. bugs-bugs)

#halla -X ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/bugs-metabolites_bugs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/bugs_bugs_clust_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a HAllA &

#halla -X ${halla_input_path}/bugs-metabolites_metabolites_residuals_${residuals_type}.tsv -Y ${halla_input_path}/bugs-metabolites_metabolites_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/metabolites_metabolites_clust_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a HAllA &

#halla -X ${halla_input_path}/ecRNAs-metabolites_ecRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecRNAs-metabolites_ecRNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecRNAs-ecRNAa_clust_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a HAllA &

#halla -X ${halla_input_path}/ecDNAs-metabolites_ecDNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/ecDNAs-metabolites_ecDNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/ecDNAs_ecDNAs_clust_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a HAllA &

#halla -X ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -Y ${halla_input_path}/koRNAs-metabolites_koRNAs_residuals_${residuals_type}.tsv -m spearman --header -o ${halla_output_path}/koRNAs_koRNAs_clust_${residuals_type} -q 1 --missing-char "NA"  --format-feature-names -a HAllA &


