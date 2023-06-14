
source("./common/pcl_utils.r")
source("./common/merge_metadata.r")
source("./common/fix_metadata.r")

# Pull dietary data for any sample with a microbial data type
hmp2_sample_metadata_ssc <- hmp2_sample_metadata[
    hmp2_sample_metadata$data_type == "metagenomics"
    | hmp2_sample_metadata$data_type == "metatranscriptomics"
    | hmp2_sample_metadata$data_type == "proteomics"
    | hmp2_sample_metadata$data_type == "metabolomics"
    | hmp2_sample_metadata$data_type == "viromics",]
hmp2_sample_metadata_ssc <- hmp2_sample_metadata_ssc[!duplicated(hmp2_sample_metadata_ssc$site_sub_coll),]
rownames(hmp2_sample_metadata_ssc) <- as.character(hmp2_sample_metadata_ssc$site_sub_coll)

hmp2_diet <- hmp2_sample_metadata_ssc[,c(
    "Soft.drinks..tea.or.coffee.with.sugar..corn.syrup..maple.syrup..cane.sugar..etc.",
    "Diet.soft.drinks..tea.or.coffee.with.sugar..Stevia..Equal..Splenda.etc.",
    "Fruit.juice..orange..apple..cranberry..prune.etc..",
    "Water",
    "Alcohol..beer..brandy..spirits..hard.liquor..wine..aperitif..etc.." ,
    "Yogurt.or.other.foods.containing.active.bacterial.cultures..kefir..sauerkraut.",
    "Dairy..milk..cream..ice.cream..cheese..cream.cheese.",
    "Probiotic",
    "Fruits..no.juice...Apples..raisins..bananas..oranges..strawberries..blueberries",
    "Vegetables..salad..tomatoes..onions..greens..carrots..peppers..green.beans..etc.",
    "Beans..tofu..soy..soy.burgers..lentils..Mexican.beans..lima.beans.etc.",
    "Whole.grains..wheat..oats..brown.rice..rye..quinoa..wheat.bread..wheat.pasta.",
    "Starch..white.rice..bread..pizza..potatoes..yams..cereals..pancakes..etc..",
    "Eggs",
    "Processed.meat..other.red.or.white.meat.such.as.lunch.meat..ham..salami..bologna",
    "Red.meat..beef..hamburger..pork..lamb.",
    "White.meat..chicken..turkey..etc..",
    "Shellfish..shrimp..lobster..scallops..etc..",
    "Fish..fish.nuggets..breaded.fish..fish.cakes..salmon..tuna..etc..",
    "Sweets..pies..jam..chocolate..cake..cookies..etc..")]

# Set appropriate factors
library(plyr)
library(dplyr)
hmp2_diet <- colwise(function(x)factor(x, levels=c(
    "No, I did not consume these products in the last 7 days",
    "Within the past 4 to 7 days",
    "Within the past 2 to 3 days",
    "Yesterday, 1 to 2 times",
    "Yesterday, 3 or more times"
)))(hmp2_diet)
rownames(hmp2_diet) <- rownames(hmp2_sample_metadata_ssc)

# Create a normalized numeric version
hmp2_diet_num <- (colwise(as.numeric)(hmp2_diet) - 1) / 4
rownames(hmp2_diet_num) <- rownames(hmp2_sample_metadata_ssc)

# Package as PCLs to allow merges etc..
hmp2_sample_metadata_ssc$pilot <- F
hmp2_sample_metadata_ssc$subject <- hmp2_sample_metadata_ssc$Participant.ID
hmp2_diet.pcl <- pcl.make(hmp2_diet, meta=hmp2_sample_metadata_ssc) %>% fix_metadata
hmp2_diet_num.pcl <- pcl.make(hmp2_diet_num, meta=hmp2_sample_metadata_ssc) %>% fix_metadata

# Clean up
rm(hmp2_sample_metadata_ssc, hmp2_diet, hmp2_diet_num)

