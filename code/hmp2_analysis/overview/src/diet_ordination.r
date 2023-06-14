
source("./common/load_diet.r")

# Filter out NA's
hmp2_diet_num_nona.pcl <- hmp2_diet_num.pcl %>%
    pcl.filter.s(!any(is.na(x)))

# Merge disease activity so we can overlay it on the ordination
source("./common/disease_activity.r")
hmp2_diet_num_nona.pcl <- hmp2_diet_num_nona.pcl %>%
    merge_disease_activity

# Get the PCoA
hmp2_diet_num_nona.ord <- pcl.pcoa(hmp2_diet_num_nona.pcl)

hmp2_diet_num_nona.ord <- pcl.tsne(hmp2_diet_num_nona.pcl)



# Plots

source("./common/disease_colors.r")
pdf("./overview/diet_ordination.pdf", 5, 6)
pcl.ordplot(hmp2_diet_num_nona.pcl, hmp2_diet_num_nona.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors,
            size_abs=2) +
    theme(legend.position="bottom", legend.direction="vertical")

pcl.ordplot(hmp2_diet_num_nona.pcl, hmp2_diet_num_nona.ord,
            colour="active", colour_override=c("TRUE"="red", "FALSE"="blue"),
            size_abs=2, sortby = "active") +
    theme(legend.position="bottom", legend.direction="vertical")

for (what in c(
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
    "Sweets..pies..jam..chocolate..cake..cookies..etc..")) {
    print(pcl.ordplot(hmp2_diet_num_nona.pcl, hmp2_diet_num_nona.ord,
                colour=what, size_abs=2, sortby = "active") +
        theme(legend.position="bottom", legend.direction="vertical"))
}

dev.off()


