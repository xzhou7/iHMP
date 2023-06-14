
include("cytoscape_functions.jl")

# Path to HMP2 dropbox folder - change as necessary
const HMP2_dropbox = "/Users/kev/Dropbox (Huttenhower Lab)/HMP2/"

const halla_dir = joinpath(HMP2_dropbox, "analysis/mechanistic_correlations/halla_output/")
const corrdir = joinpath(HMP2_dropbox, "analysis/differential_abundance/Vanilla/active_vs_inactive/")
const htxcorrdir = joinpath(HMP2_dropbox, "analysis/host_transcriptomics")

# Model to use ("full" or "empty")
const model = "empty"

#const cytoscape_output = joinpath(HMP2_dropbox, "analysis/cytoscape_network/")
const cytoscape_output = joinpath(HMP2_dropbox, "analysis/cytoscape_network/$(model)/")


# Determine and dump node activity

minq = 0.05

active_bugs = nodesactive(joinpath(corrdir, "Active_IBD_relations_bugs.tsv"), "bugs"; minq=minq, prep="DNA| ")
active_bugkoRNADNAs = nodesactive(joinpath(corrdir, "Active_IBD_relations_bugkoRNADNAs.tsv"), "bugkoRNADNAs"; minq=0.1, prep="RNA| ")
active_ecDNAs = nodesactive(joinpath(corrdir, "Active_IBD_relations_ecDNAs.tsv"), "ecDNAs"; minq=minq, prep="DNA| ")
active_ecRNAs = nodesactive(joinpath(corrdir, "Active_IBD_relations_ecRNAs.tsv"), "ecRNAs"; minq=minq, prep="RNA| ")
active_proteinECs = nodesactive(joinpath(corrdir, "Active_IBD_relations_proteinECs.tsv"), "proteinECs"; minq=minq, prep="Prot| ")
active_metabolites = nodesactive(joinpath(corrdir, "Active_IBD_relations_metabolites.tsv"), "metabolites"; minq=minq)
active_fecalcal = makenodeactive("fecalcal", "Calprotectin", "fecalcal", "UC")
active_HTXileum = nodesactive(joinpath(htxcorrdir, "DEGs_ileum.txt"), "HTXileum"; minq=minq, prep="Ileum| ")
active_HTXrectum = nodesactive(joinpath(htxcorrdir, "DEGs_rectum.txt"), "HTXrectum"; minq=minq, prep="Rectum| ")
active_serology = nodesactive(joinpath(corrdir, "#Serology_Table.txt"), "Serology"; minq=minq, prep="Sera| ")


all_edges = DataFrame()
type_prefixes = Dict(
	"bugs"        => "DNA| ",
	"bugkoRNADNAs"=> "RNA| ",
	"ecDNAs"      => "DNA| ",
	"ecRNAs"      => "RNA| ",
	"proteinECs"  => "Prot| ",
	"HTXileum"    => "Ileum| ",
	"HTXrectum"   => "Rectum| ",
	"metabolites" => "",
	"fecalcal"    => "",
	"serology"    => "Sera| "
)
type_actives = Dict(
	"bugs"        => active_bugs,
	"bugkoRNADNAs"=> active_bugkoRNADNAs,
	"ecDNAs"      => active_ecDNAs,
	"ecRNAs"      => active_ecRNAs,
	"proteinECs"  => active_proteinECs,
	"HTXileum"    => active_HTXileum,
	"HTXrectum"   => active_HTXrectum,
	"metabolites" => active_metabolites,
	"fecalcal"    => active_fecalcal,
	"serology"    => active_serology
)
dont_include = Set([
	"bugs-ecDNAs"
])
nice_dataset_names = Dict(
	"bugs"        => "MGX|Species",
	"bugkoRNADNAs"=> "MTX|Species",
	"ecDNAs"      => "MGX|EC",
	"ecRNAs"      => "MTX|EC",
	"proteinECs"  => "MPX|EC",
	"HTXileum"    => "HTX|Ileum",
	"HTXrectum"   => "HTX|Rectum",
	"metabolites" => "MBX",
	"fecalcal"    => "Calprotectin",
	"serology"    => "Serology",
)
type_edgetotals = Dict()

for pr1 in type_prefixes
	for pr2 in type_prefixes
		type1 = pr1[1]
		type2 = pr2[1]

		if type1 == type2 || in("$(type1)-$(type2)", dont_include) || in("$(type2)-$(type1)", dont_include)
			continue
		end

		# Find the association table
		q = ""
		path = joinpath(halla_dir, "$(type1)-$(type2)_residuals_$(model)$(q)_alla")
		if !isdir(path)
			q = "_q25"
			path = joinpath(halla_dir, "$(type1)-$(type2)_residuals_$(model)$(q)_alla")
			if !isdir(path)
				q = "_q05"
				path = joinpath(halla_dir, "$(type1)-$(type2)_residuals_$(model)$(q)_alla")
			end
		end

		if !isdir(path)
			tmp = type1
			type1 = type2
			type2 = tmp

			q = ""
			path = joinpath(halla_dir, "$(type1)-$(type2)_residuals_$(model)$(q)_alla")
			if !isdir(path)
				q = "_q25"
				path = joinpath(halla_dir, "$(type1)-$(type2)_residuals_$(model)$(q)_alla")
				if !isdir(path)
					q = "_q05"
					path = joinpath(halla_dir, "$(type1)-$(type2)_residuals_$(model)$(q)_alla")
				end
			end

			if !isdir(path)
				println("Warning: no output folder for $(type1)-$(type2)!")
				push!(dont_include, "$(type1)-$(type2)")
				continue
			end
		end

		prep1 = type_prefixes[type1]
		prep2 = type_prefixes[type2]
		active1 = type_actives[type1]
		active2 = type_actives[type2]

		# Read the association table
		edges = alla_associations(path, type1, type2, prep1, prep2, active1, active2)
		edges[:dataset1] = nice_dataset_names[type1]
		edges[:dataset2] = nice_dataset_names[type2]

		# Append to the giant list of associations
		if nrow(edges) > 0
			global all_edges
			if nrow(all_edges) > 0
				all_edges = [all_edges; edges]
			else
				all_edges = edges
			end

			type_edgetotals[type1] = get(type_edgetotals, type1, 0) + nrow(edges)
			type_edgetotals[type2] = get(type_edgetotals, type2, 0) + nrow(edges)
		end

		println("$(type1)-$(type2)$(q): $(nrow(edges)) edges")
		push!(dont_include, "$(type1)-$(type2)")
	end
end


save(joinpath(cytoscape_output, "edges.tsv"), all_edges, quotechar=nothing)
save(joinpath(cytoscape_output, "nodes.tsv"), [
	active_bugs; active_bugkoRNADNAs;
	active_ecDNAs; active_ecRNAs; active_proteinECs;
	active_metabolites;
	active_fecalcal;
	active_HTXileum; active_HTXrectum
	active_serology], quotechar=nothing)


println("")
println(" Edge totals")
for x in type_edgetotals
	println("$(x[1]): $(x[2])")
end
