#=
By Kevin Bonham

Here is the julia code used to filter mechanistic correlations tables for import
into cytoscape. For each metabolite, I first extract the base name (eg
"HILIC-neg-lithocholate" becomes "lithocholate"). Then, since each metabolite :
feature pair may have  multiple values, I identify the row index for the value
that has the largest absolute rho value, provided that value is significant.
=#

using DataFrames
using LightGraphs
using MetaGraphs

datadir = expanduser("~/Dropbox (Huttenhower Lab)/HMP2/analysis/mechanistic_correlations/")

# Bugs=>Metabolites
t = readtable("$datadir/bugs-metabolites_subject_corr.tsv")
t = t[t[:q_value] .< 0.2, :]

sort!(t, cols=[:q_value, :rho])

metabs = t[:Y]
metabs = match.(r"^\w+\-\w+\-(.+)$", metabs)
t[:Y] = [String(x.captures[1]) for x in metabs]

X = Set(t[:X])
Y = Set(t[:Y])

# species => metabolite => (rowindex, rho)
metabs = Dict(x => Dict(y=>(0,0.) for y in Y) for x in (X))

# iterate through rows, add if absolute rho > what's there currently
for i in 1:size(t, 1)
    sp = t[i, :X]
    met = t[i, :Y]
    ρ = t[i, :rho]

    if abs(ρ) > abs(metabs[sp][met][1])
        metabs[sp][met] = (i, ρ)
    end
end

# get row indicies
idx = [metabs[x][y][1] for x in X, y in Y if metabs[x][y] != (0, 0.)]
# subset df
filt = t[idx, :]


# EC=>Metabolites
t = readtable("$datadir/ecDNAs-metabolites_subject_corr.tsv")



# Try Graphs
using PlotRecipes

sp = unique(filt[:X])
sp = Dict(sp[i]=>i for i in 1:length(sp))

met = unique(filt[:Y])
met = Dict(met[i]=>i+length(keys(sp)) for i in 1:length(met))


G = MetaGraph(length(keys(sp))+length(keys(met)))

for i in 1:size(filt, 1)
    s = filt[i, :X]
    m = filt[i, :Y]
    ρ = filt[i, :rho]
    q = filt[i, :q_value]
    xcd = filt[i, :X_coefCD]
    xuc = filt[i, :X_coefUC]
    ycd = filt[i, :Y_coefCD]
    yuc = filt[i, :Y_coefUC]

    set_prop!(G, sp[s], :name, s)
    set_prop!(G, sp[s], :coefCD, xcd)
    set_prop!(G, sp[s], :coefUC, xuc)
    set_prop!(G, met[m], :name, m)
    set_prop!(G, met[m], :coefCD, ycd)
    set_prop!(G, met[m], :coefUC, yuc)

    # edge
    add_edge!(G, Edge(sp[s], met[m]))
    set_prop!(G, Edge(sp[s], met[m]), :weight, abs(ρ))
    set_prop!(G, Edge(sp[s], met[m]), :color, ρ < 0 ? colorant"darkred" : colorant"steelblue")
    set_prop!(G, Edge(sp[s], met[m]), :q_value, q)
end

clibrary(:colorbrewer)

graphplot(G,
    node_weights=ones(size(G)[1])./2,
    marker_z=[get_prop(G, x, :coefCD) for x in 1:size(G)[1]],
    marker=(:BrBG, [x <= length(sp) ? :circle : :rect for x in 1:size(G)[1]]),
    # names=[get_prop(G, x, :name) for x in 1:size(G)[1]],
    line=(3, 0.5, :blues),
    size=(1200,1200))

savefig(expanduser("~/Desktop/test3.png"))

g = MetaGraph(5)

add_edge!(g, 1,2)
add_edge!(g, 1,3)
add_edge!(g, 3,4)

set_prop!(g, Edge(1,2), :weight, 2.)
set_prop!(g, Edge(1,3), :weight, 3.)
set_prop!(g, Edge(3,4), :weight, 4.)

graphplot(g, line=([get_prop(g, e, :weight) for e in edges(g)], 0.5, :BrBG))
