using DataFrames
using FileIO
using CSV
using Colors

const allafiles = [#"all_association_results_one_by_one.txt",
                    "associations.txt",
                    #="similarity_table.txt"=#]
const hallafiles = ["Distance_Matrix1.tsv",
                    "Distance_Matrix2.tsv"]


function halla_associations(hallafolder::String)
    isdir(hallafolder) || error("$hallafolder is not a folder")
    for f in hallafiles
        isfile(joinpath(hallafolder, f)) || error("expected $hallafolder to contain '$f'")
    end
    df =  CSV.read(joinpath(hallafolder, "Distance_Matrix1.tsv"))
    x = reshapedist(df, symetric=true)

    df =  CSV.read(joinpath(hallafolder, "Distance_Matrix2.tsv"))
    y = reshapedist(df, symetric=true)

    df = DataFrame([x; y])
    names!(df, [:node1, :node2,:pearson])

    return df
end

function reshapedist(df::DataFrame; symetric::Bool=false)
    c1 = String[]
    c2 = String[]
    val = Float64[]
    n = replace.(String.(names(df)), "_" => " ")
    n = replace.(n, r"^x" => "")
    for i in 1:size(df,1)
        symetric ? endindx = i+1 : endindx = size(df, 1)
        for j in 2:endindx
            push!(c1, df[i,1])
            push!(c2, n[j])
            push!(val, df[i, j])
        end
    end
    return [c1 c2 val]
end

function alla_associations(hallafolder::String, xtype::String, ytype::String, prep1, prep2, xactive, yactive)
    # Note: despite naming, this function may change xactive and yactive to
    # add new nodes
    isdir(hallafolder) || error("$hallafolder is not a folder")
    for f in allafiles
        isfile(joinpath(hallafolder, f)) || error("expected $hallafolder to contain '$f'")
    end

    edges =  CSV.read(joinpath(hallafolder, "associations.txt"), delim='\t', types=[Int, String, String, BigFloat, BigFloat, Float64])
    #edges = edges[[1, 2, 4, 6, 7, 8]]
    names!(edges, [:rank, :node1, :node2, :pvalue, :qvalue, :pearson])
    edges = edges[edges[:node1] .!= edges[:node2], :]
	edges[:name1] = edges[:node1]
	edges[:name2] = edges[:node2]
    edges[:node1] = prep1.*edges[:node1]
    edges[:node2] = prep2.*edges[:node2]
    edges[:abspearson] = abs.(edges[:pearson])
    edges[:edgetype] = xtype * "-" * ytype

    edges[:rank_active] = 99999
    edges[:rank_active1] = 99999
    if !isa(xactive, Int) && !isa(yactive, Int)
        xinteresting = Set(xactive[:nodename][xactive[:ibd] .!= "None"])
        yinteresting = Set(yactive[:nodename][yactive[:ibd] .!= "None"])

        both_active = in.(edges[:node1], xinteresting) .& in.(edges[:node2], yinteresting)
        edges[:rank_active][both_active] = 1:sum(both_active)

        one_active = in.(edges[:node1], xinteresting) .| in.(edges[:node2], yinteresting)
        edges[:rank_active1][one_active] = 1:sum(one_active)
    end

    if !isa(xactive, Int)
        newedgenodes = unique(edges[:node1])
        newedgenodes = newedgenodes[.!in.(newedgenodes, Set(xactive[:nodename]))]
        if length(newedgenodes) > 0
            append!(xactive, DataFrame(
                nodename = newedgenodes,
                prevalence = 0,
                coefActiveCD = 0, coefActiveUC = 0, coefActivenonIBD = 0,
                tvalActiveCD = 0, tvalActiveUC = 0, tvalActivenonIBD = 0,
                pvalActiveCD = 1, pvalActiveUC = 1, pvalActivenonIBD = 1,
                qvalActiveCD = 1, qvalActiveUC = 1, qvalActivenonIBD = 1,
                minQ = 1,
                nicename = map(x->x[length(prep1):end], newedgenodes),
                nodetype = xtype,
                ibd = "None"
            ))
        end
    end
    if !isa(yactive, Int)
        newedgenodes = unique(edges[:node2])
        newedgenodes = newedgenodes[.!in.(newedgenodes, Set(yactive[:nodename]))]
        if length(newedgenodes) > 0
            append!(yactive, DataFrame(
                nodename = newedgenodes,
                prevalence = 0,
                coefActiveCD = 0, coefActiveUC = 0, coefActivenonIBD = 0,
                tvalActiveCD = 0, tvalActiveUC = 0, tvalActivenonIBD = 0,
                pvalActiveCD = 1, pvalActiveUC = 1, pvalActivenonIBD = 1,
                qvalActiveCD = 1, qvalActiveUC = 1, qvalActivenonIBD = 1,
                minQ = 1,
                nicename = map(x->x[length(prep2):end], newedgenodes),
                nodetype = ytype,
                ibd = "None"
            ))
        end
    end

    return edges
end



function nodesactive(infile::String, ntype::String; minq::Float64=0.2, sep::Char='\t', prep::String="")
    isfile(infile) || error("no file at $infile")
    df =  CSV.read(infile, delim=sep)
    df[1] = replace.(df[1], "_" => " ")
    nm = names(df)
    nm[1] = :nodename
    names!(df, nm)
    df[:nicename] = df[:nodename]
    df[:nodetype] = ntype
    df[:nodename] = prep.*df[:nodename]

	df[:nodename] = replace.(df[:nodename], ";" => " ")

    cs = String[]
    for i in 1:size(df, 1)
        CD = df[i, :coefActiveCD]
        UC = df[i, :coefActiveUC]
        qCD = df[i, :qvalActiveCD]
        qUC = df[i, :qvalActiveUC]

        if qCD < minq && qUC < minq
            if CD < 0 && UC < 0
                c = "nonIBD"
            elseif CD > 0 && UC > 0
                c = "UCCD"
            else
                c = "None"
            end
        elseif qCD < minq
            CD < 0 ? c = "nonIBD" : c = "CD"
        elseif qUC < minq
            UC < 0 ? c = "nonIBD" : c = "UC"
        else
            c = "None"
        end
        push!(cs, c)
    end

    df[:ibd] = cs
    return df
end

function makenodeactive(name::String, nicename::String, ntype::String, ibd::String)
    df = DataFrame()
    df[:nodename] = name

    df[:prevalence] = 1
    df[:coefActiveCD] = 0
    df[:coefActiveUC] = 0
    df[:coefActivenonIBD] = 0
    df[:tvalActiveCD] = 0
    df[:tvalActiveUC] = 0
    df[:tvalActivenonIBD] = 0
    df[:pvalActiveCD] = 0
    df[:pvalActiveUC] = 0
    df[:pvalActivenonIBD] = 0
    df[:qvalActiveCD] = 0
    df[:qvalActiveUC] = 0
    df[:qvalActivenonIBD] = 0
    df[:minQ] = 0

    df[:nicename] = nicename
    df[:nodetype] = ntype
    df[:ibd] = ibd

    return df
end
