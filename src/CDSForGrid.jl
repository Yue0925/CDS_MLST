# This file contains methods to find a connected dominating set for grid graphs. [Mujuni2013]
include("io.jl")

"""
Find a connected dominating set for grid graphs

Argument
- m: the number of lines in a grid graph
- n: the number of columns in a grid graph

Return
- CDS: a connected dominating set containing vertices' coordinates
"""
function CDSForGridGraphs(m::Int, n::Int)
    CDS = Set{Tuple{Int, Int}}()
    for k in 1:n push!(CDS, (2, k)) end # vertices on the second line

    for k in 1:floor(Int, n/3)
        for i in 3:m
            push!(CDS, (i, 3*k-1)) # columns on 3*k-1
        end
    end

    if n%3 == 1 # n ≡ 1 (mod 3)
        for i in 3:m-1
            push!(CDS, (i, n))
        end
    end

    if n%3 == 2 # n ≡ 2 (mod 3)
        for i in 3:m
            push!(CDS, (i, n-1))
        end
    end

    return CDS
end


#m=5
#n=6
#displaySolution(m, n, CDSForGridGraphs(m,n))
