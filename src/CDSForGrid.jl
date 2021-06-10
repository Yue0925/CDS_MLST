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

    cost = (3*n-7) * 2/3 + 4 * 3/2 + (3*2 + n-2 + ceil(Int, n/3)) * 5/6

    if n%3 == 0
        cost += (2 * (n/3 - 2) + 2)
    elseif n%3 == 1
        cost += (2 * (floor(Int, n/3) - 1) + 1 + 2)
    else
        cost += (2 * (floor(Int, n/3) - 1) + 1)
    end

    return CDS, m*n-length(CDS)
end
