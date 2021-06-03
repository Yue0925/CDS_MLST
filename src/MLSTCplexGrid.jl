# This file contains methods to solve an instance with CPLEX

#using CPLEX
#using JuMP

"""
Solve the MLST problem for grid with CPLEX

Argument
- m: the number of lines in a grid graph
- n: the number of columns in a grid graph

Return
- status: :Optimal if the problem is solved optimally
- x: 4-dimensional variables array such that x[i, j, i_, j_] = 1 if edge {(i, j), (i_,j_)} belongs to a spanning tree
- getsolvetime(m): resolution time in seconds
"""
function cplexSolveMLSTGrid(m::Int, n::Int)
    # the cost of edge {(i,j), (i_,j_)}
    cost = zeros(Float64, (m, n, m, n))

    # Create the model
    M = Model(CPLEX.Optimizer)

    # x[i, j, i_, j_] = 1 if edge {(i,j), (i_, j_)} belongs to a spanning tree
    @variable(M, x[1:m, 1:n, 1:m, 1:n], Bin)

    # No buckle
    @constraint(M, [i in 1:m, j in 1:n], x[i, j, i, j] == 0)

    # 4-connexity constraint
    for i in 1:m
        for j in 1:n
            for i_ in 1:m
                for j_ in 1:n
                    if abs(i-i_) + abs(j-j_) >1
                        @constraint(M, x[i, j, i_, j_] == 0)
                    else
                        @constraint(M, x[i, j, i_, j_] == x[i_, j_, i, j]) # edge {(i,j), (i_,j_)} = {(i_,j_), (i,j)}
                    end
                end
            end
        end
    end

    # No isolated vertex (interior nodes of degree-4)
    for i in 2:m-1
        for j in 2:n-1
             @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 1)
             cost[i, j, i-1, j] = cost[i, j, i+1, j] = cost[i, j, i, j-1] =  cost[i, j, i, j+1] = 2/3
        end
    end

    # No isolated vertex (boundary nodes of degree-3)
    i=1
    for j in 2:n-1
        @constraint(M, (x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 1)
        cost[i, j, i+1, j] = cost[i+1, j, i, j] = 5/6
        cost[i, j, i, j-1] = cost[i, j, i, j+1] = 1
    end
    i=m
    for j in 2:n-1
        @constraint(M, (x[i, j, i-1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >=1)
        cost[i, j, i-1, j] = cost[i-1, j, i, j] = 5/6
        cost[i, j, i, j-1] = cost[i, j, i, j+1] = 1
    end
    j=1
    for i in 2:m-1
        @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j+1]) >= 1)
        cost[i, j, i, j+1] = cost[i, j+1, i, j] = 5/6
        cost[i, j, i-1, j] = cost[i, j, i+1, j] = 1
    end
    j=n
    for i in 2:m-1
        @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1]) >= 1)
        cost[i, j, i, j-1] = cost[i, j-1, i, j] = 5/6
        cost[i, j, i-1, j] = cost[i, j, i+1, j] = 1
    end

    # No isolated vertex (corner nodes of degree-2)
    @constraint(M, (x[1, 1, 2, 1] + x[1, 1, 1, 2]) >= 1)
    cost[1, 1, 2, 1] = cost[1, 1, 1, 2] = cost[2, 1, 1, 1] = cost[1, 2, 1, 1] = 3/2
    @constraint(M, (x[1, n, 1, n-1] + x[1, n, 2, n]) >= 1)
    cost[1, n, 1, n-1] = cost[1, n, 2, n] = cost[1, n-1, 1, n] = cost[2, n, 1, n] = 3/2
    @constraint(M, (x[m, 1, m-1, 1] + x[m, 1, m, 2]) >= 1)
    cost[m, 1, m-1, 1] = cost[m, 1, m, 2] = cost[m-1, 1, m, 1] = cost[m, 2, m, 1] = 3/2
    @constraint(M, (x[m, n, m-1, n] + x[m, n, m, n-1]) >= 1)
    cost[m, n, m-1, n] = cost[m, n, m, n-1] = cost[m-1, n, m, n] = cost[m, n-1, m, n] = 3/2

    # the number of edges in a spanning tree = the number of vertices -1
    @constraint(M, [i in 1:m, j in 1:n, i_ in 1:m, j_ in 1:n], sum(x[i, j, i_, j_]) == 2*(m*n-1))

    @objective(M, Min, sum(x[i, j, i_, j_]*cost[i, j, i_, j_] for i in 1:m for j in 1:n for i_ in 1:m for j_ in 1:n))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each edge
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, round.(Int, JuMP.value.(x))

end # function


function MLSTToCDS(m::Int, n::Int, x::Array{Int64,4})
    CDS = Set{Tuple{Int, Int}}()

    # interior nodes of degree-4
    for i in 2:m-1
        for j in 2:n-1
            if sum(x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1], x[i, j, i, j+1]) > 1
                push!(CDS, (i, j))
            end
        end
    end

    # boundary nodes of degree-3
    i=1
    for j in 2:n-1
        if sum(x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) > 1
            push!(CDS, (i, j))
        end
    end
    i=m
    for j in 2:n-1
        if sum(x[i, j, i-1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) > 1
            push!(CDS, (i, j))
        end
    end
    j=1
    for i in 2:m-1
        if sum(x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j+1]) > 1
            push!(CDS, (i, j))
        end
    end
    j=m
    for i in 2:m-1
        if sum(x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1]) > 1
            push!(CDS, (i, j))
        end
    end

    # corner nodes of degree-2
    if sum(x[1, 1, 2, 1] + x[1, 1, 1, 2]) > 1 push!(CDS, (1, 1)) end

    if sum(x[1, n, 1, n-1] + x[1, n, 2, n]) > 1 push!(CDS, (1, n)) end

    if sum(x[m, 1, m-1, 1] + x[m, 1, m, 2]) > 1 push!(CDS, (m, 1)) end

    if sum(x[m, n, m-1, n] + x[m, n, m, n-1]) > 1 push!(CDS, (m, n)) end

    return CDS
end

"""
KRUSKAL algrithm finds a minimum spanning tree

"""
function KRUSKAL(m::Int, n::Int)


end # function



m=5
n=5
#isOptimal, solveTime, x = cplexSolveMLSTGrid(m, n)
#CDS = MLSTToCDS(m, n, x)
#displaySolution(m, n, CDS)
