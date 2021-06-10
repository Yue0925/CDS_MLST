# This file contains methods to solve an instance with CPLEX
include("io.jl")

#using CPLEX
#using JuMP
#using Combinatorics

TOL = 0.00001

"""
Solve the relaxed MLST problem for grid with CPLEX (minimum spanning tree problem)

Argument
- m: the number of lines in a grid graph (m>=5)
- n: the number of columns in a grid graph (n>=5)

Return
- status: Optimal if the problem is solved optimally
- x: 4-dimensional variables array such that x[i, j, i_, j_] = 1 if edge {(i, j), (i_,j_)} belongs to a spanning tree
- getsolvetime(m): resolution time in seconds
"""
function cplexSolveMLSTGrid2(m::Int, n::Int)
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
                    end
                    @constraint(M, x[i, j, i_, j_] - x[i_, j_, i, j] == 0) # symetric edge {(i,j), (i_,j_)} = {(i_,j_), (i,j)}
                end
            end
        end
    end

    # No isolated vertex (interior nodes of degree-4)
    for i in 2:m-1
        for j in 2:n-1
             @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 1)
             cost[i, j, i-1, j] = cost[i, j, i+1, j] = cost[i, j, i, j-1] = cost[i, j, i, j+1] = 2/3
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

    # Constraints acyclic for any subsets of vertices
    @constraint(M, [i_ in 2:m, j_ in 2:n], sum(x[1:i_, 1:j_, 1:i_, 1:j_]) == 2*( i_*j_-1 ) )
    #@constraint(M, [i in 1:m-1, j in 1:n-1, i_ in i+1:m, j_ in j+1:n], sum(x[i:i_, j:j_, i:i_, j:j_]) == 2*( (i_-i+1)*(j_-j+1)-1 ) )

    @objective(M, Min, sum(x[i, j, i_, j_]*cost[i, j, i_, j_] for i in 1:m for j in 1:n for i_ in 1:m for j_ in 1:n))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each edge ((Array{VariableRef, 4}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, x

end # function


"""
An integer linear programming for the maximum leaves spanning tree problem.
"""
function cplexSolveMLSTGrid3(m::Int, n::Int)

    # Create the model
    M = Model(CPLEX.Optimizer)

    # x[i, j, i_, j_] = 1 if edge {(i,j), (i_, j_)} belongs to a spanning tree
    @variable(M, x[1:m, 1:n, 1:m, 1:n], Bin)

    # y[i, j] = 1, if node (i, j) is a leaf in the spanning tree
    @variable(M, y[1:m, 1:n], Bin)

    # No buckle
    @constraint(M, [i in 1:m, j in 1:n], x[i, j, i, j] == 0)

    # 4-connexity constraint
    for i in 1:m
        for j in 1:n
            for i_ in 1:m
                for j_ in 1:n
                    if abs(i-i_) + abs(j-j_) >1
                        @constraint(M, x[i, j, i_, j_] == 0)
                    end
                    @constraint(M, x[i, j, i_, j_] - x[i_, j_, i, j] == 0) # symetric edge {(i,j), (i_,j_)} = {(i_,j_), (i,j)}
                end
            end
        end
    end

    # No isolated vertex (interior nodes of degree-4)
    for i in 2:m-1
        for j in 2:n-1
             #@constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 1)
             @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) + 3*y[i, j] <= 4)
        end
    end

    # No isolated vertex (boundary nodes of degree-3)
    i=1
    for j in 2:n-1
        #@constraint(M, (x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 1)
        @constraint(M, (x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) + 2*y[i, j] <= 3)
    end
    i=m
    for j in 2:n-1
        #@constraint(M, (x[i, j, i-1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >=1)
        @constraint(M, (x[i, j, i-1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) + 2*y[i, j] <= 3)
    end
    j=1
    for i in 2:m-1
        #@constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j+1]) >= 1)
        @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j+1]) + 2*y[i, j] <= 3)
    end
    j=n
    for i in 2:m-1
        #@constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1]) >= 1)
        @constraint(M, (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1]) + 2*y[i, j] <= 3)
    end

    # No isolated vertex (corner nodes of degree-2)
    #@constraint(M, (x[1, 1, 2, 1] + x[1, 1, 1, 2]) >= 1)
    @constraint(M, (x[1, 1, 2, 1] + x[1, 1, 1, 2]) + y[1, 1] <= 2)

    #@constraint(M, (x[1, n, 1, n-1] + x[1, n, 2, n]) >= 1)
    @constraint(M, (x[1, n, 1, n-1] + x[1, n, 2, n]) + y[1, n] <= 2)

    #@constraint(M, (x[m, 1, m-1, 1] + x[m, 1, m, 2]) >= 1)
    @constraint(M, (x[m, 1, m-1, 1] + x[m, 1, m, 2]) + y[m, 1] <= 2)

    #@constraint(M, (x[m, n, m-1, n] + x[m, n, m, n-1]) >= 1)
    @constraint(M, (x[m, n, m-1, n] + x[m, n, m, n-1]) + y[m, n] <= 2)

    # Constraints acyclic for any subsets of vertices
    @constraint(M, [i_ in 2:m, j_ in 2:n], sum(x[1:i_, 1:j_, 1:i_, 1:j_]) == 2*( i_*j_-1 ) )
    #@constraint(M, [i in 1:m-1, j in 1:n-1, i_ in i+1:m, j_ in j+1:n], sum(x[i:i_, j:j_, i:i_, j:j_]) == 2*( (i_-i+1)*(j_-j+1)-1 ) )

    @objective(M, Max, sum(y))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each edge ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, y

end # function



"""
First try to solve the minimum spanning tree problem using the acyclic constraint for all subsets of vertices
"""
function cplexSolveMLSTGrid(m::Int, n::Int)
    # Create the model
    M = Model(CPLEX.Optimizer)

    costs = preCalcul(m, n) # edge => cost
    E = length(costs)
    ID = Dict(zip(1:E, collect(keys(costs)))) # id => edge
    vertices = Set{Tuple{Int, Int}}()
    for i in 1:m
        for j in 1:n
            push!(vertices, (i, j))
        end
    end

    @variable(M, x[1:E], Bin) # if the edge x[id] belongs to a spanning tree

    # for all subsets of vertices S, the number of edges in the the induced subgraph is equal to |S|-1
    for subset in collect(powerset(collect(vertices)))
        S = length(subset)
        if S == 0
            continue
        end

        aux = Set{Int}()

        for p in 1:S-1
            (i, j) = subset[p]
            for q in p+1:S
                (i_, j_) = subset[q]
                result = filter((id, e) -> e == (i, j, i_, j_) || (i_, j_, i, j), ID)
                if length(result) != 0
                    #if length(result) >1
                    #    println("Impossible! ", result)
                    #end
                    for id in keys(result) push!(aux, id)  end
                end
            end
        end

        @constraint(M, sum(x[id] for id in collect(aux)) == S-1)
    end

    @objective(M, Min, sum(x[i]*costs[ID[i]] for i in 1:E))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each edge (Array{VariableRef, 1})
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, x
end # function



"""
Given a grid of m lines and n columns, calculate the weight of each edge.

"""
function preCalcul(m::Int, n::Int)
    # edge (i, j)-(i_, j_) => weight
    costs = Dict{Tuple{Int, Int, Int, Int}, Float64}()

    # (interior nodes of degree-4)
    for i in 2:m-1
        for j in 2:n-1
            push!(costs, (i, j, i-1, j) => 2/3)
            push!(costs, (i, j, i+1, j) => 2/3)
            push!(costs, (i, j, i, j-1) => 2/3)
            push!(costs, (i, j, i, j+1) => 2/3)
        end
    end

    # (boundary nodes of degree-3)
    i=1
    for j in 2:n-1
        costs[(i, j, i+1, j)] = 5/6
        costs[(i+1, j, i, j)] = 5/6
        push!(costs, (i, j, i, j-1) => 1)
        push!(costs, (i, j, i, j+1) => 1)
    end
    i=m
    for j in 2:n-1
        costs[(i, j, i-1, j)] = 5/6
        costs[(i-1, j, i, j)] = 5/6
        push!(costs, (i, j, i, j-1) => 1)
        push!(costs, (i, j, i, j+1) => 1)
    end
    j=1
    for i in 2:m-1
        costs[(i, j, i, j+1)] = 5/6
        costs[(i, j+1, i, j)] = 5/6
        push!(costs, (i, j, i-1, j) => 1)
        push!(costs, (i, j, i+1, j) => 1)
    end
    j=n
    for i in 2:m-1
        costs[(i, j, i, j-1)] = 5/6
        costs[(i, j-1, i, j)] = 5/6
        push!(costs, (i, j, i-1, j) => 1)
        push!(costs, (i, j, i+1, j) => 1)
    end

    #  (corner nodes of degree-2)
    costs[(1, 1, 2, 1)] = costs[(1, 1, 1, 2)] = costs[(2, 1, 1, 1)] = costs[(1, 2, 1, 1)] = 3/2
    costs[(1, n, 1, n-1)] = costs[(1, n, 2, n)] = costs[(1, n-1, 1, n)] = costs[(2, n, 1, n)] = 3/2
    costs[(m, 1, m-1, 1)] = costs[(m, 1, m, 2)] = costs[(m-1, 1, m, 1)] = costs[(m, 2, m, 1)] = 3/2
    costs[(m, n, m-1, n)] = costs[(m, n, m, n-1)] = costs[(m-1, n, m, n)] = costs[(m, n-1, m, n)] = 3/2

    for i in 1:m
        for j in 1:n
            for i_ in 1:m
                for j_ in 1:n
                    if haskey(costs, (i_, j_, i, j)) && haskey(costs, (i, j, i_, j_))
                        delete!(costs, (i_, j_, i, j))
                    end
                end
            end
        end
    end

    return costs
end

function VariablesToCDS(m::Int, n::Int, Y::Array{VariableRef, 2})
    # Convert variables
    CDS = Set{Tuple{Int, Int}}()

    for i in 1:m
        for j in 1:n
            if JuMP.value(Y[i, j]) <= TOL
                push!(CDS, (i, j))
            end
        end
    end

    return CDS, m*n-length(CDS)
end


"""
Convert the cplex solution to the corresponding connected dominating set
"""
function VariablesToCDS(m::Int, n::Int, X::Array{VariableRef, 4})
    # Convert variables
    x = zeros(Int64, (m, n, m, n))

    for i in 1:m
        for j in 1:n
            for i_ in 1:m
                for j_ in 1:n
                    if JuMP.value(X[i, j, i_, j_]) > TOL
                        x[i, j, i_, j_] = 1
                    end
                end
            end
        end
    end

    CDS = Set{Tuple{Int, Int}}()

    # interior nodes of degree-4
    for i in 2:m-1
        for j in 2:n-1
            if (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1]+ x[i, j, i, j+1]) >= 2
                push!(CDS, (i, j))
            end
        end
    end

    # boundary nodes of degree-3
    i=1
    for j in 2:n-1
        if (x[i, j, i+1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 2
            push!(CDS, (i, j))
        end
    end
    i=m
    for j in 2:n-1
        if (x[i, j, i-1, j] + x[i, j, i, j-1] + x[i, j, i, j+1]) >= 2
            push!(CDS, (i, j))
        end
    end
    j=1
    for i in 2:m-1
        if (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j+1]) >= 2
            push!(CDS, (i, j))
        end
    end
    j=n
    for i in 2:m-1
        if (x[i, j, i-1, j] + x[i, j, i+1, j] + x[i, j, i, j-1]) >= 2
            push!(CDS, (i, j))
        end
    end

    # corner nodes of degree-2
    if (x[1, 1, 2, 1] + x[1, 1, 1, 2]) >= 2 push!(CDS, (1, 1)) end

    if (x[1, n, 1, n-1] + x[1, n, 2, n]) >= 2 push!(CDS, (1, n)) end

    if (x[m, 1, m-1, 1] + x[m, 1, m, 2]) >= 2 push!(CDS, (m, 1)) end

    if (x[m, n, m-1, n] + x[m, n, m, n-1]) >= 2 push!(CDS, (m, n)) end

    return CDS
end



"""
KRUSKAL algrithm finds a minimum spanning tree

"""
function KRUSKAL(m::Int, n::Int)
    classe = Dict{Tuple{Int, Int}, Int}()
    degrees = Dict{Tuple{Int, Int}, Int}()

    c = 1
    for i in 1:m
        for j in 1:n
            push!(classe, (i, j) => c)
            c += 1
            push!(degrees, (i, j) => 0)
        end
    end

    costs = preCalcul(m, n)
    total_cost = 0

    for ((i, j, i_, j_), val)  in sort(costs; byvalue=true)
        if classe[(i, j)] != classe[(i_, j_)]
            degrees[i, j] += 1
            degrees[i_, j_] += 1
            total_cost += val
            c_ = classe[(i, j)]
            c = classe[(i_, j_)]
            for key in keys(filter(p -> p.second == c_, classe))
                classe[key] = c
            end
        end
    end

    CDS = Set(keys(filter(p -> p.second > 1, degrees)))

    return CDS, m*n-length(CDS)
end # function



m=15
n=6
isOptimal, solveTime, x = cplexSolveMLSTGrid3(m, n)
println("isOptimal: ", isOptimal)

if isOptimal
    CDS, nb_leaves = VariablesToCDS(m, n, x)
    #CDS, nb_leaves = KRUSKAL(m, n)
    displaySolution(m, n, CDS, nb_leaves)
end


#vertices = Set{Tuple{Int, Int}}()
#for i in 1:m for j in 1:n push!(vertices, (i, j)) end end
#collect(powerset(collect(vertices)))
