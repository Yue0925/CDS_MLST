# This file contains methods to solve an instance with CPLEX
include("io.jl")

#using CPLEX
#using JuMP

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



function cplexSolveMLST2(Graph::Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}})
    # node (i, j) => id
    Nodes = Dict{Tuple{Int, Int}, Int}()

    i = 1
    for node in keys(Graph)
        Nodes[node] = i
        i += 1
    end

    # Create the model
    M = Model(CPLEX.Optimizer)

    # x[u, v] = 1 if edge {u, v} belongs to a spanning tree
    @variable(M, x[1:m*n, 1:n*m], Bin)

    # z[u, v, t] = 1 if edge {u, v} is a part of the path from s to t
    @variable(M, z[1:m*n, 1:m*n, 1:m*n], Bin)

    # we choose the node (1, 1) as the source node
    s = (1, 1)

    # connexity constraint
    for (u, id_u) in Nodes
        for (v, id_v) in Nodes
            if !(v in Graph[u])
                @constraint(M, x[id_u, id_v] == 0)
            end
            @constraint(M, x[id_u, id_v] == x[id_v, id_u])
        end
    end

    # there exists a path from source node to all other nodes
    for (t, id_t) in Nodes
        if t == s continue end
        for (u, id_u) in Nodes
            if u == s
                b = 1
            elseif u == t
                b = -1
            else
                b = 0
            end
            # the flows on the path from source node to terminal node can exist if it is a part of the spanning tree
            @constraint(M, [id_v in values(Nodes)], x[id_u, id_v] >= z[id_u, id_v, id_t])
            # the incoming flows = the outgoing flows except the source and terminal nodes
            @constraint(M, sum(z[id_u, id_v, id_t] for id_v in values(Nodes)) == sum(z[id_v, id_u, id_t] for id_v in values(Nodes) )+b )
        end
    end

    @objective(M, Max, sum( (length(Graph[u]) - sum(x[id_u, Nodes[v]] for v in Graph[u]) )/(length(Graph[u])-1) for (u, id_u) in Nodes) )

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each edge ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, x, Nodes

end # function






"""
Given a graph, solve the MLST problem using the flow formulation.
"""
function cplexSolveMLST(Graph::Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}})
    # node (i, j) => id
    Nodes = Dict{Tuple{Int, Int}, Int}()

    i = 1
    for node in keys(Graph)
        Nodes[node] = i
        i += 1
    end

    # Create the model
    M = Model(CPLEX.Optimizer)

    # x[u, v] = 1 if edge {u, v} belongs to a spanning tree
    @variable(M, x[1:m*n, 1:n*m], Bin)

    # y[u] = 1, if node u is a leaf in the spanning tree
    @variable(M, y[1:m*n], Bin)

    # z[u, v, t] = 1 if edge {u, v} is a part of the path from s to t
    @variable(M, z[1:m*n, 1:m*n, 1:m*n], Bin)

    # we choose the node (1, 1) as the source node
    s = (1, 1)

    # connexity constraint
    for (u, id_u) in Nodes
        for (v, id_v) in Nodes
            if !(v in Graph[u])
                @constraint(M, x[id_u, id_v] == 0)
            end
            @constraint(M, x[id_u, id_v] == x[id_v, id_u])
        end
    end

    # there exists a path from source node to all other nodes
    for (t, id_t) in Nodes
        if t == s continue end
        for (u, id_u) in Nodes
            if u == s
                b = 1
            elseif u == t
                b = -1
            else
                b = 0
            end
            # the flows on the path from source node to terminal node can exist if it is a part of the spanning tree
            @constraint(M, [id_v in values(Nodes)], x[id_u, id_v] >= z[id_u, id_v, id_t])
            # the incoming flows = the outgoing flows except the source and terminal nodes
            @constraint(M, sum(z[id_u, id_v, id_t] for id_v in values(Nodes)) == sum(z[id_v, id_u, id_t] for id_v in values(Nodes) )+b )
        end
    end

    # the degree of vertex u in a spanning tree is no greater than the degree in the original graph
    for (u, id_u) in Nodes
        degree = length(Graph[u])
        @constraint(M, sum(x[id_u, Nodes[v]] for v in Graph[u]) + (degree-1)*y[id_u] <= degree)
    end

    @objective(M, Max, sum(y))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each node ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, y, Nodes

end # function


"""
Given a grid of m lines n columns, return the corresponding graph of a dictionnary structure
"""
function preProcessingForGrid(m::Int, n::Int)
    # a dictionnary structure for a graph, each node (i, j) => a set of neighbors
    graph = Dict((i, j) => Set{Tuple{Int, Int}}() for i in 1:m for j in 1:n)

    # (interior nodes of degree-4)
    for i in 2:m-1
        for j in 2:n-1
            push!(graph[(i, j)], (i-1, j))
            push!(graph[(i, j)], (i+1, j))
            push!(graph[(i, j)], (i, j-1))
            push!(graph[(i, j)], (i, j+1))
        end
    end

    # (boundary nodes of degree-3)
    i=1
    for j in 2:n-1
        push!(graph[(i, j)], (i+1, j))
        push!(graph[(i, j)], (i, j-1))
        push!(graph[(i, j)], (i, j+1))
    end
    i=m
    for j in 2:n-1
        push!(graph[(i, j)], (i-1, j))
        push!(graph[(i, j)], (i, j-1))
        push!(graph[(i, j)], (i, j+1))
    end
    j=1
    for i in 2:m-1
        push!(graph[(i, j)], (i, j+1))
        push!(graph[(i, j)], (i-1, j))
        push!(graph[(i, j)], (i+1, j))
    end
    j=n
    for i in 2:m-1
        push!(graph[(i, j)], (i, j-1))
        push!(graph[(i, j)], (i-1, j))
        push!(graph[(i, j)], (i+1, j))
    end

    #  (corner nodes of degree-2)
    push!(graph[(1, 1)], (2, 1))
    push!(graph[(1, 1)], (1, 2))
    push!(graph[(1, n)], (1, n-1))
    push!(graph[(1, n)], (2, n))
    push!(graph[(m, 1)], (m-1, 1))
    push!(graph[(m, 1)], (m, 2))
    push!(graph[(m, n)], (m-1, n))
    push!(graph[(m, n)], (m, n-1))

    return graph
end


"""
Given a grid of m lines and n columns, calculate the cost of each edge.

"""
function calculateCostsInGrid(m::Int, n::Int)
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


"""
Convert the linear programming variables to the corresponding connected dominating set.

"""
function variablesToCDS(m::Int, n::Int, Y::Array{VariableRef, 1}, Nodes::Dict{Tuple{Int, Int}, Int})
    # Convert variables
    CDS = Set{Tuple{Int, Int}}()

    for (u, id_u) in Nodes
        if JuMP.value(Y[id_u]) <= TOL
            push!(CDS, u)
        end
    end

    return CDS, m*n-length(CDS)
end


function variablesToCDS(m::Int, n::Int, X::Array{VariableRef, 2}, Nodes::Dict{Tuple{Int, Int}, Int})
    # Convert variables
    degrees = Dict((i, j) => 0 for i in 1:m for j in 1:n)

    for (u, id_u) in Nodes
        for (v, id_v) in Nodes
            if JuMP.value(X[id_u, id_v]) > TOL
                degrees[u] += 1
            end
        end
    end

    CDS = Set(keys(filter(p -> p.second > 1, degrees)))

    return CDS, m*n-length(CDS)
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

    costs = calculateCostsInGrid(m, n)
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



m=5
n=5
isOptimal, solveTime, x, Nodes = cplexSolveMLST(preProcessingForGrid(m, n))
println("isOptimal: ", isOptimal)

if isOptimal
    CDS, nb_leaves = variablesToCDS(m, n, x, Nodes)
    #CDS, nb_leaves = KRUSKAL(m, n)
    displaySolution(m, n, CDS, nb_leaves)
end
