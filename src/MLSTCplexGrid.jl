# This file contains methods to solve an instance with CPLEX

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


"""
Given a graph, solve the MLST problem using the multi-commodity flow formulation (MCF).
"""
function cplexSolveMLST_MCF(m::Int, n::Int, Graph::Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}})
    # node (i, j) => id
    Nodes = Dict{Tuple{Int, Int}, Int}()
    V = m*n-1

    i = 1
    for node in keys(Graph)
        Nodes[node] = i
        i += 1
    end

    # Create the model
    M = Model(CPLEX.Optimizer)
    # set cplex solver parameters
    set_optimizer_attribute(M, "CPX_PARAM_MIPINTERVAL", 10) # node interval
    set_optimizer_attribute(M, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    set_optimizer_attribute(M, "CPX_PARAM_EPAGAP", 0.99) #absolute gap

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

    for (v, id_v) in Nodes
        for (t, id_t) in Nodes
            if t == s continue end
            if v == s
                b = -1
            elseif v == t
                b = 1
            else
                b = 0
            end
            # the flows on the path from source node to terminal node can exist if it is a part of the spanning tree
            @constraint(M, [id_u in values(Nodes)], x[id_u, id_v] >= z[id_u, id_v, id_t])
            # the incoming flows = the outgoing flows except the source and terminal nodes
            @constraint(M, sum(z[id_u, id_v, id_t] for id_u in values(Nodes)) == sum(z[id_v, id_u, id_t] for id_u in values(Nodes) )+b )
        end
    end

    # the degree of vertex u in a spanning tree is no greater than the degree in the original graph
    for (u, id_u) in Nodes
        degree = length(Graph[u])
        if u==(1, 1) || u==(1, n) || u==(m, 1) || u==(m, n) # the 4 vertices in corner are leaves for sure
            @constraint(M, y[id_u] == 1)
        end
        @constraint(M, sum(x[id_u, Nodes[v]] for v in Graph[u]) + (degree-1)*y[id_u] <= degree)
        @constraint(M, sum(y[Nodes[v]] for v in Graph[u]) <= degree - 1)
    end

    @constraint(M, sum(x) == 2*V)

    @objective(M, Max, sum(y))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    nb_nodes = node_count(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each node ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, y, Nodes, nb_nodes

end # function



"""
Given a graph, solve the MLST problem using the single-commodity flow formulation (SCF).

"""
function cplexSolveMLST_SCF(m::Int, n::Int, Graph::Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}})

    # node (i, j) => id
    Nodes = Dict{Tuple{Int, Int}, Int}()
    V = m*n-1

    i = 1
    for node in keys(Graph)
        Nodes[node] = i
        i += 1
    end

    # Create the model
    M = Model(CPLEX.Optimizer)
    # set cplex solver parameters
    set_optimizer_attribute(M, "CPX_PARAM_MIPINTERVAL", 10) # node interval
    set_optimizer_attribute(M, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    set_optimizer_attribute(M, "CPX_PARAM_EPAGAP", 0.99) #absolute gap
    #set_optimizer_attribute(M, "CPX_PARAM_TILIM", 1000) # time limit

    # x[u, v] = 1 if edge {u, v} belongs to a spanning tree
    @variable(M, x[1:m*n, 1:n*m], Bin)

    # y[u] = 1, if node u is a leaf in the spanning tree
    @variable(M, y[1:m*n], Bin)

    # z[u, v, t] = 1 if edge {u, v} is a part of the path from s to t
    @variable(M, 0<= z[1:m*n, 1:m*n], Int)

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

    for (u, id_u) in Nodes
        if u == s
            b = V
        else
            b = -1
        end
        # the flows on the path from source node to terminal node can exist if it is a part of the spanning tree
        @constraint(M, [id_v in values(Nodes)], x[id_u, id_v]*V >= z[id_u, id_v])
        # the incoming flows = the outgoing flows except the source and terminal nodes
        @constraint(M, sum(z[id_u, id_v] for id_v in values(Nodes)) == sum(z[id_v, id_u] for id_v in values(Nodes) )+b )
    end


    # the degree of vertex u in a spanning tree is no greater than the degree in the original graph
    for (u, id_u) in Nodes
        degree = length(Graph[u])
        if u==(1, 1) || u==(1, n) || u==(m, 1) || u==(m, n) # the 4 vertices in corner are leaves for sure
            @constraint(M, y[id_u] == 1)
        end
        @constraint(M, sum(x[id_u, Nodes[v]] for v in Graph[u]) + (degree-1)*y[id_u] <= degree)
        @constraint(M, sum(y[Nodes[v]] for v in Graph[u]) <= degree - 1)
    end

    @constraint(M, sum(x) == 2*V)

    @objective(M, Max, sum(y))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    nb_nodes = node_count(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each node ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, y, Nodes, nb_nodes

end # function



"""
Using the MTZ subtour-elimination constraints.
"""
function cplexSolveMLST_MTZ(m::Int, n::Int, Graph::Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}})
    # node (i, j) => id
    Nodes = Dict{Tuple{Int, Int}, Int}()
    V = m*n

    i = 1
    for node in keys(Graph)
        Nodes[node] = i
        i += 1
    end

    # Create the model
    M = Model(CPLEX.Optimizer)
    # set cplex solver parameters
    set_optimizer_attribute(M, "CPX_PARAM_MIPINTERVAL", 10) # node interval
    set_optimizer_attribute(M, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    set_optimizer_attribute(M, "CPX_PARAM_EPAGAP", 0.99) # absolute gap

    # x[u, v] = 1 if edge {u, v} belongs to a spanning tree
    @variable(M, x[1:m*n, 1:n*m], Bin)

    # y[u] = 1, if node u is a leaf in the spanning tree
    @variable(M, y[1:m*n], Bin)

    #
    @variable(M, z[1:m*n], Int)


    # connexity constraint
    for (u, id_u) in Nodes
        for (v, id_v) in Nodes
            if !(v in Graph[u])
                @constraint(M, x[id_u, id_v] == 0)
                @constraint(M, x[id_v, id_u] == 0)
            else
                @constraint(M, x[id_u, id_v] + x[id_v, id_u] <=1)
            end
        end
    end


    # we choose the node (1, 1) as the source node
    s = (1, 1)
    id_s = Nodes[s]
    @constraint(M, z[id_s] == 0)

    # at least 1 out from source
    @constraint(M, sum(x[id_s, Nodes[v]] for v in Graph[s]) >= 1)

    for (u, id_u) in Nodes
        if u != s
            @constraint(M, 1 <= z[id_u] <= V)
            # for all other nodes, one incoming arc
            @constraint(M, sum(x[Nodes[v], id_u] for v in Graph[u]) == 1)

            for (v, id_v) in Nodes
                if v != s && v != u
                    @constraint(M, z[id_u] - z[id_v] + (V-1)*x[id_u, id_v] <= V-2)
                end
            end

        end
    end


    # the degree of vertex u in a spanning tree is no greater than the degree in the original graph
    for (u, id_u) in Nodes
        degree = length(Graph[u])
        if u==(1, 1) || u==(1, n) || u==(m, 1) || u==(m, n) # the 4 vertices in corner are leaves for sure
            @constraint(M, y[id_u] == 1)
        end
        @constraint(M, sum(x[id_u, Nodes[v]] for v in Graph[u]) + sum(x[Nodes[v], id_u] for v in Graph[u]) + (degree-1)*y[id_u] <= degree)
        @constraint(M, sum(y[Nodes[v]] for v in Graph[u]) <= degree - 1)
    end


    @objective(M, Max, sum(y))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    nb_nodes = node_count(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each node ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, y, Nodes, nb_nodes

end # function




"""
Using the lifed version of the MTZ subtour-elimination constraints.
"""
function cplexSolveMLST_LMTZ(m::Int, n::Int, Graph::Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}})
    # node (i, j) => id
    Nodes = Dict{Tuple{Int, Int}, Int}()
    V = m*n

    i = 1
    for node in keys(Graph)
        Nodes[node] = i
        i += 1
    end

    # Create the model
    M = Model(CPLEX.Optimizer)
    # set cplex solver parameters
    set_optimizer_attribute(M, "CPX_PARAM_MIPINTERVAL", 10) # node interval
    set_optimizer_attribute(M, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    set_optimizer_attribute(M, "CPX_PARAM_EPAGAP", 0.99) # absolute gap

    # x[u, v] = 1 if edge {u, v} belongs to a spanning tree
    @variable(M, x[1:m*n, 1:n*m], Bin)

    # y[u] = 1, if node u is a leaf in the spanning tree
    @variable(M, y[1:m*n], Bin)

    #
    @variable(M, z[1:m*n], Int)


    # connexity constraint
    for (u, id_u) in Nodes
        for (v, id_v) in Nodes
            if !(v in Graph[u])
                @constraint(M, x[id_u, id_v] == 0)
                @constraint(M, x[id_v, id_u] == 0)
            else
                @constraint(M, x[id_u, id_v] + x[id_v, id_u] <=1)
            end
        end
    end


    # we choose the node (1, 1) as the source node
    s = (1, 1)
    id_s = Nodes[s]
    @constraint(M, z[id_s] == 0)

    # at least 1 out from source
    @constraint(M, sum(x[id_s, Nodes[v]] for v in Graph[s]) >= 1)

    for (u, id_u) in Nodes
        if u != s
            @constraint(M, 1 <= z[id_u] <= V)
            # for all other nodes, one incoming arc
            @constraint(M, sum(x[Nodes[v], id_u] for v in Graph[u]) == 1)

            for (v, id_v) in Nodes
                if v != s && v != u
                    # lifted
                    @constraint(M, z[id_u] - z[id_v] + (V-1)*x[id_u, id_v] + + (V-3)*x[id_v, id_u] <= V-2)
                end
            end

        end
    end


    # the degree of vertex u in a spanning tree is no greater than the degree in the original graph
    for (u, id_u) in Nodes
        degree = length(Graph[u])
        if u==(1, 1) || u==(1, n) || u==(m, 1) || u==(m, n) # the 4 vertices in corner are leaves for sure
            @constraint(M, y[id_u] == 1)
        end
        @constraint(M, sum(x[id_u, Nodes[v]] for v in Graph[u]) + sum(x[Nodes[v], id_u] for v in Graph[u]) + (degree-1)*y[id_u] <= degree)
        @constraint(M, sum(y[Nodes[v]] for v in Graph[u]) <= degree - 1)
    end
    

    @objective(M, Max, sum(y))

    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(M)

    nb_nodes = node_count(M)

    # Return:
    # 1 - true if an optimum is found (type: Bool)
    # 2 - the resolution time (type Float64)
    # 3 - the value of each node ((Array{VariableRef, 2}))
    return JuMP.primal_status(M) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, y, Nodes, nb_nodes

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
