
TOL = 0.00001

"""
Given an integer n, construct an n-hypercube graph

"""
function preprocessingHypercube(n::Int)
    hypercube = ["0", "1"]
    Nodes = Dict{String, Int}()
    count = 1

    for i in 2:n
        hypercube0 = ["0"*v for v in hypercube]
        hypercube1 = ["1"*v for v in hypercube]
        hypercube = vcat(hypercube0, hypercube1)
    end

    Graph = Dict(v => Set{String}() for v in hypercube)

    for v in hypercube
        Nodes[v] = count
        count += 1
        for u in hypercube
            aux = 0
            for i in 1:n
                if v[i] != u[i] aux += 1 end
            end
            if aux == 1 push!(Graph[v], u) end
        end
    end

    return Graph, Nodes
end



function cplexSolveMLSTHypercube(n::Int)

    Graph, Nodes = preprocessingHypercube(n)
    V = 2^n-1

    # Create the model
    M = Model(CPLEX.Optimizer)

    # x[u, v] = 1 if edge {u, v} belongs to a spanning tree
    @variable(M, x[1:2^n, 1:2^n], Bin)

    # y[u] = 1, if node u is a leaf in the spanning tree
    @variable(M, y[1:2^n], Bin)

    # z[u, v, t] = 1 if edge {u, v} is a part of the path from s to t
    @variable(M, 0<= z[1:2^n, 1:2^n], Int)

    # we choose the source node
    s = "0"^n

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
        degree = n
        @constraint(M, sum(x[id_u, Nodes[v]] for v in Graph[u]) + (degree-1)*y[id_u] <= degree)
    end

    @constraint(M, sum(x) == 2*V)

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
Convert the linear programming variables to the corresponding connected dominating set.

"""
function variablesToCDS(n::Int, Y::Array{VariableRef, 1}, Nodes::Dict{String, Int})
    # Convert variables
    CDS = Set{String}()

    for (u, id_u) in Nodes
        if JuMP.value(Y[id_u]) <= TOL
            push!(CDS, u)
        end
    end

    return CDS
end


#n = 4
#isOptimal, solveTime, x, Nodes = cplexSolveMLSTHypercube(n)
#if isOptimal
#    CDS = variablesToCDS(n, x, Nodes)
    #writeSolution(fileName, n, CDS, nb_leaves)
#end
