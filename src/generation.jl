include("CDSForGrid.jl")
include("MLSTCplexGrid.jl")



"""
Generate all the instances and a CDS from [Mujuni2013]

Remark: a solution file is generated only if the corresponding output file does not already exist
"""
function generateResultsGrid()

    isOptimal = false
    solveTime = -1
    CDS = Set{Tuple{Int, Int}}()

    # For each grid size considered
    for m in 5:10
        for n in 5:10
            for type in ["CDS", "PL", "KRUSKAL"]
                fileName = "../res/" * type * "/Grid/instance_m" * string(m) * "x" * "_n" * string(n) * ".txt"
                #if !isfile(fileName)
                    println("-- Generating file " * fileName)

                    if type == "CDS"
                        isOptimal = true
                        start = time()
                        CDS, nb_leaves = CDSForGridGraphs(m, n)
                        solveTime = time() - start
                        writeSolution(fileName, m, n, CDS, nb_leaves)
                    end

                    if type == "PL"
                        isOptimal, solveTime, x, Nodes = cplexSolveMLST(m, n, preProcessingForGrid(m, n))
                        if isOptimal
                            CDS, nb_leaves = variablesToCDS(m, n, x, Nodes)
                            writeSolution(fileName, m, n, CDS, nb_leaves)
                        end
                    end

                    if type == "KRUSKAL"
                        isOptimal = true
                        start = time()
                        CDS, nb_leaves = KRUSKAL(m, n)
                        solveTime = time() - start
                        writeSolution(fileName, m, n, CDS, nb_leaves)
                    end

                    fout = open(fileName, "a")
                    println(fout, "solveTime = ", solveTime)
                    println(fout, "isOptimal = ", isOptimal)
                    close(fout)

                    print("Solved problem: " * type * ", optimal: ", isOptimal, ", time: "* string(round(solveTime, sigdigits=2)) * "s\n")
                    if isOptimal
                        print(" |CDS| = ", length(CDS))
                    end
                    print("\n")
                #end
            end

        end
    end
end
