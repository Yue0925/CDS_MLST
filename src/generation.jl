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

    # For each grid size considered 69
    for m in 3:10
        for n in m:10
            for type in ["heuristic", "cplex"]
                fileName = "../res/Grid/" * type * "/instance_m" * string(m) * "x" * "_n" * string(n) * ".txt"
                #if !isfile(fileName)
                    println("-- Generating file " * fileName)

                    if type == "heuristic"
                        isOptimal = true
                        start = time()
                        CDS, nb_leaves, reversed = heuristiqueCDSForGrids(m, n)
                        solveTime = time() - start
                        writeSolutionReversed(fileName, m, n, CDS, nb_leaves, reversed)
                    end

                    if type == "cplex"
                        isOptimal, solveTime, x, Nodes = cplexSolveMLST2(m, n, preProcessingForGrid(m, n))
                        if isOptimal
                            CDS, nb_leaves = variablesToCDS(m, n, x, Nodes)
                            writeSolution(fileName, m, n, CDS, nb_leaves)
                        end
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

    resultFolder = "../res/Grid/"
    resultsArray(resultFolder, resultFolder * "array.tex")
end
