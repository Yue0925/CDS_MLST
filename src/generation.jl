include("io.jl")
include("CDSForGrid.jl")
include("MLSTCplexGrid.jl")
include("resolutionHypercube.jl")



"""
Generate all the instances and a CDS from [Mujuni2013]

Remark: a solution file is generated only if the corresponding output file does not already exist
"""
function generateResultsGrid()

    isOptimal = false
    solveTime = -1
    CDS = Set{Tuple{Int, Int}}()

    # For each grid size considered 69
    for m in 4:10
        for n in m:10
            for type in ["heuristic", "cplex"]
                fileName = "../res/Grid/LIP2/" * type * "/instance_m" * string(m) * "x" * "n" * string(n) * ".txt"

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
                        logFile = "../res/Grid/LIP2/" * type * "/log/instance_m" * string(m) * "x" * "n" * string(n) * ".txt"
                        fileio = open(logFile, "w")
                        tempout = stdout # save stream
                        redirect_stdout(fileio) # redirect to fileio
                        # run solver
                        isOptimal, solveTime, x, Nodes, nb_nodes = cplexSolveMLST2(m, n, preProcessingForGrid(m, n))
                        println("the number of Nodes = ", nb_nodes)
                        close(fileio)
                        redirect_stdout(tempout) #revert back

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

    resultFolder = "../res/Grid/LIP3/"
    resultsArray(resultFolder, resultFolder * "array.tex")
end



function generateResultsHypercube()

    isOptimal = false
    solveTime = -1
    CDS = Set{String}()

    # For each grid size considered 69
    for n in 4:10
        for type in ["heuristic", "cplex"]
            fileName = "../res/Hypercube/" * type * "/Q_" * string(n)* ".txt"
            #if !isfile(fileName)
                println("-- Generating file " * fileName)

                if type == "heuristic"
                    solveTime = 0.0
                    isOptimal = true
                    writer = open(fileName, "w")
                    println(writer, "\" A ", n, "-hypercube graph ", "\"")
                    println(writer, "γ_c = ", 2^(n-2)+2)
                    println(writer, "MLST = ", 2^n - (2^(n-2)+2))
                    close(writer)
                end

                if type == "cplex"
                    isOptimal, solveTime, x, Nodes = cplexSolveMLSTHypercube(n)
                    if isOptimal
                        CDS = variablesToCDS(n, x, Nodes)
                        writeSolution(fileName, n, CDS)
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

    resultFolder = "../res/Hypercube/"
    resultsArray(resultFolder, resultFolder * "array.tex")
end
