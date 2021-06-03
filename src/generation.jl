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
    for m in 5:20
        for n in 5:20
            for type in ["CDS", "MLST"]
                fileName = "../res/" * type * "/Grid/instance_m" * string(m) * "x" * "_n" * string(n) * ".txt"
                #if !isfile(fileName)
                    println("-- Generating file " * fileName)

                    if type == "CDS"
                        isOptimal = true
                        start = time()
                        CDS = CDSForGridGraphs(m, n)
                        solveTime = time() - start
                        writeSolution(fileName, m, n, CDS)
                    end

                    if type == "MLST"
                        isOptimal, solveTime, x = cplexSolveMLSTGrid(m, n)
                        CDS = MLSTToCDS(m, n, x)
                        writeSolution(fileName, m, n, CDS)
                    end

                    fout = open(fileName, "a")
                    println(fout, "solveTime = ", solveTime)
                    println(fout, "isOptimal = ", isOptimal)
                    close(fout)

                    println("Solved problem: " * type * ", |CDS| = ", length(CDS), ", optimal: ", isOptimal, ", time: "* string(round(solveTime, sigdigits=2)) * "s\n")
                #end
            end

        end
    end
end
