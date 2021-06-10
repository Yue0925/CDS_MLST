# This file contains functions related to reading, writing and displaying a grid and experimental results

"""
Display a CDS solution for a grid graph

Arguments
- fout: the output stream (usually an output file)
- m: the number of lines in a grid graph
- n: the number of columns in a grid graph
- CDS: a connected dominating set containing vertices' coordinates
"""
function displaySolution(m::Int, n::Int, CDS::Set{Tuple{Int, Int}}, value::Float64)
    for i in 1:m
        for j in 1:n
            if (i, j) in CDS
                print('⚫')
            else
                print('⚪')
            end
        end
        print('\n')
    end
    println("A Grid graph of size ", m, "×", n, " , with γ_c(G)=", length(CDS), ", the weight of the spanning tree ", value)
end


"""
Write a solution in an output file

Arguments
- fout: the output stream (usually an output file)
- m: the number of lines in a grid graph
- n: the number of columns in a grid graph
- CDS: a connected dominating set containing vertices' coordinates
"""
function writeSolution(outputFile::String, m::Int, n::Int, CDS::Set{Tuple{Int, Int}}, value::Float64)
    # Open the output file
    writer = open(outputFile, "w")

    for i in 1:m
        for j in 1:n
            if (i, j) in CDS
                print(writer, '⚫')
            else
                print(writer, '⚪')
            end
        end
        print(writer, '\n')
    end
    println(writer, "A Grid graph of size ", m, "×", n, " , with γ_c(G)=", length(CDS), ", the weight of the spanning tree ", value)
    close(writer)
end
