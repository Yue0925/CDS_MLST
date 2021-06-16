# This file contains functions related to reading, writing and displaying a grid and experimental results

"""
Display a CDS solution for a grid graph

Arguments
- fout: the output stream (usually an output file)
- m: the number of lines in a grid graph
- n: the number of columns in a grid graph
- CDS: a connected dominating set containing vertices' coordinates
"""
function displaySolution(m::Int, n::Int, CDS::Set{Tuple{Int, Int}}, nb_leaves::Int)
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
    println("A Grid graph of size ", m, "×", n, " , with γ_c(G)=", length(CDS), " , a number of leaves ", nb_leaves)
end


"""
Write a solution in an output file

Arguments
- fout: the output stream (usually an output file)
- m: the number of lines in a grid graph
- n: the number of columns in a grid graph
- CDS: a connected dominating set containing vertices' coordinates
"""
function writeSolution(outputFile::String, m::Int, n::Int, CDS::Set{Tuple{Int, Int}}, nb_leaves::Int)
    # Open the output file
    writer = open(outputFile, "w")
    println(writer, '\"')

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

    println(writer, "A Grid graph of size ", m, "×", n, "\"")
    println(writer, "γ_c = ", length(CDS))
    println(writer, "nb_leaves = ", nb_leaves)
    close(writer)
end


function writeSolutionReversed(outputFile::String, m::Int, n::Int, CDS::Set{Tuple{Int, Int}}, nb_leaves::Int, reversed::Bool)
    if !reversed
        writeSolution(outputFile, m, n, CDS, nb_leaves)
        return
    end

    # Open the output file
    writer = open(outputFile, "w")
    println(writer, '\"')

    for i in 1:m
        for j in 1:n
            if (j, i) in CDS
                print(writer, '⚫')
            else
                print(writer, '⚪')
            end
        end
        print(writer, '\n')
    end

    println(writer, "A Grid graph of size ", m, "×", n, "\"")
    println(writer, "γ_c = ", length(CDS))
    println(writer, "nb_leaves = ", nb_leaves)
    close(writer)
end




"""
Create a latex file which contains an array with the results of the ../res folder.
Each subfolder of the ../res folder contains the results of a resolution method.

Arguments
- outputFile: path of the output file

Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
- Each text file contains a variable "solveTime", a variable "isOptimal"
"""
function resultsArray(resultFolder::String, outputFile::String)

    # Maximal number of files in a subfolder
    maxSize = 0

    # Number of subfolders
    subfolderCount = 0

    # Open the latex output file
    fout = open(outputFile, "w")

    # Print the latex file output
    println(fout, raw"""\documentclass{article}
\usepackage[french]{babel}
\usepackage [utf8] {inputenc} % utf-8 / latin1
\usepackage{multicol}
\setlength{\hoffset}{-18pt}
\setlength{\oddsidemargin}{0pt} % Marge gauche sur pages impaires
\setlength{\evensidemargin}{9pt} % Marge gauche sur pages paires
\setlength{\marginparwidth}{54pt} % Largeur de note dans la marge
\setlength{\textwidth}{481pt} % Largeur de la zone de texte (17cm)
\setlength{\voffset}{-18pt} % Bon pour DOS
\setlength{\marginparsep}{7pt} % Séparation de la marge
\setlength{\topmargin}{0pt} % Pas de marge en haut
\setlength{\headheight}{13pt} % Haut de page
\setlength{\headsep}{10pt} % Entre le haut de page et le texte
\setlength{\footskip}{27pt} % Bas de page + séparation
\setlength{\textheight}{668pt} % Hauteur de la zone de texte (25cm)
\begin{document}""")

    header = raw"""
\begin{center}
\renewcommand{\arraystretch}{1.4}
 \begin{tabular}{l"""

    # Name of the subfolder of the result folder (i.e, the resolution methods used)
    folderName = Array{String, 1}()

    # List of all the instances solved by at least one resolution method
    solvedInstances = Array{String, 1}()

    # For each file in the result folder
    for file in readdir(resultFolder)

        path = resultFolder * file

        # If it is a subfolder
        if isdir(path)

            # Add its name to the folder list
            folderName = vcat(folderName, file)

            subfolderCount += 1
            folderSize = size(readdir(path), 1)

            # Add all its files in the solvedInstances array
            for file2 in filter(x->occursin(".txt", x), readdir(path))
                solvedInstances = vcat(solvedInstances, file2)
            end

            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

    # Only keep one string for each instance solved
    unique(solvedInstances)

    # For each resolution method, add two columns in the array
    for folder in folderName
        header *= "rrr"
    end

    header *= "}\n\t\\hline\n"

    # Create the header line which contains the methods name
    for folder in folderName
        header *= " & \\multicolumn{3}{c}{\\textbf{" * folder * "}}"
    end

    header *= "\\\\\n\\textbf{Instance} "

    # Create the second header line with the content of the result columns
    for folder in folderName
        header *= " & \\textbf{Time (s)} & \\textbf{CDS} & \\textbf{Number of leaves}  "
    end

    header *= "\\\\\\hline\n"

    footer = raw"""\hline\end{tabular}
\end{center}
"""
    println(fout, header)

    # On each page an array will contain at most maxInstancePerPage lines with results
    maxInstancePerPage = 30
    id = 1

    # For each solved files
    for solvedInstance in solvedInstances

        # If we do not start a new array on a new page
        if rem(id, maxInstancePerPage) == 0
            println(fout, footer, "\\newpage")
            println(fout, header)
        end

        # Replace the potential underscores '_' in file names
        print(fout, replace(solvedInstance, "_" => "\\_"))

        # For each resolution method
        for method in folderName

            path = resultFolder * method * "/" * solvedInstance

            # If the instance has been solved by this method
            if isfile(path)

                include(path)

                println(fout, " & ", round(solveTime, digits=2), " & ", γ_c, " & ", nb_leaves)

            # If the instance has not been solved by this method
            else
                println(fout, " & - & - & - ")
            end
        end

        println(fout, "\\\\")

        id += 1
    end

    # Print the end of the latex file
    println(fout, footer)

    println(fout, "\\end{document}")

    close(fout)

end



#resultFolder = "../res/Grid/"
#include(resultFolder*"cplex/instance_m3x_n3.txt")
