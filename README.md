# Tessellation

Please, run the source code from SourceCode/MainFunc.cpp

INPUT FILE FORMAT: Input files are .vtk files with 3D particles positions only. For, figures 12-16a, and 18 are named in the format "Figure(FigureNumber)_(Number of repeated particles for periodic boundary condition)particles_(Number of non repeating particles).vtk". The rest of the files have 1000 particles as both non-repeated and total number of particles.

MAIN CODE FORMAT: At the beginning of the script, there is an option to define the input file name, input file extension, output file name, and output file extension. If the search algorithm is already performed, then neighbor_input should be defined with appropriate file name (default: "NeighborTest.txt"). Otherwise it should be kept as "NULL". Output files are: 1) capsomer file; and 2) outline file. Both can be opened in Paraview. Capsomer file should be treated as "surface" (in Paraview) and Outline file should be treated as "featured edges" (in Paraview).

Note: All the particle positions are named as figure name according to the paper. 
