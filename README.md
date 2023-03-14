# Tessellation

INPUT FILE FORMAT: Input files are .vtk files with 3D particles positions only.

[The source code will be made available once the corresponding paper is published.]

MAIN CODE FORMAT: At the beginning of the script, there is an option to define the input file name, input file extension, output file name, and output file extension. If the search algorithm is already performed, then neighbor_input should be defined with appropriate file name (default: "NeighborTest.txt"). Otherwise it should be kept as "NULL". Output files are: 1) capsomer file; and 2) outline file. Both can be opened in Paraview. Capsomer file should be treated as "surface" (in Paraview) and Outline file should be treated as "featured edges" (in Paraview).

Note: All the particle positions are named as figure name according to the paper. 
