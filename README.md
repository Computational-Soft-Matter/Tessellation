# Tessellation

INPUT FILE FORMAT: Input files are .vtk files with 3D particles positions only.

MAIN CODE FORMAT: At the beginning of the script, there is option to define the input file name, input file extension, output file name and extension. If the search algorithm is already performed, then neighbor_input should be defined with appropriate file name (default: "bigtest1.txt"). Otherwise it should be kept as "NULL". Output file is of two types- capsomer file and outline file. Both should be opened with paraview (or any visualization kit that is compatible with vtk). Capsomer file should be treated as "surface"(in paraview) and Outline file should be treated as "featured edges"(in paraview).

Note: All the particle positions are named as figure name according to the paper. 
