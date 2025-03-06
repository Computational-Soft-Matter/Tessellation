# Tessellation

Description and explanation in detail about the different sections of the source code are discussed in "https://www.sciencedirect.com/science/article/pii/S0010465524000894". Please follow this journal article. 

To run the tessellation code, please run the source code in SourceCode/MainFunc.cpp

INPUT FILE FORMAT
Input files are files with 3D particles positions only.
Input files for the figures generated in the manuscript are provided as examples and named according to the manuscript figure number.
For figures 12 to 16a and 18, the input files are named in the format "Figure(FigureNumber)_(Number of repeated particles for periodic boundary condition)particles_(Number of non repeating particles).vtk". The rest of the files have 1000 particles as both non-repeated and total number of particles.

MAIN CODE FORMAT: At the beginning of the script, there is an option to define the input file name, input file extension, output file name, and output file extension. If the search algorithm is already performed, then neighbor_input should be defined with appropriate file name (default: "NeighborTest.txt"). Otherwise it should be kept as "NULL". Output files are: 1) capsomer file; and 2) outline file. Both can be opened in Paraview. Capsomer file should be treated as "surface" (in Paraview) and Outline file should be treated as "featured edges" (in Paraview).

Particle simulation parameters such as R_e, R_c and Particle search parameters such as delr, dr, initialr, and NC should be changed by the user (if necessary) in "Tessellation.h" file. The rest of the modifications (if needed) can be implemented in "MainFunc.cpp" 
