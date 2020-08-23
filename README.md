# Turing Patterning in Discrete Heterogeneous Media

Codes we implemented to generate the results presented in my dissertation called: "Localisation of Turing Patterning in Discrete Heterogeneous Media", at the University of Oxford. Below there is a brief description of the functions that have been implemented.

makeGraph: returns the adjacency matrix of a graph with "N" nodes and topology "type" and the coordinates of its nodes. It also plots the resulting graph and computes the results using Schnakenberg kinetics with given parameters.

make2Graphs: returns the adjacency matrix of the full system composed by two graphs with "N1" and "N2" nodes and topologies "type1" and "type2", respectively. Kt also returns the coordinates of its nodes, plots the resulting graph and computes the results using Schnakenberg kinetics.

solve1Graph: solves the system of differential equations determined by a reaction diffusion system with Schnakenberg kinetics, with adjacecy matrix "A1" and kinetic parameters "param".

solve2Graphs: solves the system of differential equations determined by a reaction diffusion system with Schnakenberg kinetics on two weakly connected graphs. 

wgPlot2: plots a graph given its adjacency matrix and vertices coordinates.

ComputeGraphAsymptotics: returns the asymptotic approximations of the eigenvalues of the full system composed by two coupled graphs. Scalar case.

ComputeGraphAsymptoticsSystem: returns the asymptotic approximations of the eigenvalues of the full system composed by two coupled graphs. System case.

In order to be able to generate random graphs the package "matlab_bgl" is needed. Moreover, to generate circulant graphs we used an already implemented function called "circulant" which we also include here. 

circulant: generates a square circulant matrix.
