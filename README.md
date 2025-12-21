# BalancedColorings
Balanced coloring partitions for coupled cell networks

Inspired in:

Kamei, Cock - 2012 - Computation of balanced equivalence relations and their lattice for a coupled cell network

Aguiar, Dias - 2014 - The Lattice of Synchrony Subspaces of a Coupled Cell Network/ Characterization and Computation Algorithm

This code was created for replicating results from the literature.

Does not include state-of-the-art algorithms, so cannot go beyond ND=14.

BalancedColoringsLite.jl
This is the lightweight version does not require Oscar.jl.
Without Oscar we cannot identify the symmetry group of network of each pattern.

BalancedColorings.jl
Pro version, uses Oscar.

Things to be done:

Use Distributed.jl

See runtest.jl for examples.
