# BalancedColorings
Finds all balanced coloring partitions (also known as equitable partitions) for coupled cell networks

Its basic assumptions are: identical nodes coupled in an additive way.

Inspired in:

Stewart, Golubitsky, Pivato - 2003 - Symmetry groupoids and patterns of synchrony in coupled cell networks
Kamei, Cock - 2012 - Computation of balanced equivalence relations and their lattice for a coupled cell network
Aguiar, Dias - 2014 - The Lattice of Synchrony Subspaces of a Coupled Cell Network/ Characterization and Computation Algorithm

and several others.

This code was created for replicating different results found the literature from 2003 - 2025.

Does not include state-of-the-art algorithms or sophisticated accelerations, so it cannot go beyond N = 14 oscillators.

Its main input is the (weighted) adjacency matrix of the network.

The output are the symmetry group and all the possible synchronization patterns.

It works with networks that are:

Regular/nonregular

Symmetry/nonsymmetric

Directed/undirected

Unweighted/weighted

Nondiffusive/diffusive

We provide two versions of the module:

BalancedColoringsLite.jl

This is the lightweight version does not require Oscar.jl.
Without Oscar we cannot identify the symmetry group of the network and of each pattern.
Also the subgroups cannot be identified.
See runtest.jl for examples.

BalancedColorings.jl

This is the regular version, uses Oscar.

Future plans:

Improve lattice plot

Include Distributed.jl

Integrate with Graphs.jl

Integrate with NetworkDynamics.jl
Allow for different nodes (for instance in multiple layers)
