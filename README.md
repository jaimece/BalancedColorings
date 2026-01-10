# BalancedColorings
Finds symmetries and all balanced coloring partitions (also known as equitable partitions) for coupled cell networks.

Its basic assumptions are: identical nodes coupled in an additive way.

Inspired in:

* Stewart, Golubitsky, Pivato - 2003 - Symmetry groupoids and patterns of synchrony in coupled cell networks

* Kamei, Cock - 2012 - Computation of balanced equivalence relations and their lattice for a coupled cell network

* Aguiar, Dias - 2014 - The Lattice of Synchrony Subspaces of a Coupled Cell Network/ Characterization and Computation Algorithm

and several others.

This code was created for replicating different results found the literature from 2003 - 2025.
Exhaustivity was the goal.
Does not include state-of-the-art algorithms or sophisticated accelerations, so it cannot go beyond N = 14 oscillators.

Its main input is the (possibly weighted) adjacency matrix of the network. Self-loops are OK.
We use the convention A[source, target].

The main outputs are the symmetry group and all the possible synchronization patterns.
See the examples for additional functionalities.

It works with networks that are:

* Regular/nonregular (the vertical sums of A are/aren't equal)

* Symmetric/nonsymmetric (with respect to permutations of rows and columns)

* Undirected/directed (symmetric/unsymmetric matrix A)

* Unweighted/weighted (binary/real valued A)

* Nondiffusive/diffusive (self-loops and intra-cluster interactions are ignored for external equitable partitions)

We provide two versions of the module:

(1) BalancedColoringsLite.jl

This is the lightweight version, does not require Oscar.jl.
Without Oscar we cannot identify the symmetry group of the network and of each pattern.
Also the subgroups cannot be identified.
See runtestlite.jl for examples.

(2) BalancedColorings.jl

This is the regular version, uses Oscar.
See runtest.jl for examples.

Requirements:

* Oscar.jl

* Graphs.jl

* CairoMakie.jl

* GraphsMakie.jl

Future plans:

* Improve use of keyword arguments, namespaces, etc.

* Improve lattice plot

* Include Distributed.jl

* Integrate with Graphs.jl

* Integrate with NetworkDynamics.jl

* Allow for different types of nodes (for instance in multiple layers)

* Allow for different types of arrows
