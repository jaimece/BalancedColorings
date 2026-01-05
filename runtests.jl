#=
Test some of the functionalities in BalancedColorings.jl
using some well known examples from the literature on coupled cell networks

Let's assume Oscar is available.
Without Oscar we could not identify the symmetry group of network or of patterns.

Jaime Cisternas,Santiago, December 2025
=#

if ~("./" in LOAD_PATH)
    push!(LOAD_PATH,"./")
end

using Revise
using Combinatorics
using LinearAlgebra
using BalancedColorings
const BC = BalancedColorings

#####

# Pivato network, also known as bowtie
# directed, regular, nonsymmetric

label = "pivato"
println("\n$(label)\n")

# builds adjacency matrix, uses A[source,target]
ND = 5
A = zeros(Int64,ND,ND);
A[1,4] = 1
A[1,5] = 1
A[2,1] = 1
A[2,5] = 1
A[3,2] = 1
A[4,3] = 1
A[5,1] = 1
A[5,2] = 1
A[5,3] = 1
A[5,4] = 1

sumA = sum(A,dims=1)
@assert allequal(sumA)

xy = [(1,1), (-1,1), (-1,-1), (1,-1), (0,0)]

glabel, elments, sgs, ccs, sglabels, sgsignatures = BC.find_group(A,false)
Nelments = length(elments)

sumA = sum(A,dims=2)

minsyncpattern = BC.findminimalpattern(A)

syncpatternsk = BC.findallpatterns(A)
Nsync = length(syncpatternsk)

println("Some examples")
@show p = syncpatternsk[5]
@show R,L = BC.projectionmatrices(p)

# generate state in polydiagonal space
@show x = L*randn(length(p))

# measures distance to polydiagonal space
@show norm(x - L*R*x)

syncpatterns = BC.createdictionary(A,syncpatternsk,elments,sgsignatures,sglabels)

C, Cnonred = BC.compareall(syncpatterns,syncpatternsk);

# Finds classes of conjugate patterns
classes, classesk, invclasses = BC.findclasses(syncpatterns,syncpatternsk,elments)
Nclass = length(classes)
println("Number of sync patterns : $(Nsync), number of classes : $(Nclass)")

BC.plot_network(label,A,xy)

BC.plot_patterns(label,A,xy,syncpatterns,syncpatternsk,classes,classesk)

CC,ER = BC.latticeclasses(Cnonred,syncpatterns,syncpatternsk,classes,classesk,invclasses)

BC.plot_lattice(label,CC,ER,syncpatterns,classesk)

#####

# Ring of 6 nodes with next-nearest neighbors
# undirected, regular, symmetric

label = "nicol"
println("\n$(label)\n")

# builds adjacency matrix, uses A[source,target]
ND = 6
NN = 2
A = zeros(Int64,ND,ND);
for i in 1:ND
    for j in 1:ND
        if (min(abs(i-j+ND),abs(i-j-ND),abs(i-j))<=NN) & (i!=j)
            A[i,j] = 1
        end
    end
end

sumA = sum(A,dims=1)
@assert allequal(sumA)

xy = [(cos(2*pi*t/ND),sin(2*pi*t/ND)) for t in 0:(ND-1)]

glabel, elments, sgs, ccs, sglabels, sgsignatures = BC.find_group(A)
Nelments = length(elments)

sumA = sum(A,dims=2)

minsyncpattern = BC.findminimalpattern(A)

syncpatternsk = BC.findallpatterns(A)
Nsync = length(syncpatternsk)

syncpatterns = BC.createdictionary(A,syncpatternsk,elments,sgsignatures,sglabels)

C, Cnonred = BC.compareall(syncpatterns,syncpatternsk);

# Finds classes of conjugate patterns
classes, classesk, invclasses = BC.findclasses(syncpatterns,syncpatternsk,elments)
Nclass = length(classes)
println("Number of sync patterns : $(Nsync), number of classes : $(Nclass)")

BC.plot_network(label,A,xy)

BC.plot_patterns(label,A,xy,syncpatterns,syncpatternsk,classes,classesk)

CC,ER = BC.latticeclasses(Cnonred,syncpatterns,syncpatternsk,classes,classesk,invclasses)

BC.plot_lattice(label,CC,ER,syncpatterns,classesk)

#####

# Lodi
# directed, nonregular, nonsymmetric

label = "lodi"
println("\n$(label)\n")

# builds adjacency matrix, uses A[source,target]
ND = 5
A = zeros(Int64,ND,ND);
A[2,1] = 1
A[1,2] = 1
A[1,3] = 1
A[5,3] = 1
A[4,3] = 1
A[2,4] = 1
A[3,4] = 1
A[5,4] = 1
A[1,5] = 1
A[3,5] = 1
A[4,5] = 1

sumA = sum(A,dims=1)

xy = [(-1,1), (1,1), (-1,-1), (1,-1), (0,0)]

glabel, elments, sgs, ccs, sglabels, sgsignatures = BC.find_group(A,false)
Nelments = length(elments)

sumA = sum(A,dims=2)

minsyncpattern = BC.findminimalpattern(A)

syncpatternsk = BC.findallpatterns(A)
Nsync = length(syncpatternsk)

syncpatterns = BC.createdictionary(A,syncpatternsk,elments,sgsignatures,sglabels)

C, Cnonred = BC.compareall(syncpatterns,syncpatternsk);

# Finds classes of conjugate patterns
classes, classesk, invclasses = BC.findclasses(syncpatterns,syncpatternsk,elments)
Nclass = length(classes)
println("Number of sync patterns : $(Nsync), number of classes : $(Nclass)")

BC.plot_network(label,A,xy)

BC.plot_patterns(label,A,xy,syncpatterns,syncpatternsk,classes,classesk)

CC,ER = BC.latticeclasses(Cnonred,syncpatterns,syncpatternsk,classes,classesk,invclasses)

BC.plot_lattice(label,CC,ER,syncpatterns,classesk)

#####

# Sorrentino
# undirected, nonregular, symmetric

label = "sorrentino"
println("\n$(label)\n")

# builds adjacency matrix, uses A[source,target]
ND = 5
A = zeros(Int64,ND,ND);
A[1,2] = A[2,1] = 1
A[2,3] = A[3,2] = 1
A[3,4] = A[4,3] = 1
A[4,1] = A[1,4] = 1
A[1,5] = A[5,1] = 1
A[2,5] = A[5,2] = 1
A[3,5] = A[5,3] = 1
A[4,5] = A[5,4] = 1

sumA = sum(A,dims=1)

xy = [(-1,1),(1,1),(1,-1),(-1,-1),(0,0)]

glabel, elments, sgs, ccs, sglabels, sgsignatures = BC.find_group(A)
Nelments = length(elments)

sumA = sum(A,dims=2)

minsyncpattern = BC.findminimalpattern(A)

syncpatternsk = BC.findallpatterns(A)
Nsync = length(syncpatternsk)

syncpatterns = BC.createdictionary(A,syncpatternsk,elments,sgsignatures,sglabels)

C, Cnonred = BC.compareall(syncpatterns,syncpatternsk);

# Finds classes of conjugate patterns
classes, classesk, invclasses = BC.findclasses(syncpatterns,syncpatternsk,elments)
Nclass = length(classes)
println("Number of sync patterns : $(Nsync), number of classes : $(Nclass)")

BC.plot_network(label,A,xy)

BC.plot_patterns(label,A,xy,syncpatterns,syncpatternsk,classes,classesk)

CC,ER = BC.latticeclasses(Cnonred,syncpatterns,syncpatternsk,classes,classesk,invclasses)

BC.plot_lattice(label,CC,ER,syncpatterns,classesk)

#####

# Siddique
# undirected, nonregular, symmetric

label = "siddique"
println("\n$(label)\n")

# builds adjacency matrix, uses A[source,target]
ND = 10
A = [0 0 0 0 0 1 1 0 0 1;
0 0 0 0 0 1 1 0 0 1;
0 0 0 0 1 0 0 1 1 0;
0 0 0 0 1 0 0 1 1 0;
0 0 1 1 0 1 0 0 1 0;
1 1 0 0 1 0 0 0 0 1;
1 1 0 0 0 0 0 1 0 1;
0 0 1 1 0 0 1 0 1 0;
0 0 1 1 1 0 0 1 0 0;
1 1 0 0 0 1 1 0 0 0]

sumA = sum(A,dims=1)

xy = [(-0.5, -1),(0.5, -1),(-0.5, 1),(0.5, 1),(-1, 0.5),(-1, -0.5),(1, -0.5),(1, 0.5),(0, 0.5),(0, -0.5)]

glabel, elments, sgs, ccs, sglabels, sgsignatures = BC.find_group(A)
Nelments = length(elments)

sumA = sum(A,dims=2)

minsyncpattern = BC.findminimalpattern(A)

syncpatternsk = BC.findallpatterns(A)
Nsync = length(syncpatternsk)

syncpatterns = BC.createdictionary(A,syncpatternsk,elments,sgsignatures,sglabels)

C, Cnonred = BC.compareall(syncpatterns,syncpatternsk);

# Finds classes of conjugate patterns
classes, classesk, invclasses = BC.findclasses(syncpatterns,syncpatternsk,elments)
Nclass = length(classes)
println("Number of sync patterns : $(Nsync), number of classes : $(Nclass)")

BC.plot_network(label,A,xy)

BC.plot_patterns(label,A,xy,syncpatterns,syncpatternsk,classes,classesk)

CC,ER = BC.latticeclasses(Cnonred,syncpatterns,syncpatternsk,classes,classesk,invclasses)

BC.plot_lattice(label,CC,ER,syncpatterns,classesk)

#####

