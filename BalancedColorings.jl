#=
Balanced coloring partitions for coupled cell networks

Inspired in:
Kamei, Cock - 2012 - Computation of balanced equivalence relations and their lattice for a coupled cell network
Aguiar, Dias - 2014 - The Lattice of Synchrony Subspaces of a Coupled Cell Network/ Characterization and Computation Algorithm

This code was created for replicating results from the literature.
Does not include state-of-the-art algorithms, so cannot go beyond ND=14.

This code uses Oscar to identify the symmetry group of network of each pattern

Things to be done:
Use Distributed.jl

See runtest.jl for examples.

Jaime Cisternas, Santiago, November 2025
=#

module BalancedColorings

using Oscar
using Combinatorics
using LinearAlgebra
using ProgressBars
using Graphs
using GraphMakie
using CairoMakie
const CM = CairoMakie

#=
Finds all permutations that preserve symmetry of adjacency matrix A.
=#
function findsymmetries(A::Matrix,allelments)
    println("Finding symmetries...")
    ND = size(A,1)
    elments = Vector{Vector{Int64}}()
    for e in ProgressBar(allelments)
        ii = Vector(e)
        if all(A.==A[ii,ii])
            push!(elments,e)
        end
    end
    println("Order of the group : ",length(elments))
    return elments
end

#=
Writes a csv file
 SourceName,TargetName,Type,Source,Target
 Type can be positive, negative or dual
=#
function writecsv(A,filename)
    ND = size(A,1)
	open(filename*"_net.csv","w") do file
		write(file, "Source,Target\n")
		for i in 1:ND
			for j in 1:ND
				if A[i,j]>0
					write(file,"$(i),$(j)\n")
				end
			end
		end
	end
end

#=
Finds symmetry group of graph based on adjacency matrix A.
Based on Oscar package.
=#
function find_group(A,isundirected=true)

    ND = size(A,1)
    flag = isundirected ? Undirected : Directed

	# another way to get the same group
	g = graph_from_adjacency_matrix(flag,A) # <= KEY
	geners = automorphism_group_generators(g)
	G = permutation_group(ND,geners)
	elments = elements(G) # <= so we use right order
	
	@show glabel = describe(G)
	@show small_group_identification(G)
	
	# two key Oscar methods
	sgs = subgroups(G)
	ccs = conjugacy_classes_subgroups(G)
	
	numberelements = [Int(order(s)) for s in sgs ]
	Nsgs = length(numberelements)
	Nccs = length(ccs)
	
	println()
	println("$(Nsgs) subgroups in $(Nccs) conjugacy classes")
	
	########################################
	
	# creates an empty dictionary of strings and an empty dictionary of lists
	Nelems = Int(order(G))
	sglabels = Dict{Int64, String}()
	sgsignatures = Dict{Int64, Vector{Int64}}()
	for (i,s) = enumerate(sgs)
		local mysignature
		gs = string(gens(s))
		w = findfirst('[',gs)
		gs = describe(s)*" "*gs[w:end]
		sglabels[i] = gs
		println(i," : ",sglabels[i])
		
		mysignature = Vector{Int64}()
		for (j,e) = enumerate(elments)
			bit = e in elements(s) ? 1 : 0
			push!(mysignature, bit)
		end
		#println(signature)
		sgsignatures[i] = mysignature
	end
	
	println()
	println("Binary signatures [$(Nsgs)x$(Nelems) binary matrix] :")
	for (i,s) = enumerate(sgs)
		println(sgsignatures[i])
	end
	
	classsubgroups = [findfirst(conjugacy_class(G,s)==c for c in ccs) for s in sgs]
	
	sizesccs = [0];
	for i in 1:length(ccs)
		println(i," : ",count(classsubgroups.==i))
		push!(sizesccs,count(classsubgroups.==i))
	end

	return glabel, Vector.(elments), sgs, ccs, sglabels, sgsignatures
end

allalmostequal(x) = norm(maximum(x)-minimum(x))<1e-6

#=
Finds all synchronization patterns of an adjacency matrix A.
This is the right order.
Exhaustive search, cannot go beyond ND=14
=#
function findallpatterns(A)
    println("Finding sync patterns...")

    ND = size(A,1)
	# each node in a given partition must receive the same number of arrows from the same partitions
	
	pp = []
	
	for p in ProgressBar(Combinatorics.partitions(1:ND))
	
		isadmissible = true
		ng = length(p)
		for s in p
			if length(s)>1
				conects = []
				for i in s
				    push!(conects, [sum(A[s1,i]) for s1 in p if ~(i in s1)]) # key correction
				    #push!(conects, [sum(A[i,s1]) for s1 in p]) # key correction
					#push!(conects, [sum(A[i,s1])>0 for s1 in p])
				end
				#@show conects
				isadmissible = allalmostequal(conects) & isadmissible # <=
			end
		end
	
		if ~isadmissible
			continue
		end
		
	    push!(pp,p)
	end
	println("Number of sync patterns : ",length(pp))

    return pp
end

function quotient(A,p)
    ND = size(A,1)
    ng = length(p)
	AA = zeros(Int64,ng,ng)
	for i1 in 1:ND
		for i2 in 1:ND
			if A[i1,i2]>0
				j1 = findfirst(s->(i1 in s), p)
				j2 = findfirst(s->(i2 in s), p)
				AA[j1,j2] += 1 # (j1!=j2 ? 1 : 0)
			end
		end
	end
	return AA
end

# whether two nodes belong to same partition
function areinsame(i,j,p)
    res = false
    for s in p
        if (i in s) & (i != j)
            res = j in s
        end
    end
    return res
end

# binary matrix that represents whether two nodes belong to same partition
function synchronized(p)
    ND = maximum(maximum.(p))
    S = zeros(Int64,ND,ND)
    for i in 1:ND-1
        for j in (i+1):ND
            if areinsame(i,j,p)
                S[i,j] = S[j,i] = 1
            end
        end
    end
    return S
end

#=
Finds relations between sync patterns
Builds a (redundant) binary matrix that can be used for drawing the lattice
=#
function compareall(syncpatterns,syncpatternsk)
    Nsync = length(syncpatternsk)
    C = zeros(Int64,Nsync,Nsync)
    println("Relations between sync patterns :")
    for i1 in 1:Nsync
        p1 = replace(string(syncpatternsk[i1])," "=>"")
        m1 = syncpatterns[p1].matrix
        for i2 in 1:Nsync
            p2 = replace(string(syncpatternsk[i2])," "=>"")
            m2 = syncpatterns[p2].matrix
            C[i1,i2] = Int64(all(m1.<=m2) & ~all(m1.==m2))
        end
        println(C[i1,:]," ",p1)
    end
    
    # gets rid of redundant connections, good for drawing the lattice
    Cnonred = copy(C)
    println("Nonredundant relations between sync patterns :")
    for i in Nsync:-1:2
        for j in (i-1):-1:1
            if C[i,j]==1
                Cnonred[i,:] -= C[j,:].*(Cnonred[i,:].>=C[j,:])
            end
        end
    end
    for i1 in 1:Nsync
        p1 = replace(string(syncpatternsk[i1])," "=>"")
        m1 = syncpatterns[p1].matrix
        println(Cnonred[i1,:]," ",p1)
    end
    
    return C, Cnonred
end

# builds reduce and lift matrices
function projectionmatrices(p)
    p1 = vcat(p...)
    @assert allunique(p1)
    ND = maximum(p1)
    R = zeros(Float64,length(p),ND)
    Lt = zeros(Float64,length(p),ND)
    i = 1
    for s in p
        for s1 in s
            R[i,s1] = 1.0/length(s)
            Lt[i,s1] = 1.0
        end
        i += 1
    end
    return R, Lt'
end

roll(arr, step) = vcat(arr[end-step+1:end], arr[1:end-step])

# areconjugate(classesm[54],classesm[56])

function areconjugate(m1,m2,elments)
    ND = size(m1,1)
    @assert size(m1)==size(m2)
    arec = false
    for e in elments
        m21 = m2[e,e]
        if (m1==m21)
            arec = true
            break
        end
    end
    return arec
end

#=
Classify sync patterns based on symmetry.
Finds conjugacy classes of patterns.
=#
function classify(syncpatterns,syncpatternsk,elments)
    global ND
    Nsync = length(syncpatterns)
    classes = Dict{String,Any}()
    classesk = []
    println("Classify sync patterns in symmetrically conjugated classes")
    for kk in syncpatternsk
        k = replace(string(kk)," "=>"")
        v = syncpatterns[k]
        isnew = true
        for (k1,v1) in classes
            if areconjugate(syncpatterns[k].matrix,syncpatterns[k1].matrix,elments)
                push!(classes[k1],k)
                isnew = false
                break
            end
        end
        if isnew
            classes[k] =[k]
            push!(classesk,k)
        end
    end
    
    for c in classesk
       println(c," : ",classes[c])
    end
    
    return classes, classesk
end

# generates an example pattern from a partition
function example(p)
    ND = maximum(maximum.(p))
    ex = zeros(Int64,ND)
    for (k,s) in enumerate(p)
        for i in s
            ex[i] = k
        end
    end
    return ex
end

#=
Finds best match between symmetries present in some object and the subgroups of a given group.
=#
function findbestmatch(datum::Vector{Int64},signatures::Dict{Int64,Vector{Int64}},labels::Dict{Int64,String},elems)
    @assert length(labels)==length(signatures)
    #reddatum = datum[elems]
    s = zeros(Int64,length(elems))
    for (i,e) in enumerate(elems)
        s[i] = (datum==datum[e]) ? 1 : 0
    end
    simil = zeros(Float64,length(labels))
    for (k,lab) in labels
        #simil[k] = sum(abs,[(datum[e]-signatures[lab][e]) for e in elems])
        simil[k] = sum(abs,s.-signatures[k])
    end
    #println(simil)
    kk = argmin(simil)
    return labels[kk], simil[kk]==0
end

function createdictionary(A::Matrix,syncpatternsk,elments)
	Nsync = length(syncpatternsk)
	syncpatterns = Dict{String,Any}()
	println("Details of sync patterns : ")
	for p in syncpatternsk
		local R,L
		label = replace(string(p)," "=>"")
		R,L = projectionmatrices(p)
		this = (pattern=p, matrix=synchronized(p), rank=length(p), R=copy(R), L=copy(L))
		syncpatterns[label] = this
		ex = example(p)
		fl = [ex==ex[e] for e in elments]
		println(p," ",ex," ",fl)
	end

    return syncpatterns
end

function createdictionary(A::Matrix,syncpatternsk,elments,sgsignatures,sglabels)
	Nsync = length(syncpatternsk)
	syncpatterns = Dict{String,Any}()
	println("Details of sync patterns : ")
	for p in syncpatternsk
		local R,L
		label = replace(string(p)," "=>"")
		R,L = projectionmatrices(p)
		ex = example(p)
		sy,fl = findbestmatch(ex,sgsignatures,sglabels,elments)
		this = (pattern=p, matrix=synchronized(p), symmetry=sy, rank=length(p), R=copy(R), L=copy(L))
		syncpatterns[label] = this
		println(p," ",ex," ",sy," ",fl)
	end

    return syncpatterns
end

function findclasses(syncpatterns,syncpatternsk,elments)

	# Finds classes of conjugate patterns
	classes, classesk = classify(syncpatterns,syncpatternsk,elments)
	Nclass = length(classes)
	
	invclasses = Dict{String,String}()
	for (i,kv) in enumerate(classes)
		k,v = kv
		for k2 in classes[k]
			invclasses[k2] = k
		end
	end
	
	return classes, classesk, invclasses
end

function latticeclasses(syncpatterns,syncpatternsk,classes,classesk,invclasses)

    Nclass = length(classes)
	# Finds connections between classes
	CC = zeros(Int64,Nclass,Nclass)
	for ij in findall(C.>0)
		i1,i2 = Tuple(ij)
		label1 = replace(string(syncpatternsk[i1])," "=>"")
		label2 = replace(string(syncpatternsk[i2])," "=>"")
		class1 = invclasses[label1]
		class2 = invclasses[label2]
		j1 = findfirst(isequal.(class1,classesk))
		j2 = findfirst(isequal.(class2,classesk))
		CC[j1,j2] = 1
	end
	
	ER = zeros(Int64,Nclass,Nclass)
	for i1 in 1:Nclass
		for i2 in 1:Nclass
			ER[i1,i2] = ((length(syncpatternsk[i1])==length(syncpatternsk[i2])) & (i1>i2)) ? 1 : 0
		end
	end
	
	return CC,ER
end

function plot_network(label,A,xy)
    ND = size(A,1)
    if xy == nothing
        xy = Shell()
    end
    issym = ((A.>0)==(A'.>0))
    g = issym ? SimpleGraph(A.>0) : SimpleDiGraph(A.>0)
    
    if ND<=7
        mycolors = CM.Makie.wong_colors()
    else
        mycolors = CM.Colors.HSV.(range(0, 360, ND+1), 50, 50)
    end
    
    # just the network
    
    f = CM.Figure(size=(500,500));
    ax = f[1,1] = CM.Axis(f)
    
    cl = 1:ND
    if issym
        p = graphplot!(ax, g, layout=xy, ilabels=1:ND, node_color=mycolors[cl], node_size=60)
    else
        p = graphplot!(ax, g, layout=xy, ilabels=1:ND, node_color=mycolors[cl], node_size=60, arrow_shift=:end, arrow_size=20)
    end
    #p.node_color = mycolors[cl]
    CM.hidedecorations!(ax); CM.hidespines!(ax)
    ax.aspect = DataAspect()
    CM.xlims!(ax,-1.2,1.2)
    CM.ylims!(ax,-1.2,1.2)
    
    CM.save(label*"_network.pdf",f)
end

function nicetitle(x)
    xx = replace(x[2:end-1],"["=>"\\{","]"=>"\\}")
    return L"\bowtie = %$xx"
end

function nicetitle(x,y)
    xx = replace(x[2:end-1],"["=>"\\{","]"=>"\\}")
    yy = replace(string(y),"x"=>"\\times","C2"=>"\\mathrm{C}_2","C3"=>"\\mathrm{C}_3","C4"=>"\\mathrm{C}_4","C5"=>"\\mathrm{C}_5","C6"=>"\\mathrm{C}_6",
        "D2"=>"\\mathrm{D}_1","D4"=>"\\mathrm{D}_2","D6"=>"\\mathrm{D}_3","D8"=>"\\mathrm{D}_4","D10"=>"\\mathrm{D}_5","D12"=>"\\mathrm{D}_6",
        "S2"=>"\\mathrm{S}_2","S3"=>"\\mathrm{S}_3","S4"=>"\\mathrm{S}_4","S5"=>"\\mathrm{S}_5","S6"=>"\\mathrm{S}_6")
    return L"\bowtie = %$xx ~ ; ~ \Sigma = %$yy"
end

# all sync patterns
function plot_patterns(label,A,xy,syncpatterns,syncpatternsk,classes,classesk)
    ND = size(A,1)
    if xy == nothing
        xy = Shell()
    end
    issym = ((A.>0)==(A'.>0))
    g = issym ? SimpleGraph(A.>0) : SimpleDiGraph(A.>0)
    Nsync = length(syncpatternsk)
    Nclass = length(classesk)
    
    if ND<=7
        mycolors = CM.Makie.wong_colors()
    else
        mycolors = CM.Colors.HSV.(range(0, 360, ND+1), 50, 50)
    end
    
    # plot patterns in a more or less square array
    Nc = Int64(round(sqrt(Nclass)))
    Nr = Nclass÷Nc
    Nr += (Nclass>Nr*Nc) ? 1 : 0
    f = CM.Figure(size=(500*Nc,500*Nr));
    axs = []
    i1 = 1; i2 = 1
    for i in 1:Nclass
        if i2>Nc
            i2 = 1
            i1 += 1
        end
        ax = f[i1,i2] = CM.Axis(f; title=nicetitle(classesk[i],syncpatterns[classesk[i]].symmetry))
        push!(axs,ax)
        i2 += 1
    end
    
    for (i,k) in enumerate(classesk)
        v = syncpatterns[k]
        cl = [findfirst([i in s for s in v.pattern]) for i in 1:ND]
        #f[i,1] = ax = Axis(f) # can add title=""
        if issym
            p = graphplot!(axs[i], g, layout=xy, ilabels=1:ND, node_color=mycolors[cl], node_size=60)
        else
            p = graphplot!(axs[i], g, layout=xy, ilabels=1:ND, node_color=mycolors[cl], node_size=60,arrow_shift=:end, arrow_size=20)
        end
        #p.node_color = mycolors[cl]
        CM.hidedecorations!(axs[i]); CM.hidespines!(axs[i])
        axs[i].aspect = DataAspect()
        CM.xlims!(axs[i],-1.2,1.2)
        CM.ylims!(axs[i],-1.2,1.2)
    end
    
    CM.save(label*"_patterns.pdf",f)
end

end # module