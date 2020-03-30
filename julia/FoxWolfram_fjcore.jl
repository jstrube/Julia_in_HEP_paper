using Pkg
Pkg.activate(".")
using CxxWrap
using LinearAlgebra
using Plots
using LCIO
using BenchmarkTools
using EventShapes
using Printf
using StaticArrays
using FastJet
import Base: getindex

# to help with the threading we'll copy the data from LCIO, so the event loop can move on to the next event while we process
function getCollectionCopy(evt, name)
    col = getCollection(evt, name)
    p4vec = Vector{SVector{4, Float64}}(undef, length(col))
    @inbounds for (idx, particle) in enumerate(col)
        p3 = getMomentum(particle)
        p4vec[idx] = SVector(p3[1], p3[2], p3[3], sqrt(sum(p3.^2)))
    end
    return p4vec
end

function main()
    # read file names from command line (only argument) 
    if length(ARGS) < 1 
        println("You need to specify at least one input file")
        return 1
    end
    R = 0.7
    jet_def = JetDefinition(antikt_algorithm, R)
    h1 = Float64[]
    h2 = Float64[]
    h3 = Float64[]
    h4 = Float64[]
    hThrust = Float64[]
	sum = 0.0
    @inbounds for nfiles=1:length(ARGS)
        FILEN = ARGS[nfiles]
        println("Opening File ", nfiles, " ", FILEN)
        LCIO.open(FILEN) do reader
            for evt in reader
    			particles = PseudoJet[]
                for particle in getCollection(evt, "PandoraPFOs")
                    (E, p3) = getP4(particle)
    				push!(particles, PseudoJet(p3[1], p3[2], p3[3], E))
                end
				if length(particles) < 2 continue end
    			# run the clustering, extract the jets
   	 			cs = ClusterSequence(StdVector(particles), jet_def)
    			jets = inclusive_jets(cs, 1.0)
                for j in jets
                    #@printf("%.8f %.8f %.8f %.8f\n", px(j), py(j), pz(j), e(j))
                    sum += px(j) + py(j) + pz(j) + e(j)
                end
                col = getCollectionCopy(evt, "PandoraPFOs")
                h10, h20, h30, h40 = EventShapes.foxWolframMoments(col)
                thrust = EventShapes.Thrust(col)
				# @printf("%.8f %.8f %.8f %.8f %.8f\n", h10, h20, h30, h40, thrust)
                push!(h1, h10)
                push!(h2, h20)
                push!(h3, h30)
                push!(h4, h40)
                push!(hThrust, thrust)
            end
        end
    end
    #histogram(h1)
    #savefig("h1.pdf")
    #histogram(h2)
    #savefig("h2.pdf")
    #histogram(h3)
    #savefig("h3.pdf")
    #histogram(h4)
    #savefig("h4.pdf")
    #histogram(hThrust)
    #savefig("thrust.pdf")
end

@time main()
