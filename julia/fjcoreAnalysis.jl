using Pkg
cd(@__DIR__)
Pkg.activate("./")
Pkg.instantiate()

using CxxWrap
using LCIO
using OnlineStats
using Distributions
using BenchmarkTools
using FastJet
using Printf

function invMass(p1, p2)
    E = p1[1] + p2[1]
    p = p1[2] .+ p2[2]
    sqrt(E^2 - sum(p.^2))
end


function main(filelist)
    # choose a jet definition
    R = 0.7
    jet_def = JetDefinition(antikt_algorithm, R)
    AllEvents = 0
	sum = 0.0
    for nfiles = 1:length(filelist)
        nEvt = 0
         LCIO.open(filelist[nfiles]) do reader
            println(nfiles," ",filelist[nfiles]," ",length(reader))
            for event in reader
                nEvt += 1
    			particles = PseudoJet[]
                for particle in getCollection(event, "PandoraPFOs")
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
    			#println(length(particles), " ", length(jets))
            end
        end
        AllEvents += nEvt
    end
    println("Events : ", AllEvents)
	@printf("%.8f\n", sum)
end

@time main(ARGS)
