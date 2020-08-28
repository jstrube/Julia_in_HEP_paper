using Pkg
cd(@__DIR__)
Pkg.activate("./")
Pkg.instantiate()

using CxxWrap
using LCIO
using OnlineStats
using Distributions
using BenchmarkTools

function invMass(p1, p2)
    E = p1[1] + p2[1]
    p = p1[2] .+ p2[2]
    sqrt(E^2 - sum(p.^2))
end


function main(filelist)
    hist = Hist(85:1:95)
    normal = FitNormal()
    s = Series(hist, normal)
    AllEvents = 0
    for nfiles = 1:length(filelist)
        nEvt = 0
         LCIO.open(filelist[nfiles]) do reader
            println(nfiles," ",filelist[nfiles]," ",length(reader))
            for event in reader
                nEvt += 1
                muon1 = (0, 0)
                muon2 = (0, 0)
                for particle in getCollection(event, "PandoraPFOs")
                    if abs(getType(particle)) != 13
                        continue
                    end
                    m = getP4(particle)
                # we'll keep the two muons with the highest energies
                    if m[1] > muon1[1]
                        muon1 = m
                    elseif m[1] > muon2[1]
                        muon2 = m
                    end
                end
            # if we found two muons, we'll compute the invariant mass
                if muon1[1] > 0 && muon2[1] > 0
                    i = invMass(muon1, muon2)
                    if i < 85 || i > 95
                        continue
                    end
                # keep the inv mass in the histogram and in a likelihood fit
                    fit!(s, i)
                end
            end
        end
        AllEvents += nEvt
    end
    println("Events : ", AllEvents)
    println(normal)
    
end

@time main(ARGS)
