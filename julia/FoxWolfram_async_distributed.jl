using Pkg
cd(@__DIR__)
Pkg.activate("./")
Pkg.instantiate()

using BenchmarkTools
using Distributed
@everywhere using CxxWrap
@everywhere using LCIO
@everywhere using StaticArrays
using FastJet
using LinearAlgebra
using EventShapes

# Unfortunately, LCIO is not multi-threaded
# However, we can just start a number of independent LCIO processes, each of which reads one file and fills an event channel
# Then we can use multi-threading to process the events in the main process
# but using multi-threading 

# to help with the threading we'll copy the data from LCIO, so the event loop can move on to the next event while we process
@everywhere function readEvents(fnames, events, done)
	iEvents = 0
	while true
		try
			fn = take!(fnames)
			LCIO.open(fn) do reader
				for evt in reader
					col = getCollection(evt, "PandoraPFOs")
					collection = Vector{SVector{4, Float64}}(undef, length(col))
					@inbounds for (idx, particle) in enumerate(col)
						p3 = getMomentum(particle)
						collection[idx] = SVector(p3[1], p3[2], p3[3], sqrt(sum(p3.^2)))
					end
					put!(events, collection)
					iEvents += 1
				end
			end
		catch e
			break
		end
	end		
	put!(done, iEvents)
end

function processEvents(events, nProcessed)
	iEvents = 0
	R = 0.7
    jet_def = JetDefinition(antikt_algorithm, R)
	while true
		try
			collection = take!(events)
			h10, h20, h30, h40 = EventShapes.foxWolframMoments(collection)
			thrust = EventShapes.Thrust(collection)
			nJets = 0.0
			particles = PseudoJet[]
			for p in collection
				push!(particles, PseudoJet(p[1], p[2], p[3], p[4]))
			end
			if length(particles) >= 2 
				# run the clustering, extract the jets
				cs = ClusterSequence(StdVector(particles), jet_def)
				jets = inclusive_jets(cs, 1.0)
				nJets = 1.0*length(jets)
			end
			# we'll just use a simple representation of the results
			# h1, h2, h3, h4, thrust, nJets
			results = (h10, h20, h30, h40, thrust, nJets)
			iEvents += 1
		catch e
			put!(nProcessed, iEvents)
			break
		end
	end
end

function main()
	# let's start with reading a number of files concurrently
	fnames = RemoteChannel(()->Channel{String}(400))

	# let's make a buffer large enough for up to 2000 events concurrently
	events = RemoteChannel(()->Channel(2000))

	# the readers can signal when they are done reading events
	done = RemoteChannel(()->Channel{Int}(400))

	# for signalling the end of processing 
	nProcessed = Channel{Int}(1000)

	# spawn the readers, one per worker
	readers = [@spawnat w readEvents(fnames, events, done) for w in workers()]
	println("workers spawned", readers)
	# spawn the processors, one per thread on this process
    processors = [Threads.@spawn processEvents(events, nProcessed) for w in 2:Threads.nthreads()]
	println("processors spawned ", processors)

	# prepare the file names
	for f in ARGS
		put!(fnames, f)
	end
	close(fnames)
	println("fnames distributed")

	# wait for all readers to be done
	# then we can close the event queue and the doers can finish
	nDone = 0
	nEvents = 0
	while nDone != nworkers()
		nEvents += take!(done)
		nDone += 1
	end
	close(events)

	# wait for the processing to finish
	totalProcessed = 0
	for p in processors
		wait(p)
		totalProcessed += take!(nProcessed)
	end
	close(nProcessed) # kind of unnecessary at this point, but let's clean up
	println(totalProcessed, " events processed")
end

@time main()
