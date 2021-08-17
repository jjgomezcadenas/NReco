using DataFrames
using ATools

function recoevent!(event      ::Integer,
				   dc         ::DetConf,
				   sensor_xyz  ::DataFrame,
				   waveform    ::DataFrame,
				   n3d         ::Dict)

	n3d["total"] = n3d["total"] + 1

	#hit dataframe
	hitdf = recohits(event,sensor_xyz, waveform, dc.ecut, dc.pde, dc.sigma_tof)

	if hitdf === nothing
		n3d["empty"] = n3d["empty"] + 1
		return n3d
	end

	if nrow(hitdf) < 2
		n3d["empty"] = n3d["empty"] + 1
		return n3d
	end

	_, _, hq1df, hq2df = lor_maxq(hitdf)

	if nrow(hq1df) < 2 || nrow(hq2df) < 2
		n3d["single"] = n3d["single"] + 1

		if nrow(hq1df) < 2 
			q = sum(hq2df.q)
			if q > dc.qmin && q < dc.qmax
				n3d["single-prompt"] = n3d["single-prompt"] + 1
			end

		else
			q = sum(hq1df.q)
			if q > dc.qmin && q < dc.qmax
				n3d["single-prompt"] = n3d["single-prompt"] + 1
			end
		end

		return n3d
	end


	n3d["prompt"]=n3d["prompt"] + 1

	q1 = sum(hq1df.q)
	q2 = sum(hq2df.q)

	if q1 > dc.qmin && q1 < dc.qmax
		if q2 > dc.qmin && q2 < dc.qmax
			n3d["good-prompt"]=n3d["good-prompt"] + 1
		end
	end

	return n3d 
end


function zoo(files    ::Vector{String},
			 dconf    ::DetConf,
	         file_i   ::Integer=1,
			 file_l   ::Integer=1)

	# define data dictionary

	n3d = Dict("total" =>0, "empty"=>0,"single"=>0, "prompt"=>0, 
	           "single-prompt"=>0, "good-prompt"=>0)

	ievt = 0
	for file in files[file_i:file_l]               # loop on files
		println("reading file = ", file)
		pdf = read_abc(file)            # read file
		for event in unique(pdf.vertices.event_id)       #loop on events
			ievt+=1 
			if ievt%100 == 0
				println("reading event ", ievt, "event id =", event)
			end
			
			recoevent!(event, dconf, pdf.sensor_xyz, pdf.waveform, n3d)
    	end
	end
	n3df = DataFrame(n3d)
	return n3df
end
