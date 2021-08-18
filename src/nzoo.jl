using DataFrames
using ATools


function rhits(event::Integer, ecut::Number, pde::Number,
	Qdf::DataFrame, sxyzdf::DataFrame)

	# select the event
	qdf = select_by_column_value(Qdf, "event_id", event)
	# multiply vector of charges by PDE and add column to the DF
	Q = Float32.(qdf.charge * pde)
	qdf[!,"Q"] = Q
	# Select SiPMs with charge (q x PDE) about ecut
	qdfQ   = qdf[qdf.Q.>ecut,:]
	return sipm_xyzq(qdfQ, sxyzdf)
end


function sipm_xyzq(qdf::DataFrame, sxyz::DataFrame)
	sids = qdf.sensor_id
	pos = sipm_pos.((sxyz,),sids)
	x = [p[1] for p in pos]
	y = [p[2] for p in pos]
	z = [p[3] for p in pos]
	return DataFrame(x=x,y=y,z=z,q=qdf.Q)
end

function recoevent!(event      ::Integer,
				   dc          ::DetConf,
				   sensor_xyz  ::DataFrame,
				   waveform    ::DataFrame,
				   n3d         ::Dict)

	n3d["total"] = n3d["total"] + 1

	#hit dataframe
	hitdf = rhits(event, dc.ecut, dc.pde, waveform, sensor_xyz) 
	
	#recohits(event,sensor_xyz, waveform, dc.ecut, dc.pde, dc.sigma_tof)

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

			recoevent!(event, dconf, pdf.sensor_xyz, pdf.total_charge, n3d)
    	end
	end
	n3df = DataFrame(n3d)
	return n3df
end
