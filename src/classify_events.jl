using DataFrames
using ATools


function rhits(evt_charge::SubDataFrame, sensor_xyz::DataFrame,
			   ecut      ::Float32     , pde       ::Float32  )::DataFrame
	Q            = Float32.(evt_charge.charge * pde)
	sel_sipm     = evt_charge[Q .> ecut, :]
	hitdf        = leftjoin(sel_sipm, sensor_xyz, on=:sensor_id)
	disallowmissing!(hitdf, [:x, :y, :z])
	hitdf[!, :q] = Q[Q .> ecut]
	return hitdf[!, [:x, :y, :z, :q]]
end


function recoevent!(evt_counts::ATools.EventTypes, evt_charge::SubDataFrame,
					dconfig   ::DetConf, sensor_xyz::DataFrame   )

	#hit dataframe
	hitdf = rhits(evt_charge, sensor_xyz, dconfig.ecut, dconfig.pde)

	if nrow(hitdf) < 2
		evt_counts.empty += 1
		return nothing
	end

	_, _, hq1df, hq2df = lor_maxq(hitdf)

	qbound    = ATools.range_bound(dconfig.qmin, dconfig.qmax, ATools.OpenBound)
	q1_prompt = qbound(sum(hq1df.q))
	q2_prompt = qbound(sum(hq2df.q))

	if nrow(hq1df) < 2 || nrow(hq2df) < 2 ## This doesn't properly take into account both being < 2
		evt_counts.single += 1

		if nrow(hq1df) < 2
			if q2_prompt
				evt_counts.single_prompt += 1
			end
		else
			if q1_prompt
				evt_counts.single_prompt += 1
			end
		end
		return nothing
	end

	evt_counts.prompt += 1

	if q1_prompt && q2_prompt
		evt_counts.good_prompt += 1
	end
	nothing
end


function event_classifier(filenames ::Vector{String}, dconf    ::DetConf  ,
						  first_file::Integer=1     , last_file::Integer=1)

	evt_counts = ATools.EventTypes()

	ievt = 0
	for fn in filenames[first_file:last_file]
		println("reading file = ", fn)
		pdf = read_abc(fn)

		evt_counts.total += nrow(pdf.primaries)
		evt_counts.empty += length(setdiff(pdf.primaries.event_id   ,
										   pdf.total_charge.event_id))
		for evt_charge in groupby(pdf.total_charge, :event_id)
			ievt+=1
			if ievt%100 == 0
				println("reading event ", ievt, "event id =", evt_charge.event_id[1])
			end

			recoevent!(evt_counts, evt_charge, dconf, pdf.sensor_xyz)
    	end
	end
	return evt_counts
end
