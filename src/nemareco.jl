using DataFrames
using Distributions
using StatsBase
using ATools


"""
	true_xyz(b1::Hit, b2::Hit, df1::DataFrame, df2::DataFrame)

Return the sorted distance between baricenter and true
"""
function true_xyz(b1::Hit, df1::DataFrame, df2::DataFrame)
	# distance between baricenter and true
	d1 = dxyz([b1.x, b1.y, b1.z], [df1.x[1], df1.y[1], df1.z[1]])
	d2 = dxyz([b1.x, b1.y, b1.z], [df2.x[1], df2.y[1], df2.z[1]])

	if d2 < d1
		xt2 = [df1.x[1], df1.y[1], df1.z[1]]
		xt1 = [df2.x[1], df2.y[1], df2.z[1]]
	else
		xt1 = [df1.x[1], df1.y[1], df1.z[1]]
		xt2 = [df2.x[1], df2.y[1], df2.z[1]]
	end
	return xt1, xt2
end


"""
	check_valid_hits(hitdf::DataFrame)

Check the validity of the event hitmaps for continued processing
"""
function check_valid_hits(hitdf::Union{DataFrame, Nothing}, event::Int64)
	if hitdf === nothing
		@warn "Warning, hidtf evaluates to nothing for event = $event"
		return false
	end

	if nrow(hitdf) < 2
		@warn "Warning, hidtf is <2 for event = $event"
		return false
	end
	return true
end


"""
	check_valid_charge(h1::DataFrame, h2::DataFrame, dc::DetConf)

Check if the hemispheres have detected charge above threshold
"""
function check_valid_charge(h1::DataFrame, h2::DataFrame, dc::DetConf)
	q1 = sum(h1.q)
	q2 = sum(h2.q)
	@info " total charge: q1 = $(n3d[:q1]), q2 = $(n3d[:q2])"
	if q1 < dc.qmin || q1 > dc.qmax
		@info "Warning, q1 is $(n3d[:q1]) for event $event"
		return false
	end

	if q2 < dc.qmin || q2 > dc.qmax
		@info "Warning, q2 is $(n3d[:q2]) for event $event"
		return false
	end
	return true
end


"""
	recovent(event       ::Integer,
	dc          		 ::DetConf,
	df1         		 ::DataFrame,
	df2         		 ::DataFrame,
	primaries   		 ::SubDataFrame,
	sensor_xyz  		 ::DataFrame,
	waveform    		 ::SubDataFrame,
	lor_algo    		 ::Function)

Return a dictionary with the variables characterising the event.
"""
function recovent(event     ::Integer,
		  dc        ::DetConf,
		  df1       ::DataFrame,
		  df2       ::DataFrame,
		  primaries ::SubDataFrame,
		  sensor_xyz::DataFrame,
		  waveform  ::SubDataFrame,
		  lor_algo  ::Function)
	# hit dataframe
	hitdf = recohits(event,sensor_xyz, waveform, dc.ecut, dc.pde, dc.sigma_tof)

	if !check_valid_hits(hitdf, event)
		return nothing
	end

	# reconstruct (x,y) : barycenter
	b1, b2, hq1df, hq2df = lor_algo(hitdf)
	@info " barycenters" b1 b2

	if !check_valid_charge(hq1df, hq2df, dc)
		return nothing
	end

	# find true position (and correlate with barycenter)
	xt1, xt2 = true_xyz(b1, df1, df2)
	@info " True position in hemisphere 1" xt1
	@info " True position in hemisphere 2" xt2

	## Weights and phi positions
	int1_weights = FrequencyWeights(hq1df.q)
	phi1_values  = fphi(hq1df)
	int2_weights = FrequencyWeights(hq2df.q)
	phi2_values  = fphi(hq2df)

	true_rad1 = rxy(xt1[1], xt1[2])
	true_rad2 = rxy(xt2[1], xt2[2])
	# New (x,y) positions estimated from r1, r2
	xyz1 = radial_correction(b1, true_rad1)
	xyz2 = radial_correction(b2, true_rad2)

	ntof1 = min(dc.ntof, nrow(hq1df))
	ntof2 = min(dc.ntof, nrow(hq2df))

	# sort reco times in ascending order
	t1s = sort(hq1df.trmin)
	t2s = sort(hq2df.trmin)

	tmin1 = minimum(hq1df.tmin)
	tmin2 = minimum(hq2df.tmin)
	ht1   = select_by_column_value(hq1df, "tmin", tmin1)
	ht2   = select_by_column_value(hq2df, "tmin", tmin2)

	@info " New (x,y,z) positions estimated from r1, r2 & r1q, r2q"
	@info " from r1:  x1 = $(n3d[:x1]), y1=$(n3d[:y1]), z1=$(n3d[:z1])"
	@info " from r2:  x2 = $(n3d[:x2]), y1=$(n3d[:y2]), z1=$(n3d[:z2])"

	@info " hit dataframe: size = $size(hitdf)"
	@debug first(hitdf, 5)

	reduced_event = ATools.EventParameters(
	event_id=event,
	phot1 = df1.process_id[1] == 1,
	phot2 = df2.process_id[1] == 1,
	# Primary particle origins
	xs = primaries.x[1],
	ys = primaries.y[1],
	zs = primaries.z[1],
	ux = primaries.vx[1],
	uy = primaries.vy[1],
	uz = primaries.vz[1],
	# barycenter on SiPM plane
	xr1    = b1.x,
	yr1    = b1.y,
	zr1    = b1.z,
	nsipm1 = nrow(hq1df),
	xr2    = b2.x,
	yr2    = b2.y,
	zr2    = b2.z,
	nsipm2 = nrow(hq1df),
	# total charge
	q1 = sum(hq1df.q),
	q2 = sum(hq2df.q),
	# True interaction positions
	xt1 = xt1[1],
	yt1 = xt1[2],
	zt1 = xt1[3],
	t1  = tmin1,
	tr1 = minimum(hq1df.trmin),
	xt2 = xt2[1],
	yt2 = xt2[2],
	zt2 = xt2[3],
	t2  = tmin2,
	tr2 = minimum(hq2df.trmin),
	# find r1 and r2 (from True info)
	r1 = true_rad1,
	r2 = true_rad2,
	# Compute phistd and zstd of detected light
	phistd1  = std(phi1_values, int1_weights, corrected=true),
	zstd1     = std(hq1df.z   , int1_weights, corrected=true),
	widz1     = maximum(hq1df.z) - minimum(hq1df.z),
	widphi1   = maximum(phi1_values) - minimum(phi1_values),
	corrzphi1 = cor(hcat(hq1df.z, phi1_values), int1_weights)[1,2],
	phistd2   = std(phi2_values, int2_weights, corrected=true),
	zstd2     = std(hq2df.z   , int2_weights, corrected=true),
	widz2     = maximum(hq2df.z) - minimum(hq2df.z),
	widphi2   = maximum(phi2_values) - minimum(phi2_values),
	corrzphi2 = cor(hcat(hq2df.z, phi2_values), int2_weights)[1,2],
	# XYZ given true radius and barycentre
	x1 = xyz1[1],
	y1 = xyz1[2],
	z1 = xyz1[3],
	x2 = xyz2[1],
	y2 = xyz2[2],
	z2 = xyz2[3],
	# Average arrival time of first ntof SiPMs
	ta1 = mean(t1s[1:ntof1]),
	ta2 = mean(t2s[1:ntof2]),
	# Positions of SiPM with tmin
	xb1 = ht1.x[1],
	yb1 = ht1.y[1],
	zb1 = ht1.z[1],
	xb2 = ht2.x[1],
	yb2 = ht2.y[1],
	zb2 = ht2.z[1])

	return hq1df, hq2df, reduced_event
end


"""
	nema_analysis!(     data_vec    ::Vector{ATools.EventParameters},
						event       ::Integer,
						dc          ::DetConf,
						df1         ::DataFrame,
						df2         ::DataFrame,
						primaries   ::SubDataFrame,
						sensor_xyz  ::DataFrame,
						waveform    ::SubDataFrame,
						lor_algo    ::Function)

Fill the DataFrame for nema analysis
"""
function nema_dict!(data_vec      ::Vector{ATools.EventParameters}   ,
		    event     ::Integer     ,
		    dc        ::DetConf     ,
		    df1       ::DataFrame   ,
		    df2       ::DataFrame   ,
		    primaries ::SubDataFrame,
		    sensor_xyz::DataFrame   ,
		    waveform  ::SubDataFrame,
		    lor_algo  ::Function    )

	result = recovent(event, dc, df1, df2, primaries, sensor_xyz, waveform, lor_algo)

	if result !== nothing
		push!(data_vec, result[3])
	end
end


"""
	function nemareco(files    ::Vector{String},
					  dconf    ::DetConf,
		              file_i   ::Integer=1,
					  file_l   ::Integer=1,
					  phot     ::Bool=true
					  lor_algo ::Function=lor_maxq)

Return the evtdf DataFrame


   nsipm1  => Number of sipms in hemisphere 1 (nsipm2 for hemisphere 2)
   q1      => Total charge in hemisphere
   r1      => True radius, gamma interaction point,
   phistd  =>  Phi standard deviations (weigthed by charge)
   zstd    => Z standard deviations (weigthed by charge)

   xs, ys, zs    => Position of the pair of gammas (in target),
   ux, uy, uz    => Direction vectors of gammas

   xt1, yt1, zt1 => True position of gammas,
   x1,   y1, z1  => Reco position of gammas (from barycenter and r1),
   xr1, yr1, zr  => Baricenter
   xb1, yb1, zb1 => Position of the sipm giving time stamp

   t1            => Time stamp of first photon,
   tr1           => Time stamp of first photon, after pdf and smearing
   ta1           => Average of time stamps, first five photons

"""
function nemareco(files    ::Vector{String},
				  dconf    ::DetConf,
	              file_i   ::Integer=1,
				  file_l   ::Integer=1,
				  lor_algo ::Function=lor_maxq)

	output_vector = ATools.EventParameters[]
	total_events  = 0

	for file in files[file_i:file_l]               # loop on files
		println("reading file = ", file)
		pdf = read_abc(file)            # read file

		# count the total number of events generated.
		# This, like many other sections assumes primary generation of gammas.
		total_events += nrow(pdf.primaries)

		dfs = primary_in_lxe(pdf.vertices)       # primary photons in LXe

		## We are interested in events with two primary photons in LXe
		grp_dfs = filter(x -> any(x.track_id .== 1) && any(x.track_id .== 2),
							groupby(dfs, :event_id))
		primaries = groupby(pdf.primaries, :event_id)
		waveforms = groupby(pdf.waveform , :event_id)

		for (event, vdf) in pairs(grp_dfs)
			df1 = vdf[vdf.track_id .== 1, :]
			df2 = vdf[vdf.track_id .== 2, :]

			try
				wvf = waveforms[values(event)]
				nema_dict!(output_vector, event.event_id, dconf, df1, df2,
					   primaries[values(event)], pdf.sensor_xyz,
					   wvf, lor_algo)
			catch err
				if err isa KeyError
					continue
				else
					rethrow(err)
				end
			end
    	end
	end
	return total_events, output_vector
end
