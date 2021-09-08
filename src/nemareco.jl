using DataFrames
#using Glob
using Distributions
using ATools


"""
	recohits(         event        ::Integer,
					  sensor_xyz  ::DataFrame,
					  waveform    ::DataFrame,
					  qcut        ::Float32,
					  pde         ::Float32,
					  sigma_tof   ::Float32)

For each event, take the sensor_xyz and waveform tables
and produce a hitdf (xyzqt) data frame, where each SiPM has information of position
(x,y,z), total charge in the SiPM and time (first photon).

"""
function recohits(event       ::Integer,
				  sensor_xyz  ::DataFrame,
				  waveform    ::DataFrame,
				  qcut        ::Float32,
				  pde         ::Float32,
				  sigma_tof   ::Float32)

  	# select the waveform of this event
	wfm = select_by_column_value(waveform, "event_id", event)
	if nrow(wfm) == 0
		return nothing
	end

	# add a column with probability of surviving pdf cut (pass if prob < pde)
	wfm[!, "prob"] = rand(Float32, nrow(wfm))

	# add a column of charge
	# (each photon arriving to the SiPm in the waveform has charge of one)
	wfm[!, "q"]    = Float32.(ones(nrow(wfm)))

	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
	@debug first(wfm, 5)

	# SiPM pass PDE cut if prob < PDE
	wfmc = select_by_column_value_lt(wfm, "prob", pde)
	if nrow(wfmc) == 0
		return nothing
	end

	# add a column of gaussian random numbers
	# representing the smearing of electronics and sensor
	d    = Normal(0.0, sigma_tof)
	wfmc[!,"dt"] = Float32.(rand(d, nrow(wfmc)))

	# add column of smeared times to true times
	wfmt = transform(wfmc, [:time, :dt] => (+) => :mtime)

	@debug "waveform after prob cut and smeared time: size =$(size(wfmt))"
	@debug first(wfmt, 5)

	# group by SiPMs and take minimum time
	wtmin = combine(groupby(wfmt, :sensor_id), :time => minimum)

	@debug "waveform after grouping SiPMs and min time: size = $(size(wtmin))"
	@debug first(wtmin, 5)

	# group by SiPMs and take minimum smeared time
	wrtmin = combine(groupby(wfmt, :sensor_id), :mtime => minimum)
	@debug "waveform after grouping SiPMs and min reco time: size = $(size(wrtmin))"
	@debug first(wrtmin, 5)

	# group by SiPMs and compute the sum of the charge in the SiPMs
	wtmq = combine(groupby(wfmt, :sensor_id), :q => sum)

	@debug "waveform after grouping SiPMs and sum charge: size = $(size(wtmq))"
	@debug first(wtmq, 5)

	#construct qt dataframe
	wfmx  = DataFrame(sensor_id=wtmq.sensor_id,
					  tmin=wtmin.time_minimum,
					  trmin=wrtmin.mtime_minimum,
					  q=wtmq.q_sum)

  	@debug "DF with time and q: size = $(size(wfmx))"
	@debug first(wfmx, 5)

	# cut on total charge (> qcut)
	qdft   = wfmx[wfmx.q .>qcut,:]

	#select the id's of the sipms with charge above threshold
	sids = qdft[!,:sensor_id]

	#compute positions of the SiPMs
	pos   = sipm_pos.((sensor_xyz,),sids)
	x = [p[1] for p in pos]
	y = [p[2] for p in pos]
	z = [p[3] for p in pos]

	# Construct hit data frame
	xyzqt   = DataFrame(x=x,y=y,z=z,
					    tmin=qdft.tmin,
					    trmin=qdft.trmin,
					    q=qdft.q)

	return xyzqt
end


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
	recovent(event       ::Integer,
	dc          		 ::DetConf,
	df1         		 ::DataFrame,
	df2         		 ::DataFrame,
	primaries   		 ::DataFrame,
	sensor_xyz  		 ::DataFrame,
	waveform    		 ::DataFrame,
	lor_algo    		 ::Function)

Return a dictionary with the variables characterising the event.
"""
function recovent(event      ::Integer,
				 dc          ::DetConf,
				 df1         ::DataFrame,
				 df2         ::DataFrame,
				 primaries   ::DataFrame,
				 sensor_xyz  ::DataFrame,
				 waveform    ::DataFrame,
				 lor_algo    ::Function)
	n3d = Dict()
	# primaries
	prim = select_by_column_value(primaries, "event_id", event)
	@debug "Primaries in event:" prim

		#hit dataframe
	hitdf = recohits(event,sensor_xyz, waveform, dc.ecut, dc.pde, dc.sigma_tof)

	if hitdf === nothing
	@warn "Warning, hidtf evaluates to nothing for event = $event"
	return nothing
	end

	if nrow(hitdf) < 2
	@warn "Warning, hidtf is <2 for event = $event"
	return nothing
	end

	@info " hit dataframe: size = $size(hitdf)"
	@debug first(hitdf, 5)

	# reconstruct (x,y) : barycenter
	#b1, b2, hq1df, hq2df = lor_kmeans(hitdf)
	#b1, b2, hq1df, hq2df = lor_maxq(hitdf)
	b1, b2, hq1df, hq2df = lor_algo(hitdf)

	@info " barycenters" b1 b2

	# total charge
	q1 = sum(hq1df.q)
	q2 = sum(hq2df.q)

	@info " total charge: q1 = $q1, q2 = $q2"
	if q1 < dc.qmin || q1 > dc.qmax
	@info "Warning, q1 is $q1 for event $event"
	return nothing
	end
	if q2 < dc.qmin || q2 > dc.qmax
	@info "Warning, q2 is $q2 for event $event"
	return nothing
	end

	# Compute phistd and zstd1
	phistd1 = phistd(hq1df)
	zstd1   = xyzstd(hq1df,"z")
	phistd2 = phistd(hq2df)
	zstd2   = xyzstd(hq2df,"z")
	@info " phistd1 = $phistd1, zstd1 = $zstd1"
	@info " phistd2 = $phistd2, zstd2 = $zstd2"

	# find true position (and correlate with barycenter)
	xt1, xt2 = true_xyz(b1, df1, df2)
	# find r1 and r2 (from True info)
	r1 = rxy(xt1[1], xt1[2])
	r2 = rxy(xt2[1], xt2[2])
	@info " True position in hemisphere 1" xt1
	@info " True position in hemisphere 1" xt2

	# New (x,y) positions estimated from r1, r2
	x1, y1, z1  = radial_correction(b1, r1)
	x2, y2, z2  = radial_correction(b2, r2)

	@info " New (x,y,z) positions estimated from r1, r2 & r1q, r2q"
	@info " from r1:  x1 = $x1, y1=$y1, z1=$z1"
	@info " from r2:  x2 = $x2, y1=$y2, z1=$z2"

	# Find the sipm with the fastest time
	t1  = minimum(hq1df.tmin)
	t2  = minimum(hq2df.tmin)
	tr1 = minimum(hq1df.trmin)
	tr2 = minimum(hq2df.trmin)

	ntof1 = min(dc.ntof, nrow(hq1df))
	ntof2 = min(dc.ntof, nrow(hq2df))

	# sort reco times in ascending order
	t1s = sort(hq1df.trmin)
	t2s = sort(hq2df.trmin)
	# take average
	ta1 = mean(t1s[1:ntof1])
	ta2 = mean(t2s[1:ntof2])

	@info " time"
	@info " true:  t1 = $t1, t2=$t2"
	@info " smeared:  tr1 = $tr1, tr2=$tr2"
	@info " averaged:  ta1 = $ta1, ta2=$ta2"

	ht1  = select_by_column_value(hq1df, "tmin", t1)
	ht2  = select_by_column_value(hq2df, "tmin", t2)
	htr1 = select_by_column_value(hq1df, "trmin", tr1)
	htr2 = select_by_column_value(hq2df, "trmin", tr2)

	n3d["xs"] = prim.x[1]
    n3d["ys"] = prim.y[1]
    n3d["zs"] = prim.z[1]
	n3d["ux"] =prim.vx[1]
	n3d["uy"] =prim.vy[1]
	n3d["uz"] = prim.vz[1]

    n3d["xt1"]=xt1[1]
    n3d["yt1"]=xt1[2]
    n3d["zt1"]=xt1[3]
    n3d["t1"]=t1

    n3d["xt2"]=xt2[1]
    n3d["yt2"]=xt2[2]
    n3d["zt2"]=xt2[3]
    n3d["t2"]=t2

    n3d["x1"]=x1
    n3d["y1"]=y1
    n3d["z1"]=z1

    n3d["x2"]=x2
    n3d["y2"]=y2
    n3d["z2"]=z2

    n3d["xr1"]=b1.x
    n3d["yr1"]=b1.y
    n3d["zr1"]=b1.z
    n3d["tr1"]=tr1

    n3d["xr2"]=b2.x
    n3d["yr2"]=b2.y
    n3d["zr2"]=b2.z
    n3d["tr2"]=tr2

    n3d["ta1"]=ta1
    n3d["ta2"]=ta2

	n3d["xb1"]=ht1.x[1]
    n3d["yb1"]=ht1.y[1]
    n3d["zb1"]=ht1.z[1]
    n3d["xb2"]=ht2.x[1]
    n3d["yb2"]=ht2.y[1]
    n3d["zb2"]=ht2.z[1]

    n3d["nsipm1"]=nrow(hq1df)
    n3d["q1"]= sum(hq1df.q)
    n3d["r1"]=r1
    n3d["phistd1"]=phistd1
    n3d["zstd1"]=zstd1

    n3d["nsipm2"]=nrow(hq2df)
    n3d["q2"]= sum(hq2df.q)
    n3d["r2"]=r2
    n3d["phistd2"]=phistd2
    n3d["zstd2"]=zstd2

	return hq1df, hq2df, n3d
end


"""
	nema_analysis!(     event       ::Integer,
						dc          ::DetConf,
						df1         ::DataFrame,
						df2         ::DataFrame,
						primaries   ::DataFrame,
						sensor_xyz  ::DataFrame,
						waveform    ::DataFrame,
						lor_algo    ::Function,
						n3d         ::Dict)

Fill the DataFrame for nema analysis
"""
function nema_dict!(event       ::Integer,
						dc          ::DetConf,
						df1         ::DataFrame,
						df2         ::DataFrame,
						primaries   ::DataFrame,
						sensor_xyz  ::DataFrame,
						waveform    ::DataFrame,
						lor_algo    ::Function,
						n3d         ::Dict)


	result = recovent(event, dc, df1, df2,primaries, sensor_xyz, waveform, lor_algo)

	if result !== nothing
		hitdf1, hitdf2, evtd = result
		ks =  keys(evtd)
		#ks2 =  keys(n3d)
		#println(ks)
		#println(ks2)
		#@assert all(ks .== ks2) == true
		push!(n3d["phot1"], df1.process_id[1] == 1)
		push!(n3d["phot2"], df2.process_id[1] == 1)

		for k in ks
			push!(n3d[k],evtd[k])
		end
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

	# define data dictionary

	n3d = Dict("phot1"=>Bool[], "phot2"=>Bool[],
			   "nsipm1"=>Int64[],"nsipm2"=>Int64[],
			   "q1" =>Float32[],   "q2" =>Float32[],
	           "r1"  =>Float32[],  "r2"  =>Float32[],
			   "phistd1"=>Float32[],  "zstd1"=>Float32[],
			   "phistd2"=>Float32[],  "zstd2"=>Float32[],
			   "xs"=>Float32[], "ys"=>Float32[], "zs"=>Float32[],
		       "ux"=>Float32[], "uy"=>Float32[], "uz"=>Float32[],
	           "xt1"=>Float32[], "yt1"=>Float32[], "zt1"=>Float32[],"t1"=>Float32[],
		       "xt2"=>Float32[], "yt2"=>Float32[], "zt2"=>Float32[], "t2"=>Float32[],
               "x1"=>Float32[],   "y1"=>Float32[], "z1"=>Float32[],
               "x2"=>Float32[],   "y2"=>Float32[], "z2"=>Float32[],
			   "xr1"=>Float32[], "yr1"=>Float32[], "zr1"=>Float32[], "tr1"=>Float32[],
               "xr2"=>Float32[], "yr2"=>Float32[], "zr2"=>Float32[], "tr2"=>Float32[],
			   "xb1"=>Float32[], "yb1"=>Float32[], "zb1"=>Float32[], "ta1"=>Float32[],
			   "xb2"=>Float32[], "yb2"=>Float32[], "zb2"=>Float32[], "ta2"=>Float32[])

	# read one file to compute the radius of sipm
	#pdf = read_abc(files[1])
	#rsipm = rxy(pdf.sensor_xyz.x[2], pdf.sensor_xyz.y[2])

	for file in files[file_i:file_l]               # loop on files
		println("reading file = ", file)
		pdf = read_abc(file)            # read file
		dfs = primary_in_lxe(pdf.vertices)       # primary photons in LXe

		for event in unique(dfs.event_id)       #loop on events
			#  event DF
			vdf = select_by_column_value(dfs, "event_id", event)

			# two primary photons in LXe
			if any(vdf.track_id .== 1) && any(vdf.track_id .== 2)
				df1 = select_by_column_value(vdf, "track_id", 1)
            	df2 = select_by_column_value(vdf, "track_id", 2)

				nema_dict!(event, dconf, df1, df2,
					       pdf.primaries, pdf.sensor_xyz, pdf.waveform,
						   lor_algo, n3d)
			end
    	end
	end
	n3df = DataFrame(n3d)
end
