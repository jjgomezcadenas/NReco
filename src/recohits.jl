using DataFrames
using Distributions
using ATools


"""
	recohits(event       ::Integer,
		 sensor_xyz  ::DataFrame,
		 waveform    ::DataFrame,
		 qcut        ::Float32,
		 pde         ::Float32,
		 sigma_tof   ::Float32)

For each event, take the sensor_xyz and waveform tables
and produce a hitdf (xyzqt) data frame, where each SiPM has information of position
(x,y,z), total charge in the SiPM and time (first photon).

"""
function recohits(event     ::Integer,
		  sensor_xyz::DataFrame,
		  waveform  ::SubDataFrame,
		  qcut      ::Float32,
		  pde       ::Float32,
		  sigma_tof ::Float32)

  	# select the waveform of this event
	if nrow(waveform) == 0
		return nothing
	end

	## make a copy of the event.
	wfm = copy(waveform)

	# Filter according to probability and PDE (pass if prob < pde)
	wfm = wfm[rand(Float32, nrow(wfm)) .< pde, :]

	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
	@debug first(wfm, 5)

	if nrow(wfm) == 0
		return nothing
	end

	# add a column of gaussian random numbers
	# representing the smearing of electronics and sensor
	tdist = Normal(0.0, sigma_tof)
	wfm[!, :dt] = Float32.(rand(tdist, nrow(wfm)))

	# add column of smeared times to true times
	transform!(wfm, [:time, :dt] => (+) => :mtime)

	@debug "waveform after prob cut and smeared time: size =$(size(wfm))"
	@debug first(wfm, 5)

	# For each sensor get tmin, tmin smeared and the number of photons.
	wfm_min = combine(groupby(wfm, :sensor_id), :time => minimum => :tmin,
		:mtime => minimum => :trmin, nrow => :q)

  	@debug "DF with time and q: size = $(size(wfm_min))"
	@debug first(wfm_min, 5)

	# cut on total charge (> qcut)
	filter!(x -> x.q .> qcut, wfm_min)

	# Get the positions of each of the selected sensors.
	xyzqt = leftjoin(wfm_min, sensor_xyz, on=:sensor_id)
	disallowmissing!(xyzqt, [:x, :y, :z])

	return xyzqt[!, [:x, :y, :z, :tmin, :trmin, :q]]
end
