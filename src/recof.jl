using DataFrames
using StatsModels
using Clustering
using Statistics
using LinearAlgebra
using ATools
using HDF5

# Selection

"""
    primary_in_lxe(verticesdf::DataFrame)

Select primary photons in LXe.

Row │ event_id track_id parent_id x y z t moved pre_KE post_KE deposited process_id  volume_id
LXe           : volume_id = 0
Primary       : parent_id = 0
"""
function primary_in_lxe(verticesdf::DataFrame)
    vlxe      = select_by_column_value(verticesdf, "volume_id", 0)
    vlxephepr = select_by_column_value(vlxe, "parent_id", 0)
    return vlxephepr
end


"""
	xyz_dot(hitdf::DataFrame, simax::Hit)

Return the dot product between each SiPM in the event and the SiPM of max charge
"""
function xyz_dot(hitdf::DataFrame, simax::Hit)
	xyzmax = [simax.x, 	simax.y, simax.z]
	xyzm_dot = dot(xyzmax, xyzmax)
	return [dot(Array(hitdf[i,1:3]), xyzmax) /xyzm_dot  for i in 1:nrow(hitdf)]
end


"""
	sipm_pos(sxyz::DataFrame, index::Integer)

Return the position of the SiPMs in the sipm_xyz database for a given index

"""
function sipm_pos(sxyz::DataFrame, index::Integer)
	return Array(select_by_index(sxyz, "sensor_id", index)[1,2:end])
end


"""
	find_xyz_sipm_qmax(hitdf::DataFrame)
Return the coordinates of the SiPM with maximum charge

"""
function find_xyz_sipm_qmax(hitdf::DataFrame)
	indx, xmax, qxmax = find_max_xy(hitdf, "x", "q")
	return Hit(xmax, hitdf.y[indx], hitdf.z[indx], qxmax)
end


"""
	sipmsel(hdf::DataFrame)

Return two data frames, separating the SiPMs in two groups depending on
the sign of the angle with the SiPM of max charge
"""
function sipmsel(hdf::DataFrame)
	simax = find_xyz_sipm_qmax(hdf)
	npr   = xyz_dot(hdf, simax)
	return hdf[(npr.>0), :], hdf[(npr.<0), :]
end


"""
	ksipmsel(hdf::DataFrame, ka::Vector{Int64})

Return two data frames, separating the SiPMs in two groups depending on
the value of vector ka (1 or 2)
"""
function ksipmsel(hdf::DataFrame, ka::Vector{Int64})
	return hdf[(ka.==2), :], hdf[(ka.==1), :]
end


# kmeans algorithm
"""
	get_hits_as_matrix(hitdf::DataFrame)

Given the dataframe hitdf (with fields x,y,z), return the
underlying matrix
"""
function get_hits_as_matrix(hitdf::DataFrame)
	f = @formula(0 ~ x + y + z)
	f2 = apply_schema(f, schema(f, hitdf))
	resp, pred = modelcols(f2, hitdf)
	return transpose(pred)
end

# lors

"""
	lor_maxq(hitdf::DataFrame)

Compute lors using the SiPM of max charge (maxq algo) to divide the event
in two hemispheres and estimating the lor vertex from barycenter
"""
function lor_maxq(hitdf::DataFrame)
	hq2df, hq1df = sipmsel(hitdf)   # select using maxq
	b1 = baricenter(hq1df)          # baricenters
	b2 = baricenter(hq2df)
	return b1, b2, hq1df, hq2df
end


"""
	lor_kmeans(hitdf::DataFrame)

Compute lors using kmeans clustering.
Returns two estimations of vertices. One based in pure kmeans, the
other in barycenter.
"""
function lor_kmeans(hitdf::DataFrame)
	Mhits = get_hits_as_matrix(hitdf)  # take the underlying matrix
	kr = kmeans(Mhits, 2)              # apply kmeans
	ka = assignments(kr) # get the assignments of points to clusters

	hq2df, hq1df = ksipmsel(hitdf, ka)   # select using kmeans list
	b1 = baricenter(hq1df)     # baricenters
	b2 = baricenter(hq2df)
	return b1, b2, hq1df, hq2df
end


"""
	baricenter(hdf::DataFrame)
	returns the barycenter of a cluster of hits
"""
function baricenter(hdf::DataFrame)
	function xq(hdf::DataFrame, pos::String)
		return sum(hdf[!,pos] .* hdf.q) / qt
	end
	qt = sum(hdf.q)
	return Hit(xq(hdf, "x"), xq(hdf, "y"), xq(hdf, "z"), qt)
end


"""
	radial_correction(b::Hit, r::Float32, rsipm::Float32)

Take the estimated radius of interaction (r), and the radius of the sipm
and return the corrected positions
"""
function radial_correction(b::Hit, r::Float32)

	ϕ = atan(b.y,b.x)
    x = r .* cos.(ϕ)
    y = r .* sin.(ϕ)
    return x,y,b.z
end


"""
	phistd(hitdf::DataFrame)

Compute the std deviation in phi weighted by charge, e.g:
Sqrt(1/Q Sum_i (phi_i - phi_mean) * qi )
"""
function phistd(hitdf::DataFrame)
	phi = fphi(hitdf)
	return wstd(phi, hitdf.q)
end


"""
	xyzstd(hitdf::DataFrame)

Compute the std deviation in x weighted by charge, e.g:
Sqrt(1/Q Sum_i (phi_i - phi_mean) * qi )
"""
function xyzstd(hitdf::DataFrame, column::String="x")
	x = hitdf[!, column]
	return wstd(x, hitdf.q)
end


"""
    filter_energies
"""
function filter_energies(df::DataFrame, qmin::Float32, qmax::Float32)
    interval = ATools.range_bound(qmin, qmax, ATools.OpenBound)
    ndfq     = filter(x -> interval.(x.q1) .& interval.(x.q2), df)
end


"""
    calibration_function
"""
function calibration_function(calibFunc::NReco.CalFunction, rmin::Real, rmax::Real)
    line_pars = h5open(calibFunc.cal_file) do h5cal
        bias   = read_attribute(h5cal[calibFunc.cal_grp], calibFunc.cal_std * "-bias" )
        lconst = read_attribute(h5cal[calibFunc.cal_grp], calibFunc.cal_std * "-const")
        llin   = read_attribute(h5cal[calibFunc.cal_grp], calibFunc.cal_std * "-lin"  )
        return bias, lconst, llin
    end
    cal_func = ATools.predict_interaction_radius(ATools.gpol1(line_pars[2:end]),
        rmin, rmax, line_pars[1])
    return cal_func
end


"""
    calculate_interaction_radius!

"""
function calculate_interaction_radius!(df       ::DataFrame,
                                       predictor::Function,
                                       variable ::String,
                                       rmax     ::Union{Real, Nothing}=nothing)
    if variable == "cstd" && !(:cstd in propertynames(df))
        # We want to use the combined std which isn't necessarily saved in the H5
        comb_std(z, phi) = sqrt(z^2 + (rmax * phi)^2)
        transform!(df, [:zstd1, :phistd1] => ByRow(comb_std) => :cstd1,
            [:zstd2, :phistd2] => ByRow(comb_std) => :cstd2)
    end
    sym1 = Symbol(variable * "1")
    sym2 = Symbol(variable * "2")
    transform!(df, sym1 => predictor => :r1x, sym2 => predictor => :r2x)
    nothing
end
