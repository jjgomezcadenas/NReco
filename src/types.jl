using Configurations


"""
	Hit

Represents a hit, characterized by position and charge
"""
struct Hit
	x::Float32
	y::Float32
	z::Float32
	q::Float32
end


"""
	DetConf

Defines detector configuration
"""
@option struct DetConf
	pde      ::Float32 =    0.3
	sigma_tof::Float32 =    0.05
	ecut     ::Float32 =    2.0
	qmin     ::Float32 =  100.0
	qmax     ::Float32 = 5000.0
	ntof     ::Integer =    7
	nsigma   ::Integer =    3
end


"""
	CalConfig

Configuration variables for DOI/radius of interaction calibration.
"""
@option struct CalFunction
	cal_file::String  = "../config/radius_calibration.h5"
	cal_grp ::String  = "n3-40mm"
	cal_std ::String  = "cstd"
	nLXe    ::Float32 = 1.69
end

@option struct CalConfig
  input_dir::String
  conf_dir ::String
  plot_dir ::String
  qmin     ::Float32
  qmax     ::Float32
  save_cal ::Bool = false
  cal_func ::CalFunction = CalFunction()
end
