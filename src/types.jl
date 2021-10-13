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
end
