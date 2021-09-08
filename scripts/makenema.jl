#push!(LOAD_PATH,"../src/")
#using DrWatson
#@quickactivate(@__DIR__)
#@quickactivate "JPetalo"
#quickactivate("../.", "JPetalo")
#include("../src/JPetalo.jl")
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using DataFrames
using CSV
using Glob
using ArgParse
using Logging
using Printf
using ATools
using NReco

logger = SimpleLogger(stdout, Logging.Warn)
old_logger = global_logger(logger)


function makenema(args)

	lorf    = args["loralgo"]
	detconf = args["detconf"]
	dr      = args["dir"]
	outd    = args["odir"]
	outf    = args["ofile"]
	file_i  = args["filei"]
	file_l  = args["filel"]
	qxmin   = Float32(args["qmin"])
	qxmax   = Float32(args["qmax"])

	lor_algo = NReco.lor_kmeans
	if lorf == "lor_qmax"
		lor_algo = NReco.lor_maxq
	end

	if isdir(outd) == false
		mkdir(outd)
	end
	output = string(outd,"/", outf)
	files = glob("*.h5",dr)

	if detconf == "pde_1_sigmatof_1ps"
		pde  = 1.0f0
		sigma_tof = 0.001f0
		ecut = 10.0f0
		qmin = 100.0f0
		qmax = 50000.0f0
		ntof =5
		dconf = NReco.DetConf(pde, sigma_tof, ecut, qmin, qmax,  ntof)
	elseif detconf == "pde_0.3_sigmatof_1ps"
		pde  = 0.3f0
		sigma_tof = 0.001f0
		ecut = 3.0f0
		qmin = 100.0f0
		qmax = 5000.0f0
		ntof =5
	elseif detconf == "pde_0.3_sigmatof_85ps"
		pde  = 0.3f0
		sigma_tof = 0.085f0
		ecut = 3.0f0
		qmin = 100.0f0
		qmax = 5000.0f0
		ntof =5
	elseif detconf == "pde_0.3_sigmatof_50ps"
		pde  = 0.3f0
		sigma_tof = 0.05f0
		ecut = 3.0f0
		qmin = 100.0f0
		qmax = 5000.0f0
		ntof =5
	else                # by default, 40 ps jitter  & 30 ps electronics = 50 ps, harder cuts
		pde  = 0.3f0
		sigma_tof = 0.050f0
		ecut = 2.0f0           # best performance for small  ecut
		qmin = 100.0f0
		qmax = 5000.0f0
		ntof =7                 # increase sipms for average
	end

	if qxmin > 0
		qmin = qxmin
	end

	if qxmax > 0
		qmax = qxmax
	end

	dconf = NReco.DetConf(pde, sigma_tof, ecut, qmin, qmax, ntof)

	println("makenema configuration")
	println("detector configuration", dconf)
	println(" lor_algo  = $lorf")

	println("number of files in data dir = $(length(files))")
	println("reading = $(file_l - file_i + 1) files")
	println("output file  = $output")

	n3df = NReco.nemareco(files, dconf, file_i, file_l, lor_algo)

	CSV.write(output, n3df)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dir", "-d"
            help = "directory with nema simulations"
            arg_type = String
		"--odir", "-o"
            help = "output directory"
            arg_type = String
		"--ofile", "-x"
            help = "output file"
            arg_type = String
            default = "evtdf.csv"
		"--filei", "-i"
	        help = "number of initial file in glob list"
	        default  = 1
			arg_type = Int
		"--filel", "-l"
		    help = "number of last file in glob list"
		    default  = 1
			arg_type = Int
		"--loralgo", "-g"
			help = "algorithm to use for LOR reconstruction "
			arg_type = String
			default = "lor_kmeans"
		"--detconf", "-c"
			help = "Detector configuration"
			arg_type = String
			default = "all"
		"--qmin", "-q"
			help = "minimum charge (negative to set by detconf)"
			arg_type = Float64
			default = -1.0
		"--qmax", "-Q"
			help = "maximum charge (negative to set by detconf)"
			arg_type = Float64
			default = -1.0
		#"--phot", "-p"
		#	help = "Select photoelectric if 1"
		#	action = :store_true
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makenema with arguments", parsed_args)
	makenema(parsed_args)
end

@time main()
