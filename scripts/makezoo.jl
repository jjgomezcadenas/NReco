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

function makezoo(args)
	dr      = args["dir"]
	outd    = args["odir"]
	outf    = args["ofile"]
	file_i  = args["filei"]
	file_l  = args["filel"]
	
	if isdir(outd) == false
		mkdir(outd)
	end

	output = string(outd,"/", outf)
	files = glob("*.h5",dr)

	output = string(outd,"/", outf)
	files = glob("*.h5",dr)

    pde  = 0.3f0
    sigma_tof = 0.050f0
    ecut = 2.0f0           # best performance for small  ecut
    qmin = 1900.0f0
    qmax = 2200.0f0
    ntof =7                 # increase sipms for average
	
	dconf = NReco.DetConf(pde, sigma_tof, ecut, qmin, qmax, ntof)

	println("number of files in data dir = $(length(files))")
	println("reading = $(file_l - file_i + 1) files")
	println("output file  = $output")

	n3df = NReco.zoo(files, dconf, file_i, file_l)
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
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makezoo with arguments", parsed_args)
	makezoo(parsed_args)
end

@time main()
