if abspath(PROGRAM_FILE) == @__FILE__
	using Pkg
	Pkg.activate(normpath(joinpath(@__DIR__, "..")))
end

using DataFrames
using Configurations
using HDF5
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
	detconf = args["config"]

	if isdir(outd) == false
		mkdir(outd)
	end

	output_path = joinpath(outd, outf)
	filenames   = glob("*.h5",dr)

	if detconf != "default"
		dconf = from_toml(NReco.DetConf, detconf)
	else
		dconf = NReco.DetConf()
	end

	println("number of files in data dir = $(length(filenames))")
	println("reading = $(file_l - file_i + 1) files")
	println("output file  = $output_path")

	h5open(output_path, "w") do h5out
		grp        = create_group(h5out, "EventCounts")
		dtype      = ATools.generate_hdf5_datatype(ATools.EventTypes)

		evt_counts = NReco.event_classifier(filenames, dconf, file_i, file_l)

		dspace     = dataspace([evt_counts])
		dset       = create_dataset(grp, "counts", dtype, dspace)
		write_dataset(dset, dtype, [evt_counts])
	end
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
		"--config", "-c"
			help = "detector configuration"
			default = "default"
			arg_type = String
    end

    return parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
	parsed_args = parse_commandline()
	println("Running makezoo with arguments", parsed_args)
	makezoo(parsed_args)
end
