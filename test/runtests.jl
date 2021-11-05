using NReco
using ATools
using Test
using DataFrames
using Distributions
using Statistics
using Logging
#using StatsModels
# Lower function verbosity
logger = global_logger(SimpleLogger(stdout, Logging.Warn))

fname = "../data/testdata/nema3-window-1m-LXe-20mm-1-999.h5"
pdf   = NReco.read_abc(fname)
dconf = NReco.DetConf(0.3f0, 0.05f0, 2.0f0, 100.0f0, 5000.0f0, 7)
sxyz  = pdf.sensor_xyz
wfm   = pdf.waveform
vdf   = pdf.vertices




@testset "util" begin
    @test NReco.select_by_index(sxyz,
    "sensor_id", 0) == NReco.select_by_column_value(sxyz, "sensor_id", 0)

    wevt = NReco.select_event(wfm, 4995000)
    @test mean(wevt.time) ≈ 42.39652f0
    @test NReco.select_by_column_value(wevt, "sensor_id", 1783).sensor_id[1] == 1783
    @test mean(NReco.select_by_column_value_lt(wevt, "time", 5.0).time) < 5.0
    @test mean(NReco.select_by_column_value_gt(wevt, "time", 5.0).time) > 5.0
    @test mean(NReco.select_by_column_value_interval(wevt, "time", 5.0, 10.0).time) >5.0
    @test mean(NReco.select_by_column_value_interval(wevt, "time", 5.0, 10.0).time) <10.0

    _, imx, _ = NReco.find_max_xy(wevt, "sensor_id", "time")
    m, i      = findmax(wevt.time)
    @test wevt.sensor_id[i] == imx
end

@testset "recof" begin
    @test all(NReco.primary_in_lxe(vdf).parent_id .== 0)

    df          = DataFrame(:q1=>Float32[100.0, 120.5, 132.6],
                            :q2=>Float32[123.6, 122.9,  99.9])
    filtered_df = NReco.filter_energies(df, 100.0f0, 150.0f0)
    @test nrow(filtered_df) == 1
    @test all(in(Array(filtered_df[1, :])).([120.5f0, 122.9f0]))

    df[!, :zstd1] = Float32[2.5, 3.2, 4.9]
    df[!, :zstd2] = Float32[7.9, 4.5, 2.3]
    NReco.calculate_interaction_radius!(df, x -> 2*x, "zstd")
    @test all(in(propertynames(df)).([:r1x, :r2x]))
    @test all(isapprox.(2 * df[!, :zstd1], df[!, :r1x]))
    @test all(isapprox.(2 * df[!, :zstd1], df[!, :r1x]))
end

@testset "nemareco" begin
    exp_keys = [:event_id, :phot1, :phot2, :nsipm1, :nsipm2, :q1, :q2,
	            :r1,  :r2, :r1x, :r2x,
                :phistd1, :zstd1, :widz1, :widphi1, :corrzphi1,
                :phistd2, :zstd2, :widz2, :widphi2, :corrzphi2,
			    :xs, :ys, :zs, :ux, :uy, :uz, :xt1, :yt1, :zt1,
                :t1, :xt2, :yt2, :zt2, :t2, :x1, :y1, :z1,
                :x2, :y2, :z2, :xr1, :yr1, :zr1, :tr1,
                :xr2, :yr2, :zr2, :tr2, :xb1, :yb1, :zb1, :ta1,
			    :xb2, :yb2, :zb2, :ta2]
    _, result = NReco.nemareco([fname], dconf)
    result_fields = fieldnames(typeof(result[1]))
    @test length(result_fields) == length(exp_keys)
    @test all(in(exp_keys).(result_fields))
    filter_func = ismissing, isnothing, isnan
    corrzphi1   = filter(c -> !any(f -> f(c), filter_func), getfield.(result, :corrzphi1))
    corrzphi2   = filter(c -> !any(f -> f(c), filter_func), getfield.(result, :corrzphi2))
    @test all((corrzphi1 .<= 1.0) .& (corrzphi1 .>= -1.0))
    @test all((corrzphi2 .<= 1.0) .& (corrzphi2 .>= -1.0))
end
