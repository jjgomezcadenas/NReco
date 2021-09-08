using NReco
using ATools
using Test
using Distributions
using Statistics
#using StatsModels

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

    mx, imx = NReco.find_max_xy(wevt, "sensor_id", "time")
    m, i    = findmax(wevt.time)
    @test wevt.sensor_id[i] == imx
end

@testset "recof" begin
    @test all(NReco.primary_in_lxe(vdf).parent_id .== 0) == 1
end

@testset "nemareco" begin
    exp_keys = ["phot1", "phot2", "nsipm1", "nsipm2", "q1",   "q2",
	            "r1",  "r2", "phistd1", "zstd1", "corrzphi1", "phistd2",  "zstd2", "corrzphi2",
			    "xs", "ys", "zs", "ux", "uy", "uz", "xt1", "yt1", "zt1",
                "t1", "xt2", "yt2", "zt2", "t2", "x1", "y1", "z1",
                "x2", "y2", "z2", "xr1", "yr1", "zr1", "tr1",
                "xr2", "yr2", "zr2", "tr2", "xb1", "yb1", "zb1", "ta1",
			    "xb2", "yb2", "zb2", "ta2"]
    result = NReco.nemareco([fname], dconf)
    @test length(names(result)) == length(exp_keys)
    @test all(in(exp_keys).(names(result)))
end
