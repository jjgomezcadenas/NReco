using NReco
using ATools
using Test
using Distributions
using Statistics
#using StatsModels

pdf  = NReco.read_abc("../data/testdata/nema3-window-1m-LXe-20mm-1-999.h5")
sxyz = pdf.sensor_xyz
wfm  = pdf.waveform
vdf  = pdf.vertices




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
