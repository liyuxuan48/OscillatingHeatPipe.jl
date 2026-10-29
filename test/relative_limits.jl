using Test, Random, OscillatingHeatPipe
const OHP = OscillatingHeatPipe
Random.seed!(7320)

@testset "Axial scales follow L_newbubble" begin
    @test (OHP.DEFAULT_LIQUID_MERGE_FRAC,OHP.DEFAULT_VAPOR_MERGE_FRAC,
        OHP.DEFAULT_L0_FILM_FRAC,OHP.DEFAULT_L0_PURE_VAPOR_FRAC,
        OHP.DEFAULT_BOIL_CLEARANCE_FRAC,OHP.DEFAULT_NEWBUBBLE_FILM_FRAC,
        OHP.DEFAULT_L_ADJUST_MAX_FACTOR) == (0.5,2,0.15,0.5,2,0.125,2)
    @test OHP.DEFAULT_L_NEWBUBBLE == 6e-3
    mock=(tube=(L=1.0,closedornot=true),liquid=(Xp=[(0.00075,0.00375),(0.01,0.9999)],))
    @test !any(OHP.getmerge_flags(OHP.DEFAULT_LIQUID_MERGE_FRAC*1e-3,mock))
    @test any(OHP.getmerge_flags(OHP.DEFAULT_LIQUID_MERGE_FRAC*2e-3,mock))
    @test !any(OHP.getVaporMergeFlags(OHP.DEFAULT_VAPOR_MERGE_FRAC*1e-3,mock))
    @test any(OHP.getVaporMergeFlags(OHP.DEFAULT_VAPOR_MERGE_FRAC*2e-3,mock))
end

@testset "Relative-only film formulas" begin
    for d in (1e-3,2e-4)
        @test Catoδ(d,0.) ≈ 0.002d
        @test Catoδ(d,1e6) ≈ 0.1d
        @test Catoδ(d,0.;δmin_over_r=0.006) ≈ 0.003d
        @test Catoδ(d,1e6;δmax_over_r=0.16) ≈ 0.08d
        v=OHP.Vapor(k=0.1,P=[1e5],δstart=[0.02d],δend=[0.02d],
            Lfilm_start=[1e-3],Lfilm_end=[1e-3],δmin=0.002d,
            δthreshold=0.005d,δmax=0.1d,δnucleated=0.02d)
        @test Hfilm(0.005d,v) ≈ v.k/(0.005d)+1e-6
        @test Hfilm(0.02d,v) ≈ v.k/(0.02d)
        @test Hfilm(0.002d,v)==0
    end
    @test_throws ArgumentError Catoδ(0.,0.)
    @test_throws ArgumentError Catoδ(2e-4,0.;δmin_over_r=-1.)
    @test_throws ArgumentError Catoδ(2e-4,0.;δmax_over_r=NaN)
    @test_throws MethodError Catoδ(2e-4,0.;δmin=2e-6)
    @test_throws MethodError Catoδ(2e-4,0.;film_scale=:absolute)
end

fixture=first(split(read(joinpath(@__DIR__,"thermomodel.jl"),String),"@testset"))
include_string(Main,fixture,joinpath(@__DIR__,"thermomodel.jl"))
@testset "Relative-only initialization" begin
    for d in (1e-3,2e-4)
        sys=initialize_ohpsys(sys_plate,p_fluid,power;d=d,slugnum=4)
        @test length(sys.wall.Xstations) == 250
        @test length(sys.wall.boiltime_stations) == 250
        @test issorted(sys.wall.Xstations)
        @test all(x -> 0 <= x < sys.tube.L,sys.wall.Xstations)
        @test sys.vapor.δmin ≈ 0.002d
        @test sys.vapor.δthreshold ≈ 0.005d
        @test sys.vapor.δmax ≈ 0.1d
        @test sys.vapor.δnucleated ≈ 0.02d
        @test all(isapprox.(sys.vapor.δstart,0.04*d/2))
        @test all(isapprox.(sys.vapor.δend,0.04*d/2))
    end
    d=2e-4
    make(;kwargs...)=initialize_ohpsys(sys_plate,p_fluid,power;d=d,slugnum=4,kwargs...)
    @test length(make(nucleatenum=100).wall.Xstations) == 100
    @test length(make(nucleatenum=(2*nturn+2)*20).wall.Xstations) == 680
    @test isempty(make(nucleatenum=0).wall.Xstations)
    @test_throws MethodError make(nucleate_per_channel=20)
    sys=make(δmin_over_r=0.006,δthreshold_over_r=0.012,
        δmax_over_r=0.16,δnucleated_over_r=0.06,δfilm_relative=0.06)
    @test sys.vapor.δmin ≈ 0.003d
    @test sys.vapor.δthreshold ≈ 0.006d
    @test sys.vapor.δmax ≈ 0.08d
    @test sys.vapor.δnucleated ≈ 0.03d
    @test all(isapprox.(sys.vapor.δstart,0.03d))
    L=sys.tube.L
    mid=mod(sys.liquid.Xp[1][1]+OHP.XptoLliquidslug(sys.liquid.Xp,L)[1]/2,L)
    lb=sys.wall.L_newbubble
    inserted=nucleateboiling(sys,(mid-lb/2,mid+lb/2),sys.vapor.P[1])
    @test length(inserted.vapor.P)==length(sys.vapor.P)+1
    @test all(isapprox.(inserted.vapor.δstart,0.03d))
    @test all(iszero,make(δfilm_relative=0.).vapor.δstart)
    @test_throws MethodError make(δmax=2e-5)
    @test_throws MethodError make(δfilm=4e-6)
    @test_throws MethodError make(film_scale=:absolute)
end
