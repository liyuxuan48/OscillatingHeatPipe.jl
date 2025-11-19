using Statistics

@testset "Liquid and vapor arrays" begin

    L, Lmin, ϕ0 = 4.0 + rand(), 0.001, 0.6 + 0.1*rand()
    for closedornot in [true,false]
        X, dXdt, liquid_realratio = randomXp(L,Lmin,closedornot;chargeratio=ϕ0)

        # if this is a closed case, then shift it so that the last
        # slug crosses the line. Ensures more robust testing
        xshift = closedornot ? L - mean(X[end]) : 0.0
        X .= map(u -> (mod(u[1]+xshift,L),mod(u[2]+xshift,L)),X)        

        lslug = mod.(map(u -> u[2],X) - map(u -> u[1],X),L)
        @test all(lslug .> 0)

        lslug2 = OscillatingHeatPipe.XptoLliquidslug(X,L)
        @test lslug == lslug2

        ϕ = sum(lslug)/L
        @test abs(ϕ - ϕ0) < 0.1
        @test ϕ ≈ liquid_realratio

        Xvapor = OscillatingHeatPipe.getXpvapor(X,closedornot)
        lvap = mod.(map(u -> u[2],Xvapor) - map(u -> u[1],Xvapor),L)
        @test all(lvap .> 0)

        lvap2 = OscillatingHeatPipe.XptoLvaporplug(X,L,closedornot)
        
        @test lvap == lvap2

        # Continuity test
        @test sum(lvap) + sum(lslug) ≈ L

        # Test checking positions in liquid slugs.
        np = length(X)

        # Liquid point should return true
        ir = rand(1:np-1)
        @test OscillatingHeatPipe.ifamong(mean(X[ir]),X)

        # Vapor point should return false
        ir = rand(2:np)
        Xi, Xf = X[ir-1][2], X[ir][1]
        @test !OscillatingHeatPipe.ifamong(0.5*(Xi+Xf),X)

        # Test creation of liquid slug arrays
        N = rand(200:300)
        θi = rand()
        Xarray, θarray = OscillatingHeatPipe.constructXarrays(X,N,θi,L)
        # check that all points lie between 0 and L
        @test all(map(u -> all(0 .<= u .<= L),Xarray))
        # both arrays for a slug are the same length
        @test all(map((u,v) -> length(u)==length(v),Xarray,θarray))


        ## test assembly and disassembly 
        
        dXdt .= [(rand(),rand()) for j in 1:np]
        u = OscillatingHeatPipe.Xtovec(X,dXdt)

        ir = rand(1:np)
        Xi, Xf = X[ir]
        dXi, dXf = dXdt[ir]
        @test u[2*ir] == Xf && u[2*ir-1] == Xi
        @test u[2*np + 2*ir] == dXf && u[2*np + 2*ir - 1] == dXi

        M = rand(np)
        δstart = rand(np)
        δend = rand(np)
        Lfilm_start = rand(np)
        Lfilm_end = rand(np)
        
        # Test that assembly and disassembly are inverses of each other
        u = OscillatingHeatPipe.XMδLtovec(X,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end)

        X2,dXdt2,M2,δstart2,δend2,Lfilm_start2,Lfilm_end2 = OscillatingHeatPipe.vectoXMδL(u)

        @test X2 == X
        @test dXdt2 == dXdt
        @test M2 == M
        @test δstart2 == δstart
        @test δend2 == δend
        @test Lfilm_start2 == Lfilm_start
        @test Lfilm_end2 == Lfilm_end

        # test a case in which there is one more X and dXdt element than the others
        # This occurs in an open tube, e.g., vapor plugs at ends and one slug between them
        # so that Xp contains (0,Xp1) and (Xpn,L) to describe the end vapor plugs.
        X = [(rand(),rand()) for j in 1:np+1]
        dXdt = [(rand(),rand()) for j in 1:np+1]
        u = OscillatingHeatPipe.XMδLtovec(X,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end)
        X2,dXdt2,M2,δstart2,δend2,Lfilm_start2,Lfilm_end2 = OscillatingHeatPipe.vectoXMδL(u)

        @test X2 == X
        @test dXdt2 == dXdt
        @test M2 == M
        @test δstart2 == δstart
        @test δend2 == δend
        @test Lfilm_start2 == Lfilm_start
        @test Lfilm_end2 == Lfilm_end

        

    end

end

ρₛ = 2730 # material density [kg/m^3]
cₛ  = 8.93e02 # material specific heat [J/kg K]
kₛ  = 1.93e02 # material heat conductivity
αₛ = kₛ/ρₛ/cₛ
    
dₛ = 1.5e-3

Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref) # This function relies on CoolProp.jl package

power = 70 # [W], total power
areaheater_area = 50e-3 * 50e-3 # [m] total area

function get_qbplus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbplus = zeros_surface(base_cache)
    return qbplus
end
    
function get_qbminus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbminus = zeros_surface(base_cache)
    # qbminus .= nrm.u
    return qbminus
end

bcdict = Dict("exterior" => get_qbplus,"interior" => get_qbminus)
    

phys_params = Dict( "diffusivity"              => αₛ,
                    "flux_correction"          => ρₛ*cₛ*dₛ,
                    # "angular velocity"         => 0.0,
                    "Fourier"                  => 1.0,
                    "ohp_flux"                 => [NaN], # initial value, the value here is useless
                    "areaheater_power"         => power, # total power
                    "areaheater_area"          => areaheater_area, # total area
                    "areaheater_temp"          => 0.0,   # relative temperature compared with "background temperature"
                    "areaheater_coeff"         => 4000.0,
                    "background temperature"   => Tref
                     )

Δx = 0.0007 # [m] # grid size, at the same order of 1D OHP channel node spacing ~ 0.001[m]

Lx = 6*INCHES*1.02 # plate size x [m]
Ly = 2*INCHES*1.05 # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits
    
g = PhysicalGrid(1.03 .* xlim,1.1 .* ylim,Δx)

Δs = 1.4*cellsize(g) # 1D OHP node spacing, here it is 1.4Δx

xbound = [ -Lx/2,-Lx/2, 
             Lx/2, Lx/2] # x coordinates of the shape

ybound = [  Ly/2,-Ly/2, 
            -Ly/2, Ly/2] # y coordinates of the shape

body = Polygon(xbound,ybound,Δs)
    
X = MotionTransform([0,0],0) # move the plate or rotate the plate
joint = Joint(X)
m = RigidBodyMotion(joint,body)
x = zero_motion_state(body,m)
update_body!(body,x,m)

function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["flux_correction"] 
end

function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["areaheater_coeff"]
    corr = phys_params["flux_correction"] 
    
    σ .= h*(T0 - T) / corr
end

fregion1_h = Rectangle(25e-3,25e-3,1.4*Δx)
tr1_h = RigidTransform((0.0,-0.0),0.0)
heater1 = AreaForcingModel(fregion1_h,tr1_h,heatermodel!)

fregion1_c = Rectangle(15e-3,1.0INCHES,1.4*Δx)
tr1_c = RigidTransform((2.4INCHES,-0.0),0.0)
cond1 = AreaForcingModel(fregion1_c,tr1_c,condensermodel!)

ds = 1.5Δx
nturn = 13
width_ohp = 46.25*1e-3
length_ohp = 147.0*1e-3
gap = 3e-3
pitch = width_ohp/(2*nturn+1)
x0, y0 = length_ohp/2 +2e-3, width_ohp/2

x, y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,3pi/2)
ohp = BasicBody(x,y) # build a BasicBody based on x,y
tr_ohp = RigidTransform((0.0,0.0),0.0)

function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
    σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"] 
end
ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!)

forcing_dict = Dict("heating models" => [heater1,cond1,ohp_linesource])
    
timestep_fixed(u,sys) = tstep

prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             # timestep_func=timestep_fourier
                                             timestep_func=timestep_fixed)

sys_plate = construct_system(prob)



@testset "Plate" begin

    ohp = OscillatingHeatPipe._get_ohp_from_forcing_list(sys_plate)
    
    @test ohp == ohp_linesource

                                                 
end


# @testset "Plate ADI solver" begin

#     sys_plate_test = deepcopy(sys_plate)
#     u_plate = newstate(sys_plate_test) .+ Tref # initialize plate T field to uniform Tref
#     integrator_plate = init(u_plate,tspan,sys_plate_test) # construct integrator_plate

#     ohp = sys_plate.qline
    
#     @test ohp[1].body == ohpgeom.body

#     #	Single eigenmode decay

#     xylim = sys_plate_test.grid.xlim
#     Δx = sys_plate_test.grid.Δx
#     I0 = sys_plate_test.grid.I0

#     xgrid = zero(temperature(integrator_plate))
#     ygrid = zero(xgrid)
#     for i in eachindex(xgrid[:,1])
#         xgrid[i,:] .= xylim[1][1] .+ Δx*(i-1)
#     end
#     for i in eachindex(ygrid[1,:])
#         ygrid[:,i] .= xylim[2][1] .+ Δx*(i-1)
#     end

#     temperature(integrator_plate) .= cos.(π .* xgrid ./ (xylim[1][2] .+ 0.5Δx)) .* cos.(2π .* ygrid ./ (xylim[2][2] .+ 0.5Δx))
#     α =integrator_plate.p.params.α
#     ρₛ = integrator_plate.p.params.ρ
#     cₛ = integrator_plate.p.params.c

#     T_hist = []
#     Tgrid = temperature(integrator_plate);
#     Tgrid_ini = deepcopy(Tgrid)
#     tstep = 5e-3
#     tend = 1.0
#     ts = tstep:tstep:tend

#     for i in ts
#         push!(T_hist,deepcopy(ADI_newT!(Tgrid,sys_plate_test,tstep)))
#     end


#     j = rand(1:size(T_hist,1))
#     Tratio_max = maximum(T_hist[j] ./ Tgrid_ini)
#     Tratio_min = minimum(T_hist[j] ./ Tgrid_ini)

#     # verify they all decay at the same rate
#     @test isapprox(Tratio_max,Tratio_min,rtol=1e-10)

#     Lx = xylim[1][2] .+ 0.5Δx
#     Ly = xylim[2][2] .+ 0.5Δx
    
#     # then verify maximum magnitude decay exponentially
#     @test isapprox(log.(maximum.(T_hist) ./ maximum(Tgrid_ini)),-α * ((π/Lx)^2 + (2π/Ly)^2) .* ts,rtol=2e-3)

#     # small heater case for one step
#     T0 = Tref + 0.0 # uniform initial temperature

#     temperature(integrator_plate) .= T0 # reset to uniform Tref
#     T_hist = []
#     Tgrid = temperature(integrator_plate);
#     Tgrid_ini = deepcopy(Tgrid)
#     tstep = 1e-2
#     ts = tstep:tstep:tstep

#     for i in ts
#         push!(T_hist,deepcopy(ADI_timemarching!(Tgrid,sys_plate_test,tstep)))
#     end

#     # test mean temperature (energy conservation)

#     T_hist_mean = mean(T_hist[1])
#     ΔT_mean = T_hist_mean - T0
#     ΔT_mean_analytical = power * tstep / (ρₛ * cₛ * 2Lx * 2Ly * plate_d)

#     @test isapprox(ΔT_mean,ΔT_mean_analytical,atol=1e-12,rtol=1e-4)

#     # test T rise/drop for one step with heater and condenser
#     T0 = Tref + 10.0 # uniform initial temperature

#     temperature(integrator_plate) .= T0 # reset to uniform Tref
#     T_hist = []
#     Tgrid = temperature(integrator_plate);
#     Tgrid_ini = deepcopy(Tgrid)
#     tstep = 1e-2
#     ts = tstep:tstep:tstep

#     for i in ts
#         push!(T_hist,deepcopy(ADI_timemarching!(Tgrid,sys_plate_test,tstep)))
#     end

#     maximum_T = maximum(T_hist[1])
#     minimum_T = minimum(T_hist[1])
#     ΔT_max = maximum_T - T0
#     ΔT_min = minimum_T - T0
#     total_heater_area = 0.5inches*0.5inches;
#     qe = power/total_heater_area
#     ΔT_max_analytical = qe * tstep / (ρₛ * cₛ * plate_d)


#     hc = sys_plate_test.qhdT[2].hc
#     qc = -hc*(T0-Tref)
#     ΔT_min_analytical = qc * tstep / (ρₛ * cₛ * plate_d)

#     @test isapprox(ΔT_max,ΔT_max_analytical,atol=1e-12,rtol=5e-3)
#     @test isapprox(ΔT_min,ΔT_min_analytical,atol=1e-12,rtol=5e-3)


#     # test energy rise/drop for one step with small heater and line source

#     sys_plate_test2 = ComputationalHeatTransfer.HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom)
#     # sys_tube = initialize_ohpsys(sys_plate_test2,p_fluid,power)
#     u_plate2 = newstate(sys_plate_test2) .+ Tref # initialize plate T field to uniform Tref
#     integrator_plate2 = init(u_plate2,tspan,sys_plate_test2) #

#     xylim = sys_plate_test2.grid.xlim
#     Lx = xylim[1][2] .+ 0.5Δx
#     Ly = xylim[2][2] .+ 0.5Δx

#     T0 = Tref + 0.0 # uniform initial temperature
#     q1D_value = -2.0 # W/m line heat flux
#     q1D = zeros(length(sys_plate_test2.qline[1].arccoord)) .+ q1D_value # W/m line heat flux
#     # q1D = deepcopy(zero(sys_plate_test2.qline[1].arccoord) .+ 1.0) # W/m line heat flux
#     L = sys_plate_test2.qline[1].arccoord[end]
#     set_linesource_strength!(sys_plate_test2,q1D)

#     temperature(integrator_plate2) .= T0 # reset to uniform Tref
#     T_hist = []
#     Tgrid = temperature(integrator_plate2);
#     Tgrid_ini = deepcopy(Tgrid)
#     tstep = 1e-2
#     ts = tstep:tstep:tstep

#     for i in ts
#         push!(T_hist,deepcopy(ADI_timemarching!(Tgrid,sys_plate_test2,tstep)))
#     end

#     T_hist_mean = mean(T_hist[1])
#     ΔT_mean = T_hist_mean - T0
#     ΔT_mean_analytical = (-L*q1D_value) * tstep / (ρₛ * cₛ * 2Lx * 2Ly * plate_d)

#     @test isapprox(ΔT_mean,ΔT_mean_analytical,atol=1e-12,rtol=1e-4)
                                                 
# end

sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

@testset "Film areas" begin

    d = sys_tube.tube.d
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    dXdt = sys_tube.liquid.dXdt
    Ac = sys_tube.tube.Ac

    δdep = 0.05*rand()*d

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2)
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2)
    δarea_dep = Ac .* (1 .- ((d .- 2*δdep) ./ d) .^ 2)
    @test all(δarea_start .== getδarea.(Ac,d,δstart))
    @test all(δarea_end .== getδarea.(Ac,d,δend))

    # No slugs are moving. Should only be equal to existing films
    Adep = getAdeposit(sys_tube,δdep)
    @test all(map((u,v) -> u[1]==v,Adep,δarea_end))
    @test all(map((u,v) -> u[2]==v,Adep,circshift(δarea_start,-1)))

    # Left slug interfaces are advancing. Should be equal to deposited film
    dXdt .= [(0.1,0.0) for i in eachindex(dXdt)]
    Adep = getAdeposit(sys_tube,δdep)
    @test all(map(u -> u[1]==δarea_dep,Adep))
    @test all(map((u,v) -> u[2]==v,Adep,circshift(δarea_start,-1)))

    # Right slug interfaces are advancing. Should be equal to deposited film
    dXdt .= [(0.0,-0.1) for i in eachindex(dXdt)]
    Adep = getAdeposit(sys_tube,δdep)
    @test all(map((u,v) -> u[1]==v,Adep,δarea_end))
    @test all(map(u -> u[2]==δarea_dep,Adep))

    vol = getVolumevapor(sys_tube)
    @test all(vol .> 0)

    ρv = sys_tube.propconvert.PtoD.(sys_tube.vapor.P)
    @test getMvapor(sys_tube) ≈ ρv.*vol

end


@testset "Mass functions" begin

    numofslugs = length(sys_tube.liquid.Xp)

    sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
    sys_tube.vapor.δend .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
    

    d = sys_tube.tube.d
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end
    Xp = sys_tube.liquid.Xp
    dXdt = sys_tube.liquid.dXdt
    Ac = sys_tube.tube.Ac
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρₗ
    closedornot = sys_tube.tube.closedornot

    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)
    Lliquidslug = XptoLliquidslug(Xp,L)
    Astart = getδarea.(Ac,d,δstart)
    Aend = getδarea.(Ac,d,δend)

    ρv = sys_tube.propconvert.PtoD.(sys_tube.vapor.P)

    # mass of vapor
    vol_vapor_analytical = Lvaporplug .* Ac .- Lfilm_start .* Astart .- Lfilm_end .* Aend
    M_vapor_analytical = ρv .* vol_vapor_analytical

    @test OscillatingHeatPipe.getMvapor(sys_tube) ≈ M_vapor_analytical

    # mass of liquid
    vol_liquid_analytical = Lliquidslug .* Ac
    M_liquid_analytical = ρₗ .* vol_liquid_analytical

    @test OscillatingHeatPipe.getMliquid(sys_tube) ≈ M_liquid_analytical

    # mass of films
    vol_film_start_analytical = Lfilm_start .* Astart
    M_film_start_analytical = ρₗ .* vol_film_start_analytical
    vol_film_end_analytical = Lfilm_end .* Aend
    M_film_end_analytical = ρₗ .* vol_film_end_analytical

    Mfilm_start,Mfilm_end = OscillatingHeatPipe.getMfilm(sys_tube)

    @test Mfilm_start ≈ M_film_start_analytical
    @test Mfilm_end ≈ M_film_end_analytical

    @test getMtotal(sys_tube) ≈ sum(M_vapor_analytical) + sum(M_liquid_analytical) + sum(M_film_start_analytical) + sum(M_film_end_analytical)

end

# some threshold values need to be changed for other applications
@testset "Hfilm" begin

    δmin = sys_tube.vapor.δmin;
    δthreshold = 5e-6
    δmax = 1e-4

    kₗ = sys_tube.vapor.k
    Hᵥ = sys_tube.vapor.Hᵥ
    δs = [1e-6;3e-6;1e-5;1.5e-4;3e-4]

    @test OscillatingHeatPipe.Hfilm.(δs,[sys_tube]) ≈ [0.0;(δs[2]-δmin)*(kₗ/δthreshold - Hᵥ)/(δthreshold-δmin);kₗ/δs[3];0.5*kₗ/δmax+1e-6;0.0]
end

