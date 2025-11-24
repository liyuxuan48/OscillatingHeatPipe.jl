using Statistics

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))


using OscillatingHeatPipe # our main package

ρₛ = 2730 # material density [kg/m^3]
cₛ  = 8.93e02 # material specific heat [J/kg K]
kₛ  = 1.93e02 # material heat conductivity
αₛ = kₛ/ρₛ/cₛ
    
dₛ = 1.5e-3

Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref) # This function relies on CoolProp.jl package

power = 30 # [W], total power
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


    tspan = (0.0, 5.0); # start time and end time
    dt_record = 0.01   # saving time interval

    tstep = 1e-3     # actrual time marching step
    
timestep_fixed(u,sys) = tstep

prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             # timestep_func=timestep_fourier
                                             timestep_func=timestep_fixed)

sys_plate = construct_system(prob)


# test 1, boiling callbacks
@testset "merging callbacks" begin

    A_SMALL_FRAC = 1e-2

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)
    u_tube_1 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_1 = init(u_tube_1,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test merging_condition(integrator_tube_1.u,integrator_tube_1.t,integrator_tube_1) == false

    L_newbubble= sys_tube.wall.L_newbubble
    L = sys_tube.tube.L

    Xp_old = sys_tube.liquid.Xp
    Xp_new = deepcopy(Xp_old)

    dXp = mod(Xp_new[2][1] - Xp_new[1][end],L) - 0.4*L_newbubble
    @test  dXp > 0.0

    Xp_new[2] =  (Xp_new[2][1] - dXp,Xp_new[2][2] - dXp) # make two vapor plugs close enough to merge

    Lfilm_start= sys_tube.vapor.Lfilm_start
    Lfilm_end= sys_tube.vapor.Lfilm_end

    Lfilm_start[2] = A_SMALL_FRAC*L_newbubble
    Lfilm_end[2] = A_SMALL_FRAC*L_newbubble

    sys_tube.liquid.Xp = deepcopy(Xp_new)
    sys_tube.vapor.Lfilm_start = deepcopy(Lfilm_start)
    sys_tube.vapor.Lfilm_end = deepcopy(Lfilm_end)

    u_tube_2 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_2 = init(u_tube_2,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test merging_condition(integrator_tube_2.u,integrator_tube_2.t,integrator_tube_2) == true

    # println("Before modifying Xp: ", sys_tube.liquid.Xp)
    # println("After modifying Xp: ", Xp_new)

    #   ### combine inner tube and plate together

    u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

    SimuResult = SimulationResult(integrator_tube,integrator_plate);


    # println(integrator_tube.p.liquid.Xp)


    p_old = deepcopy(integrator_tube.p)

    Mvapor_old = sum(getMvapor(p_old))
    Mfilm_old = sum(sum.(getMfilm(p_old)))
    Mliquid_old = sum(getMliquid(p_old))

    merging_affect!(integrator_tube)
    getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

    p_new = deepcopy(integrator_tube.p)

    @test length(p_new.liquid.Xp) == length(p_old.liquid.Xp) - 1
    @test isapprox(sum(XptoLliquidslug(p_new.liquid.Xp,p_new.tube.L)), sum(XptoLliquidslug(p_old.liquid.Xp,p_old.tube.L)), rtol=2e-3)
    @test p_new.liquid.Xp[2:end] == p_old.liquid.Xp[3:end]
    @test p_new.liquid.dXdt[2:end] == p_old.liquid.dXdt[3:end]


    L_liquidslug_old = XptoLliquidslug(p_old.liquid.Xp,p_old.tube.L)
    L_liquidslug_new = XptoLliquidslug(p_new.liquid.Xp,p_new.tube.L)
    @test p_new.liquid.dXdt[1][1] ≈ (p_old.liquid.dXdt[1][1] * L_liquidslug_old[1] + p_old.liquid.dXdt[2][1] * L_liquidslug_old[2])/(L_liquidslug_old[1]+L_liquidslug_old[2]) 

    Mvapor_new = sum(getMvapor(p_new))
    Mfilm_new = sum(sum.(getMfilm(p_new)))
    Mliquid_new = sum(getMliquid(p_new))

    @test isapprox(Mvapor_new+Mfilm_new+Mliquid_new, Mvapor_old+Mfilm_old+Mliquid_old, rtol=1e-10)



    # @showprogress for t in tspan[1]:tstep:tspan[1]

    #     timemarching!(integrator_tube,integrator_plate,tstep)

    #     if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
    #         store!(SimuResult,integrator_tube,integrator_plate)
    #     end

    # end
end

@testset "boiling callbacks" begin

    numofslugs = 2

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,nucleatenum=2,
        boil_waiting_time=1e-2)

    Xp = sys_tube.liquid.Xp
    L = sys_tube.tube.L

    Xp1_mid = mod(Xp[1][1] + mod((Xp[1][2]- Xp[1][1]),L)/2,L)
    Xp2_vapor_mid = mod(Xp[1][2] + mod((Xp[2][1]- Xp[1][2]),L)/2,L)

    sys_tube.wall.Xstations=[Xp1_mid,Xp2_vapor_mid]


    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = OscillatingHeatPipe.DEFAULT_BOIL_SCAN_INTERVAL
    
    Rn = sys_tube.wall.Rn
    d = sys_tube.tube.d
    TtoP = sys_tube.propconvert.TtoP
    PtoT = sys_tube.propconvert.PtoT


    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    Tavg = mean(PtoT(integrator_tube.p.vapor.P)) # initial vapor average temperature equals Tref
    ΔT = RntoΔT(Rn,Tavg,fluid_type,d,TtoP)
    

    u_plate = init_sol(sys_plate) .+ 10.0ΔT # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate


    SimuResult = SimulationResult(integrator_tube,integrator_plate);

 
    timemarching!(integrator_tube,integrator_plate,tstep)


    @test integrator_tube.p.wall.boiltime_stations == [0.0,tstep]

    println

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test integrator_tube.p.wall.boiltime_stations == [2*tstep,2*tstep]


    u_plate = init_sol(sys_plate) .+ 1.1ΔT# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate


    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube
    integrator_tube.p.wall.θarray = OscillatingHeatPipe.temperature_linesource(integrator_plate)
    

    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    boiling_affect!(integrator_tube)
    getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

    @test integrator_tube.p.liquid.Xp == sys_tube.liquid.Xp

    integrator_tube.t = integrator_tube.p.wall.boil_interval + tstep

    p_old = deepcopy(integrator_tube.p)
    Mvapor_old = sum(getMvapor(p_old))
    Mfilm_old = sum(sum.(getMfilm(p_old)))
    Mliquid_old = sum(getMliquid(p_old))

    boiling_affect!(integrator_tube)
    getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

    @test length(integrator_tube.p.liquid.Xp) == length(sys_tube.liquid.Xp) + 1


    p_new = deepcopy(integrator_tube.p)
    Mvapor_new = sum(getMvapor(p_new))
    Mfilm_new = sum(sum.(getMfilm(p_new)))
    Mliquid_new = sum(getMliquid(p_new))
    @test isapprox(Mvapor_new+Mfilm_new+Mliquid_new, Mvapor_old+Mfilm_old+Mliquid_old, rtol=1e-4)

    u_plate = init_sol(sys_plate) .+ 0.9ΔT# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate
    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube
    integrator_tube.p.wall.θarray = OscillatingHeatPipe.temperature_linesource(integrator_plate)
    integrator_tube.t = integrator_tube.p.wall.boil_interval + tstep

    boiling_affect!(integrator_tube)
    getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

    @test integrator_tube.p.liquid.Xp == sys_tube.liquid.Xp

end



@testset "fixdx callbacks" begin



    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-2     # actrual time marching step


 sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

@test fixdx_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == false

sys_tube.liquid.Xarrays[1] = OscillatingHeatPipe.constructoneXarray(integrator_tube.p.liquid.Xp[1],3*length(integrator_tube.p.liquid.Xarrays[1]),integrator_tube.p.tube.L)
sys_tube.liquid.θarrays[1] = sys_tube.liquid.Xarrays[1] .* 0 .+ Tref
u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube


@test fixdx_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == true

fixdx_affect!(integrator_tube)
getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

@test length(integrator_tube.p.liquid.Xarrays[1]) < length(sys_tube.liquid.Xarrays[1])
@test fixdx_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == false

end

@testset "slugbc callbacks" begin



    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-2     # actrual time marching step


 sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)
 
u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube
# for i in eachindex(integrator_tube.p.vapor.P)
#     integrator_tube.p.liquid.θarrays[i] = integrator_tube.p.liquid.θarrays[i] .+ rand() .* 1e3
# end

integrator_tube.p.vapor.P = integrator_tube.p.vapor.P + rand(length(integrator_tube.p.vapor.P)) .* 1e3


@test slugbc_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == true # always true


T_first_pt = [x[1] for x in integrator_tube.p.liquid.θarrays]
T_last_pt = [x[end] for x in integrator_tube.p.liquid.θarrays]

PtoT = integrator_tube.p.propconvert.PtoT

@test !isapprox(PtoT.(integrator_tube.p.vapor.P),T_first_pt, rtol=1e-10)
@test !isapprox(PtoT.(integrator_tube.p.vapor.P),circshift(T_last_pt,1), rtol=1e-10)

slugbc_affect!(integrator_tube)
getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

@test isapprox(PtoT.(integrator_tube.p.vapor.P),T_first_pt, rtol=1e-10)
@test isapprox(PtoT.(integrator_tube.p.vapor.P),circshift(T_last_pt,1), rtol=1e-10)

end

# test 1, boiling callbacks
@testset "vapormerging callbacks" begin

    tspan = (0.0, 5.0); # start time and end time
    dt_record = 0.01   # saving time interval

    tstep = 1e-3     # actrual time marching step


    A_SMALL_FRAC = 1e-2

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)
    u_tube_1 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_1 = init(u_tube_1,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test vaporMergingCondition(integrator_tube_1.u,integrator_tube_1.t,integrator_tube_1) == false

    L_newbubble= sys_tube.wall.L_newbubble
    L = sys_tube.tube.L

    Xp_old = sys_tube.liquid.Xp
    Xp_new = deepcopy(Xp_old)

    dXp = mod(Xp_new[2][2] - Xp_new[2][1],L) - 0.4*L_newbubble
    Xp_new[2] =  (Xp_new[2][1], Xp_new[2][2] - dXp) # make two liquid slugs close enough to merge

    sys_tube.liquid.Xp = deepcopy(Xp_new)

    u_tube_2 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_2 = init(u_tube_2,tspan,deepcopy(sys_tube)); # construct integrator_tube

    @test vaporMergingCondition(integrator_tube_2.u,integrator_tube_2.t,integrator_tube_2) == true

#     # println("Before modifying Xp: ", sys_tube.liquid.Xp)
#     # println("After modifying Xp: ", Xp_new)

#     #   ### combine inner tube and plate together

#     u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
#     integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

#     SimuResult = SimulationResult(integrator_tube,integrator_plate);


#     # println(integrator_tube.p.liquid.Xp)


    p_old = deepcopy(integrator_tube.p)

    Mvapor_old = sum(getMvapor(p_old))
    Mfilm_old = sum(sum.(getMfilm(p_old)))
    Mliquid_old = sum(getMliquid(p_old))

    vaporMergingAffect!(integrator_tube)
    getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)

    p_new = deepcopy(integrator_tube.p)

    @test length(p_new.liquid.Xp) == length(p_old.liquid.Xp) - 1
    @test isapprox(sum(XptoLliquidslug(p_new.liquid.Xp,p_new.tube.L)), sum(XptoLliquidslug(p_old.liquid.Xp,p_old.tube.L)), rtol=2e-3)

    Mvapor_new = sum(getMvapor(p_new))
    Mfilm_new = sum(sum.(getMfilm(p_new)))
    Mliquid_new = sum(getMliquid(p_new))

    @test isapprox(Mvapor_new+Mfilm_new+Mliquid_new, Mvapor_old+Mfilm_old+Mliquid_old, rtol=1e-10)

end

@testset "weakly coupled time marching, Jacobian style" begin

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-3     # actrual time marching step

    u_plate = init_sol(sys_plate) .+ 10.0# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    @test all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-10))

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test !all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-10))

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test !all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-10))

    store!(SimuResult,integrator_tube,integrator_plate)

    #make sure after time marching, the system state is up-to-date with the latest u
    p_old = deepcopy(integrator_tube.p)
    OscillatingHeatPipe.exchangepinfo!(integrator_tube,integrator_plate)

    L = p_old.tube.L
    randX = rand()*L
    @test p_old.mapping.θ_interp_walltoliquid(randX) == integrator_tube.p.mapping.θ_interp_walltoliquid(randX)
    @test all(p_old.vapor.P .== integrator_tube.p.vapor.P)
    @test SimuResult.tube_hist_θwall[end] == integrator_tube.p.wall.θarray
    @test SimuResult.tube_hist_u[end] == integrator_tube.u
    @test SimuResult.tube_hist_t[end] == integrator_tube.t
end


@testset "Multithreading" begin

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-3     # actrual time marching step

    u_plate = init_sol(sys_plate) .+ 10.0# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube


    integrator_plate_parallel = deepcopy(integrator_plate)
    integrator_tube_parallel = deepcopy(integrator_tube)

    @test integrator_tube_parallel.u == integrator_tube.u
    @test integrator_plate_parallel.u == integrator_plate.u

    # SimuResult = SimulationResult(integrator_tube,integrator_plate);

    for t in tspan[1]:tstep:tspan[1]
    
        timemarching!(integrator_tube,integrator_plate,tstep,force_sequential=true)
        timemarching!(integrator_tube_parallel,integrator_plate_parallel,tstep)

    end

    @test integrator_tube_parallel.u == integrator_tube.u
    @test integrator_plate_parallel.u == integrator_plate.u
end