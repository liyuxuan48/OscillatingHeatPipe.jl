using Statistics
using JLD2

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))


using OscillatingHeatPipe # our main package

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

    dXp = mod(Xp_new[2][1] - Xp_new[1][end],L) - 0.8*OscillatingHeatPipe.DEFAULT_LIQUID_MERGE_FRAC*L_newbubble
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
    @test integrator_plate.t == tstep

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
    u_tube_0 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_0 = init(u_tube_0,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test vaporMergingCondition(integrator_tube_0.u,integrator_tube_0.t,integrator_tube_0) == false

    L_newbubble= sys_tube.wall.L_newbubble
    L = sys_tube.tube.L

    Xp_old = sys_tube.liquid.Xp
    Xp_new = deepcopy(Xp_old)

    dXp = mod(Xp_new[2][2] - Xp_new[2][1],L) - 1.1*OscillatingHeatPipe.DEFAULT_VAPOR_MERGE_FRAC*L_newbubble
    Xp_new[2] =  (Xp_new[2][1], Xp_new[2][2] - dXp) # make two liquid slugs not close enough to merge

    sys_tube.liquid.Xp = deepcopy(Xp_new)

    u_tube_1 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_1 = init(u_tube_1,tspan,deepcopy(sys_tube)); # construct integrator_tube

    @test vaporMergingCondition(integrator_tube_1.u,integrator_tube_1.t,integrator_tube_1) == false

    dXp = mod(Xp_new[2][2] - Xp_new[2][1],L) - 0.9*OscillatingHeatPipe.DEFAULT_VAPOR_MERGE_FRAC*L_newbubble
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

    u_plate = init_sol(sys_plate)# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    @test all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-12))

    timemarching!(integrator_tube,integrator_plate,tstep)

    println(integrator_plate.t)

    # because they exchange data after time marching, the wall temperature should no longer be Tref
    @test !all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-12))
    @test integrator_plate.t == tstep

    # record system state at t_n
    integrator_tube_n = deepcopy(integrator_tube)
    save("Coupling_int_plate.jld2","integrator_plate",integrator_plate)
    integrator_plate_n = load("Coupling_int_plate.jld2")["integrator_plate"]
    t_n = deepcopy(integrator_tube.t)

    # change plate temperature and create a new integrator_plate
    integrator_plate_n.u .+= 20.0 # increase plate temperature by 20K
    timemarching!(integrator_tube,integrator_plate,tstep)
    timemarching!(integrator_tube_n,integrator_plate_n,tstep)

    # because they exchange data after time marching, this change of plate 
    # temperature should not affect the wall temperature of the tube
    @test !all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-12))
    @test all(isapprox.(integrator_tube_n.u,integrator_tube.u, rtol=1e-12))

    store!(SimuResult,integrator_tube,integrator_plate)

    #make sure after time marching, the system state is up-to-date with the latest u
    p_old = deepcopy(integrator_tube.p)
    OscillatingHeatPipe.exchangepinfo!(integrator_tube,integrator_plate)

    L = p_old.tube.L
    randX = rand()*L
    @test p_old.mapping.θ_interp_walltoliquid(randX) == integrator_tube.p.mapping.θ_interp_walltoliquid(randX)
    @test all(p_old.vapor.P .== integrator_tube.p.vapor.P)
    @test all(p_old.liquid.θarrays .== integrator_tube.p.liquid.θarrays)
    @test p_old.wall.θarray == integrator_tube.p.wall.θarray
    @test SimuResult.tube_hist_θwall[end] == integrator_tube.p.wall.θarray
    @test SimuResult.tube_hist_u[end] == integrator_tube.u
    @test SimuResult.tube_hist_t[end] == integrator_tube.t
end

@testset "Coupling Scheme Identification" begin
    # 1. Setup a fresh system for identification
    # We reuse the variables defined in the previous scope if available, 
    # or you can re-initialize here to ensure total independence.

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-3     # actrual time marching step

    u_plate = init_sol(sys_plate)# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    tspan = (0.0, 1.0); # start time and end time

    # 2. Advance a few steps to reach a non-zero state
    t_step_id = 1e-3
    for _ in 1:3
        timemarching!(integrator_tube, integrator_plate, t_step_id)
    end

    # 3. Record baseline state at t_n
    save("Coupling_int_plate.jld2","integrator_plate",integrator_plate)
    integrator_tube_n = deepcopy(integrator_tube)
    # integrator_plate_n = deepcopy(integrator_plate)
    t_n = deepcopy(integrator_tube.t)

    # Get reference pressure at t_{n+1} without any perturbation
    timemarching!(integrator_tube, integrator_plate, t_step_id)
    p_reference = deepcopy(getcurrentsys_nowall!(integrator_tube.u, integrator_tube.p).vapor.P)

    # 4. Reset to t_n for the perturbation test
    integrator_tube = deepcopy(integrator_tube_n)
    integrator_plate = load("Coupling_int_plate.jld2")["integrator_plate"]

    # 5. Apply thermal perturbation (+50K) to the plate
    # This happens immediately before the marching call
    integrator_plate.u .+= 50.0 
    # 6. Run the marching step with the perturbation
    timemarching!(integrator_tube, integrator_plate, t_step_id)
    p_perturbed = getcurrentsys_nowall!(integrator_tube.u, integrator_tube.p).vapor.P

    # 7. Identify coupling logic
    # Jacobian (Synchronous/Parallel): Fluid solver uses values from t_n, ignoring the +50K jump.
    # Gauss-Seidel (Sequential): Fluid solver sees the +50K jump and adjusts P immediately.
    is_jacobian = all(isapprox.(p_reference, p_perturbed, atol=1e-12))

    if is_jacobian
        println("\n[Coupling Test]: Jacobian Scheme Identified")
        println(" - Logic: Synchronous/Parallel execution.")
        println(" - Observation: Fluid solver is independent of plate changes within the same step.")
    else
        println("\n[Coupling Test]: Gauss-Seidel Scheme Identified")
        println(" - Logic: Sequential/Successive execution.")
        println(" - Observation: Fluid solver responded immediately to the plate perturbation.")
    end

    #  Explicitly test based on your current goal (Optional)
    @test is_jacobian 
end


@testset "Multithreading" begin

    @test Threads.nthreads() > 1
    println("This test requires multiple threads (>1) to run.")
    println("num of threads= ",Threads.nthreads())

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    tspan = (0.0, 2.0); # start time and end time
    dt_record = 2.0   # saving time interval

    tstep = 4e-4     # actrual time marching step

    u_plate = init_sol(sys_plate) .+10.0# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    integrator_tube_parallel = deepcopy(integrator_tube)
    save("Coupling_int_plate.jld2","integrator_plate",integrator_plate)

    for t in tspan[1]:tstep:tspan[2]
    
        timemarching!(integrator_tube,integrator_plate,tstep,force_sequential=true)

    end

    integrator_plate_parallel = load("Coupling_int_plate.jld2")["integrator_plate"]
    

    for t in tspan[1]:tstep:tspan[2]
    
        timemarching!(integrator_tube_parallel,integrator_plate_parallel,tstep)

    end


    # SimuResult = SimulationResult(integrator_tube,integrator_plate);

    @test integrator_tube_parallel.u == integrator_tube.u
    @test integrator_plate_parallel.u == integrator_plate.u
    @test integrator_tube.p.cache.boil_hist == integrator_tube_parallel.p.cache.boil_hist
end