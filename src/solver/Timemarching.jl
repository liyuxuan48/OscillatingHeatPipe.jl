import ConstrainedSystems: init, solve

export ODE_innertube,ODE_steadyfilm,timemarching!

"""
    ODE_innertube(u,p,t)
The ODE function for the tube system, to be used in `ODEProblem`.
It takes the current state vector `u`, the parameter struct `p` (which is a `PHPSystem`), and the current time `t`.
It returns the time derivative of the state vector `du`.
`p` is used as the input for the dynamicsmodel and liquidmodel to get the current system state. Needs to be updated
    everytime this function is called. When called by the ODE solver, `p` is not updated automatically when `u` is 
    changed, this can be because of callbacks or substeps of the ODE solver.
getcurrentsys_nowall! is used to update the system state from `u` and `p`, without updating wall temperature.
"""

function ODE_innertube(u,p,t)

    sys_init = p

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    # These updated don't have wall temperature updated yet so getcurrentsys_nowall! can be used.
    newsys = getcurrentsys_nowall!(u,sys_init)

    dynamicsdu = dynamicsmodel(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]
    
    return(du)

end

# function ODE_steadyfilm(u,p,t)

#     index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

#     newsys = getcurrentsys_nowall!(u,p)

#     dynamicsdu = dynamicsmodel_steadyfilm(u[1:index_dynamics_end-1],newsys)

#     liquiddu = duliquidθtovec(liquidmodel(newsys))

#     du = [dynamicsdu;liquiddu]

#     return(du)

# end

"""
    temperature_linesource(integrator_plate::ODEIntegrator)

Get the temperature distribution on the plate from the plate integrator, to be used as wall temperature boundary condition for the tube system.
"""

function temperature_linesource(integrator_plate::ODEIntegrator)

    T = ComputationalHeatTransfer.temperature(integrator_plate)
    ohp_E = integrator_plate.p.extra_cache.fcache[end].region_cache.cache.E;
    ohp_E*T
end

"""
    timemarching!(integrator_tubez::ODEIntegrator,integrator_plate::ODEIntegrator,tstep::Float64)

Perform one time step of size `tstep` for both the tube and plate integrators, exchanging information between them.
This is a Jacobi style weakly coupled alternate time marching (exchange info at the same time after each step)
"""

function timemarching!(integrator_tube::ODEIntegrator,integrator_plate::ODEIntegrator,tstep::Float64)

    # initial exchange info at t=0
    if isapprox(integrator_tube.t,0.0)
        exchangepinfo!(integrator_tube,integrator_plate) # initial exchange info
        println("Initial exchange info done.")
    end

    # use two threads to step the two integrators together
    @sync begin
         Threads.@spawn step!(integrator_tube,tstep,true);
         Threads.@spawn step!(integrator_plate,tstep,true);
    end

    # exchange info after stepping at each time step
    exchangepinfo!(integrator_tube,integrator_plate)

    integrator_tube,integrator_plate
end

function init(u_tube::Vector{Float64},tspan::Tuple{Any, Any},sys_tube::PHPSystem,kwargs...)
    
    cb_boiling =  DiscreteCallback(boiling_condition,boiling_affect!)
    cb_liquidmerging =  DiscreteCallback(merging_condition,merging_affect!)
    cb_vapormerging = DiscreteCallback(vaporMergingCondition,vaporMergingAffect!)
    cb_fixdx =  DiscreteCallback(fixdx_condition,fixdx_affect!)
    cb_slugbc = DiscreteCallback(slugbc_condition,slugbc_affect!)
    cbst = CallbackSet(cb_fixdx,cb_boiling,cb_vapormerging,cb_liquidmerging,cb_slugbc);
    
    # cb_liquidmerging =  DiscreteCallback(merging_condition,merging_affect!)
    # cbst = CallbackSet(cb_liquidmerging);

    prob = ODEProblem(ODE_innertube, u_tube, tspan, sys_tube,kwargs...) # construct integrator_tube problem
    return init(prob, alg=RK4(),dt=1e-3,save_on=false, callback=cbst,maxiters=1e10,kwargs...)
    # return init(prob, alg=RK4(),dt=1e-3,save_on=false, maxiters=1e10,kwargs...)
end

function exchangepinfo!(integrator_tube::ODEIntegrator,integrator_plate::ODEIntegrator)

    # update interpolated wall temperature to tube system
    # currentsys.wall.θarray = temperature_linesource(integrator_plate)
    integrator_tube.p.wall.θarray = temperature_linesource(integrator_plate)
    # update the system of integrator_tube.p from the latest integrator_tube.u, without updating wall temperature
    currentsys = getcurrentsys_nowall!(integrator_tube.u,integrator_tube.p)


    # update heat flux to plate system
    qtmp = sys_to_heatflux(currentsys)
    integrator_plate.p.phys_params["ohp_flux"] = -qtmp # q' source per length
    
    integrator_tube,integrator_plate
end