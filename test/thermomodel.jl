using Statistics
using Interpolations
using LinearAlgebra

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))


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
    
timestep_fixed(u,sys) = tstep

prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             # timestep_func=timestep_fourier
                                             timestep_func=timestep_fixed)

sys_plate = construct_system(prob)

numofslugs = 30
sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=rand(),g = 0 .* [-9.8,0.0])
sys_tube.liquid.dXdt = [zero.(X) .+ 1.0 for X in sys_tube.liquid.dXdt]# initial velocity of the slugs
sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 # initial velocity of the slugs

u_tube = newstate(sys_tube) # initialize OHP tube
# integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

@testset "dynamicsmodel test1" begin
    σ = sys_tube.liquid.σ
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρₗ
    μₗ = sys_tube.liquid.μₗ
    ad_fac = sys_tube.vapor.ad_fac
    d = sys_tube.tube.d
    Ac = sys_tube.tube.Ac
    peri = sys_tube.tube.peri
    Xp = sys_tube.liquid.Xp
    dXdt = sys_tube.liquid.dXdt
    δend = sys_tube.vapor.δend
    # characteristic bulk velocities for each liquid slug
    V = [mean(elem) for elem in sys_tube.liquid.dXdt]
    # get a characteristic Capilarry number based on the average velocities
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)

    δdep = Catoδ(d,Ca,adjust_factor=ad_fac)

    Lliquidslug = XptoLliquidslug(Xp,L)

    Adeposit = getAdeposit(sys_tube,δdep)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]
    
# The first test is to check the case where there is no heat transfer and an initial velocity
    uu_test1 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)

    @test all(uu_test1[1:2:2*numofslugs-1] .≈ Ac ./ (Ac .- Adeposit_left) .*  V) # velocity should remain the same
    @test all(uu_test1[2:2:2*numofslugs]   .≈ Ac ./ (Ac .- Adeposit_right) .*  V) # velocity should remain the same

    # get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug
    Re = ρₗ .* abs.(V) .* d ./ μₗ
    f_coefficient = f_churchill.(Re)
    dXdt_to_stress = -0.125 .* f_coefficient .* ρₗ .* V .* abs.(V)
    dVdt = dXdt_to_stress*peri .* Lliquidslug ./ lhs

    rhs_dLdt = V .* ((Ac ./ (Ac .- Adeposit_right) - Ac ./ (Ac .- Adeposit_left)) .*  V ) .* ρₗ.* Ac ./ lhs
    # @test uu_test1[2:2:2*numofslugs]   .== Ac ./ (Ac .- Adeposit_right) .*  V # velocity should remain the same
    @test all(uu_test1[2*numofslugs+1:2:4*numofslugs-1] .≈ dXdt_to_stress*peri .* Lliquidslug ./ lhs .- rhs_dLdt)

    @test all(uu_test1[2*numofslugs+2:2:4*numofslugs]   .== uu_test1[2*numofslugs+1:2:4*numofslugs-1])

    @test all(isapprox.(uu_test1[4*numofslugs+1:5*numofslugs],0.0,atol=1e-12)) #dMdt = 0

    @test all(isapprox.(uu_test1[5*numofslugs+1:6*numofslugs],0.0,atol=1e-12)) #dδstart/dt = 0

    F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
    F2_end = ρₗ .* Ac .* 4 .* δdep .* (d .- δdep) ./ (d^2)
    C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)
    Lfilm_end = sys_tube.vapor.Lfilm_end

    @test all(uu_test1[6*numofslugs+1:7*numofslugs] .≈ (F2_end .- F_end) .* Ac ./ (Ac .- Adeposit_left) .*  V ./ C_end ./ Lfilm_end)

    @test all(uu_test1[7*numofslugs+1:8*numofslugs] .≈ -uu_test1[2:2:2*numofslugs])

    @test all(uu_test1[8*numofslugs+1:9*numofslugs] .≈ uu_test1[1:2:2*numofslugs-1])
                                                 
end

@testset "dynamicsmodel test2" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue,g = 0 .* [-9.8,0.0])

    ΔT = rand()
    ΔT_array = [ΔT,-ΔT]
    sys_tube.vapor.P[1] = sys_tube.propconvert.TtoP(Tref-ΔT_array[1])
    sys_tube.vapor.P[2] = sys_tube.propconvert.TtoP(Tref-ΔT_array[2])

    P1 = sys_tube.vapor.P[1]
    P2 = sys_tube.vapor.P[2]
    L = sys_tube.tube.L
    Ac = sys_tube.tube.Ac
    d = sys_tube.tube.d
    ρₗ = sys_tube.liquid.ρₗ
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.propconvert.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    u_tube = newstate(sys_tube) # initialize OHP tube

    uu_test2 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)


    @test all(isapprox.(uu_test2[1:2:2*numofslugs-1],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test2[2:2:2*numofslugs],0.0,atol=1e-12)) # velocity should remain the same

    Lliquidslug = XptoLliquidslug(Xp,L)
    lhs = ρₗ*Ac .* Lliquidslug
    @test all(uu_test2[2*numofslugs+1:2:4*numofslugs-1] .≈ [1,-1] .* (P1-P2)*Ac ./ lhs) # velocity should remain the same
    @test all(uu_test2[2*numofslugs+2:2:4*numofslugs]  .== uu_test2[2*numofslugs+1:2:4*numofslugs-1])

    dMdt_latent_start = zeros(numofslugs)
    dMdt_latent_end = zeros(numofslugs)
    dMdt_latent_start_positive = zeros(numofslugs)
    dMdt_latent_end_positive = zeros(numofslugs)

    dMdt_latent_start .= Lfilm_start .* peri .* k ./ δstart ./ Hfg .* ΔT_array
    dMdt_latent_end .= Lfilm_end .* peri .* k ./ δend ./ Hfg .* ΔT_array
    dMdt_latent_start_positive = heaviside.(dMdt_latent_start) .* dMdt_latent_start
    dMdt_latent_end_positive = heaviside.(dMdt_latent_end) .* dMdt_latent_end

    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative = dMdt_latent_end .- dMdt_latent_end_positive

    dMdt_latent = dMdt_latent_start .+ dMdt_latent_end

    @test all(isapprox.(uu_test2[4*numofslugs+1:5*numofslugs],dMdt_latent,atol=1e-12)) #dMdt

    F_start = ρₗ .* Ac .* 4 .* δstart .* (d .- δstart) ./ (d^2)
    F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
    dLdt_start = -(ηplusvalue.*dMdt_latent_start_positive .+ ηminusvalue.*dMdt_latent_start_negative) ./ F_start
    dLdt_end = -(ηplusvalue.*dMdt_latent_end_positive .+ ηminusvalue.*dMdt_latent_end_negative) ./ F_end

    @test all(uu_test2[7*numofslugs+1:8*numofslugs] .≈ dLdt_start)

    @test all(uu_test2[8*numofslugs+1:9*numofslugs] .≈ dLdt_end)

    C_start = ρₗ .* Ac .* 4 .* (d .- 2δstart) ./ (d^2)
    C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)
    dδdt_start = -(dMdt_latent_start .+ dLdt_start.*F_start) ./ C_start ./ Lfilm_start
    dδdt_end = -(dMdt_latent_end .+ dLdt_end.*F_end) ./ C_end ./ Lfilm_end

    @test all(uu_test2[5*numofslugs+1:6*numofslugs] .≈ dδdt_start)

    @test all(uu_test2[6*numofslugs+1:7*numofslugs] .≈ dδdt_end)

end

@testset "gravitational force" begin

    angle = 0.5*pi
    gvec = sin(angle).* [-9.8,0.0]
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,g = gvec)

    Xwallarray = sys_tube.wall.Xarray

    x_inter = linear_interpolation(Xwallarray,ohp.x,extrapolation_bc = Line())
    y_inter = linear_interpolation(Xwallarray,ohp.y,extrapolation_bc = Line())

    L = sys_tube.tube.L
    Ac = sys_tube.tube.Ac
    d = sys_tube.tube.d
    ρₗ = sys_tube.liquid.ρₗ
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.propconvert.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    u_tube = newstate(sys_tube) # initialize OHP tube
    uu_test3 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)


    @test all(isapprox.(uu_test3[1:2:2*numofslugs-1],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test3[2:2:2*numofslugs],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test3[5*numofslugs+1:6*numofslugs],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test3[6*numofslugs+1:7*numofslugs],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test3[7*numofslugs+1:8*numofslugs],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test3[8*numofslugs+1:9*numofslugs],0.0,atol=1e-12)) # velocity should remain the same


    X1 = [elem[1] for elem in Xp]
    X2 = [elem[2] for elem in Xp]
    hdiff = x_inter.(X1) .- x_inter.(X2)

    Lliquidslug = XptoLliquidslug(Xp,L)
    lhs = ρₗ*Ac .* Lliquidslug

    @test all(uu_test3[2*numofslugs+1:2:4*numofslugs-1] .≈ ρₗ .* Ac .* norm(gvec) .* hdiff ./ lhs) # velocity should remain the same
    @test all(uu_test3[2*numofslugs+2:2:4*numofslugs]  .== uu_test3[2*numofslugs+1:2:4*numofslugs-1])
end
