#= 
simulation.jl
    Create a simulation
    Run using
        julia ... -- simulation.jl [mode folder] [paramfile] [init time] [fixed evolution time] [runtime] [output folder]
=#
using Oceananigans
include("InputHelper.jl")
using .InputHelper

(mode_folder, param_file, init_time, fix_evo_time, run_time, output_folder) = ARGS
const (simulation_parameters, mp) = readparams(param_file)

# Include the components
include("parameters.jl")
include("filament_state.jl")
include("sponge_layer.jl")
include("boundary_conditions.jl")

include("$mode_folder/outputs.jl")
include("$mode_folder/tracers.jl")
# Build parameters
sp = (; simulation_parameters, create_simulation_parameters(simulation_parameters))
(filament_b, filament_v) = get_filament_state(sp)

# To ensure that the grids are isotropic in horizontal,
horizontal_aspect_ratio = sp.Ny / sp.Nx

# Calculate the z spacing
K = 1.05
z_faces(k) = -sp.Lz * (K - K^k) / (K - K^sp.Nz);

grid = RectilinearGrid(GPU(),
        size=(sp.Nx, sp.Ny, sp.Nz),
        x=(-2sp.L, 2sp.L),
        y=(-2horizontal_aspect_ratio * sp.L, 2horizontal_aspect_ratio * sp.L),
        z=z_faces,
        topology=(Periodic, Periodic, Bounded));

@info grid

@info "Building velocity and tracer forcing functions"
tracers = (:b, additional_tracers(sp, mp)...)
boundary_conditions = (get_boundary_conditions(sp)..., additional_tracer_bcs(sp, mp)...)
forcing = (get_sponge_layer_forcing(sp; σ=sp.σ).., additional_tracer_forcing(sp, mp)...)
@info tracers

# Turbulence closure
closure = if sp.ν == nothing
    Oceananigans.TurbulenceClosures.SmagorinskyLilly(Pr=sp.Pr)
else
    Oceananigans.TurbulenceClosures.ScalarDiffusivity(sp.ν; b=sp.ν / sp.Pr)
end

@info "Creating model"
model = NonhydrostaticModel(;
        grid,
        coriolis = FPlane(sp.f),
        closure,
        buoyancy = BuoyancyTracer(),
        tracers,
        boundary_conditions,
        forcing,
    )

@info model

# set initial conditions as the conditions far from the filament
# and a random secondary circulation
u₀(x, y, z) = 1e-5*rand()
v₀(x, y, z) = v_filament(100sp.L, z)
w₀(x, y, z) = 1e-5*rand()
b₀(x, y, z) = b_filament(100sp.L, z)

set!(model; u=u₀, v=v₀, w=w₀, b=b₀, additional_tracer_initial_conditions(sp, mp)...)

# How to initialise