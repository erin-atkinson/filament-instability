#= 
simulation.jl
    Create a simulation
    Run using
        julia ... -- simulation.jl [mode folder] [paramfile] [init time] [runtime] [writefreq] [output folder]
=#

ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0

using Oceananigans
using JLD2
using Oceananigans.Operators
using SpecialFunctions: erf
using Statistics
using CUDA

include("InputHelper.jl")
using .InputHelper

(mode_folder, param_file, init_time, run_time, write_freq, output_folder) = ARGS
init_time = Meta.parse(init_time)
run_time = Meta.parse(run_time)
write_freq = Meta.parse(write_freq)

const (simulation_parameters, mp) = readparams(param_file)

# Include the components
include("parameters.jl")
include("filament_state.jl")
include("sponge_layer.jl")
include("boundary_conditions.jl")
include("z_faces.jl")

include("$mode_folder/outputs.jl")
include("$mode_folder/tracers.jl")
# Build parameters
sp = (; simulation_parameters..., create_simulation_parameters(simulation_parameters)...)
@info sp
(b_filament, v_filament) = get_filament_state(sp)

# To ensure that the grids are isotropic in horizontal,
horizontal_aspect_ratio = sp.Ny / sp.Nx

grid = RectilinearGrid(GPU(),
        size=(sp.Nx, sp.Ny, sp.Nz),
        x=(-5sp.L, 5sp.L),
        y=(-5horizontal_aspect_ratio * sp.L, 5horizontal_aspect_ratio * sp.L),
        z=get_z_faces(sp),
        topology=(Periodic, Periodic, Bounded));

@info grid

@info "Building velocity and tracer forcing functions"
tracers = (:b, additional_tracers(sp, mp)...)
boundary_conditions = (; get_boundary_conditions(sp, init_time)..., additional_tracer_bcs(sp, mp, init_time)...)
forcing = (; get_sponge_layer_forcing(sp; σ=sp.σ, c=sp.c)..., additional_tracer_forcing(sp, mp, init_time)...)
@info tracers

# Turbulence closure
closure = if sp.ν == nothing
    Oceananigans.TurbulenceClosures.SmagorinskyLilly(Pr=sp.Pr)
else
    # This actually doesn't really work due to the changing grid size, will probably replace
    Oceananigans.TurbulenceClosures.ScalarDiffusivity(ν=sp.ν, κ=sp.ν / sp.Pr)
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

@info "Setting model state"
# set initial conditions as the conditions far from the filament
# and a random vertical velocity at the surface
u₀(x, y, z) = 0 #1e-2*randn() * (tanh((z + sp.H) / (sp.λ * sp.H)) + 1)
v₀(x, y, z) = v_filament(100sp.L, z) #+ 1e-2*randn() * (tanh((z + sp.H) / (sp.λ * sp.H)) + 1)
w₀(x, y, z) = 1e-4*randn() * (tanh(z / (sp.λ * sp.H)) + 1)
b₀(x, y, z) = b_filament(100sp.L, z)

set!(model; u=u₀, v=v₀, w=w₀, b=b₀, additional_tracer_initial_conditions(sp, mp)...)

# Create a default output that saves the average state of the simulation u, v, w, b fields
@info "Creating simulation"
simulation = Simulation(model, Δt=1/sp.Np, stop_time=init_time)

# Progress and timestepper wizards
progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", time(sim), ", step: ", sim.Δt)
wizard = TimeStepWizard(cfl=0.5, diffusive_cfl=0.5, min_Δt = 0.1/sp.Np, max_Δt = 10/sp.Np)
simulation.callbacks[:progress] = Callback(progress, TimeInterval(10/(sp.f*write_freq)))
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

(u, v, w) = model.velocities
b = model.tracers.b

# Create output folder
output_folder = if isdir("$output_folder")
    i = [1]
    while isdir("$(output_folder)_$(i[1])")
        i[1] = i[1] + 1
    end
    "$(output_folder)_$(i[1])"
else
    output_folder
end

mkdir(output_folder)
@info "Creating output writers at $output_folder"

# Mean fields
u_dfm = Field(Average(u; dims=2))
v_dfm = Field(Average(v; dims=2))
w_dfm = Field(Average(w; dims=2))
b_dfm = Field(Average(b; dims=2))

simulation.output_writers[:mean_state] = JLD2OutputWriter(
    model,
    (; u_dfm, v_dfm, w_dfm, b_dfm),
    filename = "$output_folder/down_front_mean.jld2",
    schedule = TimeInterval(1/(sp.f*write_freq)),
    overwrite_existing = true
)

# Save parameters to a file
jldopen("$output_folder/parameters.jld2", "a") do file
    file["parameters/simulation"] = sp
    file["parameters/init_time"] = init_time
    file["parameters/run_time"] = run_time
end

# Insert additional output code
@additional_outputs! simulation
@info simulation
# run the simulation for initialisation phase
@info "Initialising boundary layer turbulence"
run!(simulation)
if run_time == 0
    return
end

# set the filament
@info "Setting filament state"
# Need to edit GPU contents
CUDA.@allowscalar let xsᶜᶠᶜ = Array(xnodes(grid, Center(), Face(), Center())),
    zsᶜᶠᶜ = Array(znodes(grid, Center(), Face(), Center())),
    xsᶜᶜᶜ = Array(xnodes(grid, Center(), Center(), Center())),
    zsᶜᶜᶜ = Array(znodes(grid, Center(), Center(), Center()))
    
    v_mean_new = [v_filament(x, z) for x in xsᶜᶠᶜ, y in [1], z in zsᶜᶠᶜ]
    b_mean_new = [b_filament(x, z) for x in xsᶜᶜᶜ, y in [1], z in zsᶜᶜᶜ]
    
    # Get the current state
    u_new = Array(interior(model.velocities.u))
    v_new = Array(interior(model.velocities.v))
    w_new = Array(interior(model.velocities.w))
    b_new = Array(interior(model.tracers.b))
    
    # Remove any mean flow
    u_new .-= mean(u_new; dims=2)
    v_new .-= mean(v_new; dims=2)
    w_new .-= mean(w_new; dims=2)
    b_new .-= mean(b_new; dims=2)
    
    # Add a filament on
    v_new .+= v_mean_new
    b_new .+= b_mean_new
    
    set!(model; u=u_new, v=v_new, w=w_new, b=b_new)
    # Should also remove tendencies?
    #=
    model.timestepper.Gⁿ.u .= 0
    model.timestepper.Gⁿ.v .= 0
    model.timestepper.Gⁿ.w .= 0
    model.timestepper.Gⁿ.b .= 0
    
    model.timestepper.G⁻.u .= 0
    model.timestepper.G⁻.v .= 0
    model.timestepper.G⁻.w .= 0
    model.timestepper.G⁻.b .= 0
    =#
    # DEBUG
    # Save the resultant state
    #=
    jldopen("$output_folder/DEBUG_filament_state.jld2", "w") do file
        file["u_new"] = u_new
        file["v_new"] = v_new
        file["w_new"] = w_new
        file["b_new"] = b_new
    end
    =#
end

# Insert additional initialisation code
@additional_post_init! model

# run for the rest of the simulation
@info "Running filament simulation"
simulation.stop_time += run_time
run!(simulation)