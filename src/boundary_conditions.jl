#=
boundary_conditions.jl
    Set surface gradient boundary conditions for b and v
=#

using Oceananigans

include("filament_state.jl")

@inline function get_boundary_conditions(simulation_parameters, init_time)
    (b, v) = get_filament_state(simulation_parameters; verbose=false)
    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(sp.B₀))
    #v_bc_func(x, y, t) = t <= init_time ? 0 : (v(x+5e-7, 0) - v(x-5e-7, 0)) / 1e-6
    #v_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(v_bc_func))
    return (; b=b_bcs)
end