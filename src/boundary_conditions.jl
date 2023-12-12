#=
boundary_conditions.jl
    Set surface gradient boundary conditions for b and v
=#

using Oceananigans

@inline function get_boundary_conditions(simulation_parameters, init_time)
    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(sp.B₀))
    return (; b=b_bcs)
end