#=
boundary_conditions.jl
    Set surface gradient boundary conditions for b and v
=#

using Oceananigans

FieldBoundaryConditions(top = GradientBoundaryCondition((x, y, t)->(v_front(x, 0) - v_front(x, -0.001))/0.002))

@inline function get_boundary_conditions(simulation_parameters)
    (b, v) = get_filament_state(simulation_parameters; verbose=false)
    b_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(simulation_parameters.Nb/10))
    v_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition((x, y, t)->(v(x, 0) - v(x, -1e-4))/1e-4))
    return (; b=b_bcs, v=v_bcs)
end