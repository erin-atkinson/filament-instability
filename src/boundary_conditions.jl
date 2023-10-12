#=
boundary_conditions.jl
    Set surface gradient boundary conditions for b and v
=#

using Oceananigans

@inline function get_boundary_conditions(simulation_parameters, init_time)
    (b, v) = get_filament_state(simulation_parameters; verbose=false)
    v_grad(x) = (v(x, 1e-4) - v(x, -1e-4)) / 2e-4
    v_bc_func(x, y, t, p) = t >= p.init_time ? p.v_grad(x) : 0
    b_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(-simulation_parameters.cool * simulation_parameters.Nb^2))
    v_bcs = FieldBoundaryConditions(
        top = GradientBoundaryCondition(v_bc_func, parameters=(; init_time, v_grad))
    )
    return (; b=b_bcs, v=v_bcs)
end