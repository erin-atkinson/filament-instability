#=
default/tracers.jl
    Definition, forcing and boundary conditions of additional tracers
=#
using Oceananigans

@inline function additional_tracers(simulation_parameters, mode_parameters)
    return (; )
end

@inline function additional_tracer_bcs(simulation_parameters, mode_parameters)
    return (; )
end

@inline function additional_tracer_forcing(simulation_parameters, mode_parameters)
    return (; )
end