#=
sponge_layer.jl
    Create a forcing function that relaxes the fields to zero in the bottom σ of the domain
    at a rate c * internal wave freq that decays quadratically away from the boundary
=#

using Oceananigans

@inline function sponge_layer_damping_profile(simulation_parameters; σ=0.2)
    (x, y, z) -> let a = z / simulation_parameters.Lz
        a > (σ-1) ? 0 : ( (a + 1 - σ) / σ)^2
    end
end

@inline function sponge_layer_damping(simulation_parameters; σ=0.2, c=0.2, target=0)
    rate = c * simulation_parameters.N₀ / 2π
    mask=sponge_layer_damping_profile(simulation_parameters; σ)
    Relaxation(; rate, mask, target)
end