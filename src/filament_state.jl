#=
filament_state.jl
    Functions that define the nominal buoyancy state and associated thermal wind
=#
using SpecialFunctions

# Transition function between the two stratifications
# Must satisfy g(∞) → x and g(-∞) → 0 and g''(x) ∼ 1
# i.e. it doesn't happen slowly...
g(s) = (s + log(2*cosh(s)))/2

# Filament shape function
γ(s, δ, α) = -1 + (δ / 2) * (erf(s + 1/(2α)) - erf(s - 1/(2α)))
∂γ∂s(s, δ, α) = (δ / sqrt(π)) * (exp(-(s + 1/(2α))^2) - exp(-(s - 1/(2α))^2))

# You'd want this to be a macro if these functions weren't just used for setup
@inline function get_filament_state(simulation_parameters; verbose=true)
    sp = simulation_parameters
    let δH=sp.δH, L=sp.L, ℓ=sp.ℓ, H=sp.H, N₀=sp.N₀, Nb=sp.Nb, λ=sp.λ, f=sp.f, δ=sp.δ, α=sp.α
        # Create filament profile
        filament_h(x) = H*γ(x/ℓ, δ, α)
        filament_∂xh(x) = H*∂γ∂s(x/ℓ, δ, α)/ℓ

        b(x, z) = N₀^2 * z + H* λ * (Nb^2 - N₀^2) * g((z - filament_h(x)) / (λ * H))
        v(x, z) =  -(H*λ / f) * filament_∂xh(x) * (Nb^2 - N₀^2) * g((z - filament_h(x)) / (λ * H))
        if verbose
            ζ_max = ζ_bar(simulation_parameters; v)
            ζ_min = ζ_bar(simulation_parameters; v=(x, z)->-v(x, z))
            @info "Filament state created: Ro=$(ζ_max/f), Ro_min=$(ζ_min/f), Fr₀=$(ζ_max/N₀), Frb=$(ζ_max/Nb), Ri_min=$(Ri_min(simulation_parameters; v, b))"
        end
        (; b, v)
    end
end



@inline function ζ_bar(simulation_parameters; v)
    # Maximum vertical vorticity
    @inline ζ(x) = ( v(x+5e-4, 0) - v(x-5e-4, 0) ) / 1e-3
    xs = range(-5simulation_parameters.L, 5simulation_parameters.L, 10000)
    return maximum(ζ.(xs))
end

@inline function Ri_min(simulation_parameters; v, b)
    # The minimum Richardson number
    xs = range(-simulation_parameters.L, simulation_parameters.L, 100)
    zs = range(-simulation_parameters.H, 0, 100)
    @inline ∂v∂z(x, z) =  (v(x, z+5e-4) - v(x, z-5e-4)) / 1e-3
    @inline ∂b∂z(x, z) =  (b(x, z+5e-4) - b(x, z-5e-4)) / 1e-3
    @inline Ri(x, z) = abs(∂b∂z(x, z)) / ∂v∂z(x, z)^2
    return minimum([Ri(x, z) for x in xs, z in zs])
end