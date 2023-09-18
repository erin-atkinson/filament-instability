#=
    Functions that define the nominal buoyancy state and associated thermal wind
=#


# Transition function between the two stratifications
# Must satisfy g(∞) → x and g(-∞) → 0 and g''(x) ∼ 1
# i.e. it doesn't happen slowly...
@inline g(x) = (x + log(2*cosh(x)))/2

# Front profile
@inline h(x) =  (erf(x) + 1) /2
# Derivative
@inline ∂xh(x) = exp(-x^2) / sqrt(π)

# You'd want this to be a macro if these functions weren't just used for setup
function get_filament_state(simulation_parameters; verbose=true)
    sp = simulation_parameters
    let δH=sp.δH, L=sp.L, ℓ=sp.ℓ, H=sp.H, N₀=sp.N₀, Nb=sp.Nb, λ=sp.λ, f=sp.f
        filament_h(x) = δH * (h((x + L/2)/ℓ) - h((x - L/2)/ℓ)) - H
        filament_∂xh(x) = δH * (∂xh((x + L/2)/ℓ) - ∂xh((x - L/2)/ℓ)) / ℓ
        b(x, z) = N₀^2 * z + λ * (Nb^2 - N₀^2) * g((z - filament_h(x)) / λ)
        v(x, z) =  -(λ / f) * filament_∂xh(x) * (Nb^2 - N₀^2) * g((z - filament_h(x)) / λ)
        if verbose
            ζ_max = ζ_bar(simulation_parameters; v)
            @info "Filament state created: Ro=$(ζ_max/f), Fr₀=$(ζ_max/N₀), Frb=$(ζ_max/Nb)"
        end
        (; b, v)
    end
end

function ζ_bar(simulation_parameters; v)
    # Maximum vertical vorticity
    @inline ζ(x) = ( v(x+5e-4, 0) - v(x-5e-4, 0) ) / 1e-3
    xs = range(-simulation_parameters.L, simulation_parameters.L, 1000)
    return maximum(ζ.(xs))
end