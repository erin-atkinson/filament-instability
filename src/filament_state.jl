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
    filament_h(x) = δH * (h((x + L/2)/ℓ) - h((x - L/2)/ℓ)) - H
    filament_∂xh(x) = δH * (∂xh((x + L/2)/ℓ) - ∂xh((x - L/2)/ℓ)) / ℓ
    b(x, z) = N₀^2 * z + λ * (Nb^2 - N₀^2) * g((z - filament_h(x)) / λ)
    v(x, z) =  -(λ / f) * filament_∂xh(x) * (Nb^2 - N₀^2) * g((z - filament_h(x)) / λ)
    
end

function Ro(simulation_paramters; v)
    # Maximum vertical vorticity
end

function Fr₀(simulation_paramters; v)
    # Maximum vertical vorticity
end