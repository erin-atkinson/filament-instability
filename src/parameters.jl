#=
parameters.jl
    Input parameters are non-dimensional numbers, here any dependent parameters are derived
    Ro=1: Maximum magnitude of Rossby number
    Ri=0.6: Minimum Richardson number
    Ek=nothing: (Turbulent) Ekman number, for closure. If nothing then a Smagorinsky-Lilly default is used
    Pr=1: Prandtl number for buoyancy diffusion
    α=0.25: Front width - separation ratio
    λ=0.05: Fractional width of transition to deep water
    δ=-0.25: Fractional height change of transition to deep water across filament

    λ ≪ δ ≪ 1 is required for the correct computation of parameters. (factors of two is enough)
=#
using SpecialFunctions
# Filament shape function
γ(s, δ, α) = -1 + (δ / 2) * (erf(s + 1/(2α)) - erf(s - 1/(2α)))
# Function to maximise for Ro
square_curvature(x, δ, α) = 2*abs((γ(x, δ, α)^2 - (γ(x+1e-5, δ, α)^2 + γ(x-1e-5, δ, α)^2) / 2) / 1e-10)
# Function to maximise for Ri
grad_squared(x, δ, α) = ((γ(x+1e-5, δ, α) - γ(x-1e-5, δ, α)) / (2e-5))^2

default_inputs = (; Ro=1, Ri=0.6, Frb=0.1, Ek=nothing, Pr=1, α=1/4, λ=0.05, δ=-1/4, β=0.1, cool=0.01)

@inline function create_simulation_parameters(input_parameters=(; ); verbose=true)
    ip = (; default_inputs..., input_parameters...)
    let Ro=ip.Ro, Ri=ip.Ri, Ek=ip.Ek, Pr=ip.Pr, α=ip.α, λ=ip.λ, δ=ip.δ, β=ip.β, cool=ip.cool
        # Setting variables
        # Distance between fronts
        L = 1
        # Height of the boundary layer
        H = β*L
        # Coriolis parameter
        f=1

        # Derived variables
        # Front width
        ℓ = α * L
        # Front height
        δH = δ * H
        # Difference in stratification
        ΔN² = - (2Ro * f^2 * ℓ^2) / (H^2 * maximum([square_curvature(x, δ, α) for x in range(-3L, 3L, 1000)])) 
        # Boundary layer stratification
        Nb² = ((Ri * ΔN²^2 * H^2) / (f^2 * ℓ^2)) * maximum([grad_squared(x, δ, α) for x in range(-3L, 3L, 1000)])
        # Turbulent viscosity scale
        ν = Ek == nothing ? nothing : (Ek * f * H^2)
        Nb = sqrt(Nb²)
        N₀ = sqrt(Nb² - ΔN²)
        # These aren't physical parameters but will be necessary for later
        # Defined here to avoid hard-coding elsewhere...
        # Simulation vertical extent
        Lz = 2.5H
        verbose && @info "Created simulation parameters\
            \nInput:\n Ro=$Ro\n Ri=$Ri\n Ek=$Ek\n α=$α\n λ=$λ\n δ=$δ\
        \nOutput:\n L=$L\n f=$f\n H=$H\n δH=$δH\n N₀=$N₀\n Nb=$Nb\n ℓ=$ℓ\n ν=$ν\n Lz=$Lz"
        (; Ro, Ri, Ek, α, λ, δ, L, f, H, δH, N₀, Nb, ℓ, ν, Lz, β, cool)
    end
end

@inline function create_simulation_parameters(; verbose=true, input_parameters...)
    create_simulation_parameters(input_parameters; verbose)
end
    