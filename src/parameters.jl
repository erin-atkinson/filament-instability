#=
    Input parameters are non-dimensional numbers, here any dependent parameters are derived
    Ro=1: Rossby number
    Fr₀=1: Deep Froude number
    Frb=10: Boundary layer Froude number
    Ek=nothing: (Turbulent) Ekman number, for closure. If nothing then a 
    Pr=1: Prandtl number for buoyancy diffusion
    α=0.1: Boundary layer aspect ratio
    λ=0.05: Fractional width of transition to deep water
    δ=-0.25: Fractional height change of transition to deep water across filament
=#
default_inputs = (; Ro=1, Fr₀=1, Frb=10, Ek=nothing, Pr=1, α=0.1, λ=0.05, δ=-1/4)

function create_simulation_parameters(input_parameters; verbose=true)
    ip = (; default_inputs..., input_parameters...)
    simulation_parameters = let Ro=ip.Ro, Fr₀=ip.Fr₀, Frb=ip.Frb, Ek=ip.Ek, Pr=ip.Pr, α=ip.α, λ=ip.λ, δ=ip.δ
        # Impose length and frequency scales
        # Distance between fronts
        L = 1
        # Coriolis parameter
        f=1
        
        # Derived values
        # Height of the boundary layer
        H = α * L
        # Front height
        δH = δ * H
        # Vorticity scale
        ζ = Ro * f
        # Buoyancy frequencies
        N₀ = ζ / Fr₀
        Nb = ζ / Frb
        # Front width is then imposed by desired Rossby number
        # this value is an estimate, valid for small δ, α, do a better deriv
        ℓ = sqrt(0.560 * H * abs(δH) * abs(N₀^2 - Nb^2) / (f^2 * Ro))
        # Turbulent viscosity scale
        ν = Ek == nothing ? nothing : Ek * f * H^2
        # These aren't physical parameters but will be necessary for later
        # Defined here to avoid hard-coding elsewhere...
        # Simulation vertical extent
        Lz = 2.5H
        verbose && @info "Created simulation parameters\
            \nInput:\n Ro=$Ro\n Fr₀=$Fr₀\n Frb=$Frb\n Ek=$Ek\n α=$α\n λ=$λ\n δ=$δ\
        \nOutput:\n L=$L\n f=$f\n H=$H\n δH=$δH\n ζ=$ζ\n N₀=$N₀\n Nb=$Nb\n ℓ=$ℓ\n ν=$ν"
        simulation_parameters = (; Ro, Fr₀, Frb, Ek, α, λ, δ, L, f, H, δH, ζ, N₀, Nb, ℓ, ν)
    end
end
