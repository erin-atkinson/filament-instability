using Oceananigans.Operators

# This calculates all of the "nice" terms in the streamfunction evolution equation
# (no non-linear, turbulence, or sgs)
# doesn't compute time derivatives

# Need to grab the base state
include("../filament_state.jl")
base_state = get_filament_state(sp; verbose=false)

v₀, b₀ = let
    xs = xnodes(grid, Center(); with_halos=true)
    zs = znodes(grid, Center(); with_halos=true)
    v_data = [base_state.v(x, z) for x in xs, y in 1:1, z in zs]
    b_data = [base_state.b(x, z) for x in xs, y in 1:1, z in zs]
    Field{Center, Nothing, Center}(grid; data=v_data), Field{Center, Nothing, Center}(grid; data=b_data)
end 

const f = sp.f

@inline dfm(a) = Field(Average(a; dims=2))

u_dfm = dfm(u)
ψ = Field(@at (Center, Nothing, Center) CumulativeIntegral(-u_dfm; dims=3))

mean_fields = (; u_dfm, ψ)

# Kernel functions 
function Lψ_func(i, j, k, grid, v₀, b₀, ψ)
    ζ₀ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v₀)
    S₀ = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v₀)
    M₀² = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b₀)
    N₀² = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b₀)
    
    ψ_zz = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ψ)
    ψ_xz = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, ∂zᶜᶜᶠ, ψ)
    ψ_xx = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ψ)
    
    return @inbounds -f^2 * ψ_zz[i, j, k] - f * ζ₀[i, j, k] * ψ_zz[i, j, k] + f * S₀[i, j, k] * ψ_xz[i, j, k] + M₀²[i, j, k] * ψ_xz[i, j, k] - N₀²[i, j, k] * ψ_xx[i, j, k]
end
function ∇²ψ_func(i, j, k, grid, ψ)
    return ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ψ) + ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ψ)
end

Lψ = Field(KernelFunctionOperation{Center, Nothing, Center}(Lψ_func, grid, v₀, b₀, ψ))
∇²ψ = Field(KernelFunctionOperation{Center, Nothing, Center}(∇²ψ_func, grid, ψ))


outputs = (; Lψ)
temp_outputs = (; ∇²ψ)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, temp_outputs)
    map(compute!, outputs)
end

# Custom script for processing-after-processing?
macro postpostprocess()
    @info "Post-post-processing for time derivative of ψ"
    # time derivative tt of ψ
    # continue at boundaries
    Δt = ts[2] - ts[1]
    jldopen(temp_filename, "r") do temp_file
        @inline ∇²ψ(iter) = temp_file["timeseries/∇²ψ/$iter"]
        for i in 1:length(iters)
            p_iter = iters[max(1, i-1)]
            iter = iters[i]
            n_iter = iters[min(length(iters), i+1)]
            ∇²ψ_tt = -(2 / Δt^2) * (∇²ψ(iter) - (∇²ψ(n_iter) + ∇²ψ(p_iter)) / 2)
            jldopen(filename, "a") do file
                file["timeseries/∇²ψ_tt/$iter"] = ∇²ψ_tt
            end
        end
    end
end