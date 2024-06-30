using Oceananigans.Operators

# All of the not-nice terms in the streamfunction evolution equation\
# All on C N C
# Some are first-order time derivatives, these are all computed then
# differentiated later.


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
const Pr = sp.Pr

@inline dfm(a) = Field(Average(a; dims=2))

# Perturbation fields

u_bar = dfm(u)
v_dfm = dfm(v)
v_bar = dfm(v - v₀)
w_bar = dfm(w)
b_bar = dfm(b - b₀)
b_dfm = dfm(b)

mean_fields = (; u_bar, v_bar, w_bar, b_bar, v_dfm, b_dfm)

# Due to Oceananigans being weird when I use a unreduced-reduced field I need to do it like this
uu = dfm(u*u)
wu = dfm(w*u)

uv = dfm(u*v)
wv = dfm(w*v)

ww = dfm(w*w)

ub = dfm(u*b)
wb = dfm(w*b)

var_fields = (; uu, wu, uv, wv, ww, ub, wb)

w′u′ = Field(@at (Center, Nothing, Center) wu - w_bar * u_bar)
u′u′ = Field(@at (Face, Nothing, Face) uu - u_bar * u_bar)

w′v′ = Field(@at (Center, Nothing, Center) wv - w_bar * v_dfm)
u′v′ = Field(@at (Face, Nothing, Face) uv - u_bar * v_dfm)

w′w′ = Field(@at (Face, Nothing, Face) ww - w_bar * w_bar)

w′b′ = Field(@at (Face, Nothing, Face) wb - w_bar * b_dfm)
u′b′ = Field(@at (Center, Nothing, Center) ub - u_bar * b_dfm)

fluc_fields = (; w′u′, u′u′, w′v′, u′v′, w′w′, u′b′, w′b′)

# Kernel functions 
@inline ab(i, j, k, grid, a, b) = @inbounds a[i, j, k] * b[i, j, k]

# Non-linear (advection) terms
@inline function F_nl_func(i, j, k, grid, u, w, v, b)
    vᶠⁿᶜ = ℑxᶠᵃᵃ(i, j, k, grid, v)
    vᶜⁿᶠ = ℑzᵃᵃᶠ(i, j, k, grid, v)
    bᶠⁿᶜ = ℑxᶠᵃᵃ(i, j, k, grid, b)
    bᶜⁿᶠ = ℑzᵃᵃᶠ(i, j, k, grid, b)
    
    v1 = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂xᶜᶜᶜ, ab, u, vᶠⁿᶜ)
    v2 = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂zᶜᶜᶜ, ab, w, vᶜⁿᶠ)
    
    b1 = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂xᶜᶜᶜ, ab, u, bᶠⁿᶜ)
    b2 = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂zᶜᶜᶜ, ab, w, bᶜⁿᶠ)
    
    return @inbounds -f*v1[i, j, k] - f*v2[i, j, k] + b1[i, j, k] + b2[i, j, k]
end
# Non-linear (advection) terms that need to be differentiated later
@inline function t_F_nl_func(i, j, k, grid, u, w)
    
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    
    u1 = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂xᶜᶜᶜ, ab, u, u)
    u2 = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ab, uᶜᶜᶜ, wᶜᶜᶜ)
    
    w1 = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ab, uᶜᶜᶜ, wᶜᶜᶜ)
    w2 = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂xᶜᶜᶜ, ab, w, w)
    
    return @inbounds -u1[i, j, k,] - u2[i, j, k] + w1[i, j, k] + w2[i, j, k]
end
# turbulent production terms (that aren't a la VSP)
@inline function F_other_func(i, j, k, grid, u′v′, u′b′)
    
    a1 = ∂zᶜᶜᶜ(i, j, k, grid, ∂xᶜᶜᶠ, u′v′)
    a2 = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, u′b′)
    
    return @inbounds -f*a1[i, j, k] + a2[i, j, k]
end
# turbulent production terms that need to be differentiated later
@inline function t_F_other_func(i, j, k, grid, w′u′, u′u′, w′w′)
    
    a1 = ∂zᶜᶜᶜ(i, j, k, grid, ∂xᶜᶜᶠ, u′u′)
    a2 = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, w′u′)
    a3 = ∂zᶜᶜᶜ(i, j, k, grid, ∂xᶜᶜᶠ, w′w′)
    a4 = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, w′u′)
    
    return @inbounds -a1[i, j, k] - a2[i, j, k] + a3[i, j, k] + a4[i, j, k]
end
# I seperate the w'v' "VSP" term for illustrative reasons
@inline function F_VSP_func(i, j, k, grid, w′v′)
    
    a = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, w′v′)
    
    return @inbounds -f*a[i, j, k]
end

@inline function F_BFLUX_func(i, j, k, grid, w′b′)
    
    a = ∂zᶜᶜᶜ(i, j, k, grid, ∂xᶜᶜᶠ, w′b′)
    
    return @inbounds a[i, j, k]
end
# Sub-grid-scale viscosity
# This doesn't work and I don't know why
@inline function F_sgs_func(i, j, k, grid, ν, v, b)
    
    ∂xv = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)
    ∂zv = ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᶜᶠᶠ, v)
    
    ∂z∂xν∂xv = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, ∂xᶠᶜᶜ, ab, ν, ∂xv)
    ∂z∂zν∂zv = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ab, ν, ∂zv)
    
    ∂xb = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    ∂zb = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
    
    ∂x∂xν∂xb = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ab, ν, ∂xb)
    ∂x∂zν∂zb = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, ∂zᶜᶜᶠ, ab, ν, ∂zb)
    
    return @inbounds f*∂z∂xν∂xv[i, j, k] + f*∂z∂zν∂zv[i, j, k] - (∂x∂xν∂xb[i, j, k] + ∂x∂zν∂zb[i, j, k]) / Pr
end

@inline function t_F_sgs_func(i, j, k, grid, ν, u, w)
    
    ∂xu = ∂xᶜᶜᶜ(i, j, k, grid, u)
    ∂zu = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u)
    
    ∂z∂xν∂xu = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, ∂xᶠᶜᶜ, ab, ν, ∂xu)
    ∂z∂zν∂zu = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ab, ν, ∂zu)
    
    ∂xw = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w)
    ∂zw = ∂zᶜᶜᶜ(i, j, k, grid, w)
    
    ∂x∂xν∂xw = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ab, ν, ∂xw)
    ∂x∂zν∂zw = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, ∂zᶜᶜᶠ, ab, ν, ∂zw)
    
    return @inbounds (∂z∂xν∂xu[i, j, k] + ∂z∂zν∂zu[i, j, k]) - (∂x∂xν∂xw[i, j, k] + ∂x∂zν∂zw[i, j, k])
end

# Outputs
F_nl = Field(KernelFunctionOperation{Center, Nothing, Center}(F_nl_func, grid, u_bar, w_bar, v_bar, b_bar))
t_F_nl = Field(KernelFunctionOperation{Center, Nothing, Center}(t_F_nl_func, grid, u_bar, w_bar))

#F_sgs = dfm(KernelFunctionOperation{Center, Center, Center}(F_sgs_func, grid, ν, v, b))
#t_F_sgs = dfm(KernelFunctionOperation{Center, Center, Center}(t_F_sgs_func, grid, ν, u, w))

F_other = Field(KernelFunctionOperation{Center, Nothing, Center}(F_other_func, grid, u′v′, u′b′, w′b′))
t_F_other = Field(KernelFunctionOperation{Center, Nothing, Center}(t_F_other_func, grid, w′u′, u′u′, w′w′))

F_VSP = Field(KernelFunctionOperation{Center, Nothing, Center}(F_VSP_func, grid, w′v′))
F_BFLUX = Field(KernelFunctionOperation{Center, Nothing, Center}(F_BFLUX_func, grid, w′b′))

outputs = (; F_VSP, F_BFLUX)
temp_outputs = (; F_other, t_F_other, F_nl, t_F_nl) # , F_sgs, t_F_sgs

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, var_fields)
    map(compute!, fluc_fields)
    map(compute!, temp_outputs)
    map(compute!, outputs)
end

# Custom script for processing-after-processing?
macro postpostprocess()
    @info "Post-post-processing for time derivative of forcing terms"
    
    # continue at boundaries
    Δt = ts[2] - ts[1]
    jldopen(temp_filename, "r") do temp_file
        
        @inline F_nl(iter) = temp_file["timeseries/F_nl/$iter"]
        @inline t_F_nl(iter) = temp_file["timeseries/t_F_nl/$iter"]

        #@inline F_sgs(iter) = temp_file["timeseries/F_sgs/$iter"] 
        #@inline t_F_sgs(iter) = temp_file["timeseries/t_F_sgs/$iter"] 

        @inline F_other(iter) = temp_file["timeseries/F_other/$iter"] 
        @inline t_F_other(iter) = temp_file["timeseries/t_F_other/$iter"]
        
        for i in 1:length(iters)
            p_iter = iters[max(1, i-1)]
            iter = iters[i]
            n_iter = iters[min(length(iters), i+1)]
            jldopen(filename, "a") do file
                file["timeseries/F_nl/$iter"] = F_nl(iter) + (t_F_nl(n_iter) - t_F_nl(p_iter)) / (2Δt)
                #file["timeseries/F_sgs/$iter"] = F_sgs(iter) + (t_F_sgs(n_iter) - t_F_sgs(p_iter)) / (2Δt)
                file["timeseries/F_other/$iter"] = F_other(iter) + (t_F_other(n_iter) - t_F_other(p_iter)) / (2Δt)
            end
        end
    end
end
