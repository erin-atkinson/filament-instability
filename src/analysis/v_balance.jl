using Oceananigans.Operators
# Need to grab the base state velocity
include("../filament_state.jl")
base_state = get_filament_state(sp; verbose=false)

v₀ = let
    xs = xnodes(grid, Center(); with_halos=true)
    zs = znodes(grid, Center(); with_halos=true)
    data = [base_state.v(x, z) for x in xs, y in 1:1, z in zs]
    Field{Center, Nothing, Center}(grid; data)
end 

@inline dfm(a) = Field(Average(a; dims=2))

u_bar = dfm(u)
v_bar = dfm(v - v₀)
w_bar = dfm(w)
b_bar = dfm(b)

mean_fields = (; u_bar, v_bar, w_bar, b_bar)

wv = dfm(w*v)
uv = dfm(u*v)

fluc_fields = (; wv, uv)

w′v′ = wv - w_bar * v_bar
u′v′ = uv - u_bar * v_bar

F_nl = Field(@at (Center, Nothing, Center) -∂x(u_bar * v_bar)-∂z(w_bar * v_bar))

F_VSP = Field(@at (Center, Nothing, Center) -∂x(u′v′))
F_LSP = Field(@at (Center, Nothing, Center) -∂z(w′v′))

uζ₀ = Field(@at (Center, Nothing, Center) u_bar*∂x(v₀))
wS₀ = Field(@at (Center, Nothing, Center) w_bar*∂z(v₀))

∇²v = ∂x(∂x(v)) + ∂y(∂y(v)) + ∂z(∂z(v))

# Need to make an operator for viscosity
@inline ab(i, j, k, grid, a, b) = @inbounds a[i, j, k] * b[i, j, k]

@inline function ν∇²vᶜᶜᶜ(i, j, k, grid, ν, v)
    
    ∂xv = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)
    ∂yv = ∂yᶜᶜᶜ(i, j, k, grid, v)
    ∂zv = ℑxyᶜᶜᵃ(i, j, k, grid, ∂zᶜᶠᶠ, v)
    
    ∂xν∂xv = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ab, ν, ∂xv)
    ∂yν∂yv = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ab, ν, ∂yv)
    ∂zν∂zv = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ab, ν, ∂zv)
    
    return @inbounds ∂xν∂xv[i, j, k] + ∂yν∂yv[i, j, k] + ∂zν∂zv[i, j, k]
end

F_sgs = dfm(KernelFunctionOperation{Center, Center, Center}(ν∇²vᶜᶜᶜ, grid, ν, v))

outputs = (; uζ₀, wS₀, F_nl, F_VSP, F_LSP, F_sgs)

#update_outputs!(outputs) = map(compute!, (; mean_fields, outputs))

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, fluc_fields)
    map(compute!, outputs)
end