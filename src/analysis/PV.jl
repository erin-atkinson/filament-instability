using Oceanostics
@inline dfm(a) = Field(Average(a; dims=2))

b_dfm = dfm(b)
v_dfm = dfm(b)
q = Field(@at (Center, Center, Center) KernelFunctionOperation{Face, Face, Face}(Oceanostics.FlowDiagnostics.ertel_potential_vorticity_fff, grid, u, v, w, b, 0, 0, sp.f))

mean_fields = (; b_dfm, v_dfm, q)

q_dfm = dfm(q)

q_vol_func(i, j, k, grid, q) = @inbounds q[i, j, k] <= 0 ? grid.Δxᶜᵃᵃ * grid.Δyᵃᶜᵃ * grid.Δzᵃᵃᶜ[k] : 0
q_vol_op = KernelFunctionOperation{Center, Center, Center}(q_vol_func, grid, q)

q_from_dfm_op = @at (Center, Nothing, Center) ∂z(b_dfm)*(∂x(v_dfm)+sp.f)-∂x(b_dfm)*∂z(v_dfm)
q_from_dfm = Field(q_from_dfm_op)

q_vol = Field(Reduction(sum!, q_vol_op; dims=(1, 2, 3)))

q_vol_from_dfm_op = KernelFunctionOperation{Center, Nothing, Center}(q_vol_func, grid, q_from_dfm)
q_vol_from_dfm = Field(Reduction(sum!, q_vol_from_dfm_op; dims=(1, 3)))

outputs = (; q_dfm, q_from_dfm, q_vol, q_vol_from_dfm)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, outputs)
end