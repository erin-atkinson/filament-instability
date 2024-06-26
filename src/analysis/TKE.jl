@inline dfm(a) = Field(Average(a; dims=2))

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

wv = dfm(w*v)
wb = dfm(w*b)
uv = dfm(u*v)

fluc_fields = (; wv, wb, uv)

w′v′ = wv - w_dfm * v_dfm
w′b′ = wb - w_dfm * b_dfm
u′v′ = uv - u_dfm * v_dfm

VSP = Field(@at (Center, Nothing, Center) -w′v′ * ∂z(v_dfm))
LSP = Field(@at (Center, Nothing, Center) -u′v′ * ∂x(v_dfm))
BFLUX = Field(@at (Center, Nothing, Center) w′b′)

outputs = (; VSP, LSP, BFLUX)

#update_outputs!(outputs) = map(compute!, (; mean_fields, outputs))

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, fluc_fields)
    map(compute!, outputs)
end