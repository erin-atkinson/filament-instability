@inline dfm(a) = Field(Average(a; dims=2))

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

# Actually, I want to compute the streamfunction here too
ψ = Field(@at (Center, Nothing, Center) CumulativeIntegral(-u_dfm; dims=3))

outputs = (; u_dfm, v_dfm, w_dfm, b_dfm, ψ)

#update_outputs!(outputs) = map(compute!, (; mean_fields, outputs))

function update_outputs!(outputs)
    map(compute!, outputs)
end