#=
default/outputs.jl
    Additional fields to output, none in this case
=#
using Oceananigans

# This macro should define and create output writers for any desired saved fields
macro additional_outputs!(simulation)
    return quote
        # Define the lateral and geostrophic shear production and buoyancy flux
        
        BFLUX = Average((w - w_dfm) * (b - b_dfm); dims=2)
        
        LSP = Field(Average((u - u_dfm) * (v - v_dfm); dims=2)) * ∂x(v_dfm)
        
        GSP = Field(Average((w - w_dfm) * (v - v_dfm); dims=2)) * ∂z(v_dfm)
        
        MFLUX = Field(Average(w * v; dims=2)) - w_dfm * v_dfm
        
        $simulation.output_writers[:shearproduction] = JLD2OutputWriter(
            model,
            (; BFLUX, MFLUX, LSP, GSP),
            filename = "$output_folder/shear_production.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        
        q = (∂x(v_dfm) + sp.f) * ∂z(b_dfm) - ∂z(v_dfm) * ∂x(b_dfm)
        invRi = ∂z(v_dfm) * ∂z(v_dfm) / ∂z(b_dfm)
        
        $simulation.output_writers[:qRi] = JLD2OutputWriter(
            model,
            (; q, invRi),
            filename = "$output_folder/qRi.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        
        ω = Field(∂z(u_dfm) - ∂x(w_dfm))
        ω_depth = Average(ω * ω; dims=(1, 2))
        
        $simulation.output_writers[:secondarycirculation] = JLD2OutputWriter(
            model,
            (; ω, ω_depth),
            filename = "$output_folder/secondary_circulation.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        
    end
end