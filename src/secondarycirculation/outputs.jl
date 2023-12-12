#=
secondarycirculation/outputs.jl
    Additional fields to output, none in this case
=#
using Oceananigans
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, ∇_dot_qᶜ
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators

# This macro should define and create output writers for any desired saved fields
dfm(a) = Field(Average(a; dims=2))
macro additional_outputs!(simulation)
    return quote
        model_fields = merge(model.velocities, model.tracers, model.auxiliary_fields)
        # Terms in the down-front vorticity evolution equation
        ω₁ = ∂y(w) - ∂z(v)
        ω₂ = ∂z(u) - ∂x(w)
        ζ = ∂x(v) - ∂y(u)
        ω_dfm = dfm(ω₂)
        
        # Vortex tilting terms
        tiltx = dfm(ω₁) * dfm(∂x(v))
        tiltz = dfm(ζ) * dfm(∂z(v))
        tiltx′ = dfm(ω₁ * ∂x(v)) - tiltx
        tiltz′ = dfm(ζ * ∂z(v)) - tiltz
        
        # Non-thermal-wind torque
        ATW = dfm(sp.f * ∂z(v) - ∂x(b))
        
        # Torques from friction, gonna not split them
        τ₃₁ = dfm(KernelFunctionOperation{Center, Center, Center}(ℑxzᶜᵃᶜ, grid, ∂zᶠᶜᶠ, ∂ⱼ_τ₁ⱼ,
            model.closure, model.diffusivity_fields, model.clock, model_fields, model.buoyancy))
        τ₁₃ = -dfm(KernelFunctionOperation{Center, Center, Center}(ℑxzᶜᵃᶜ, grid, ∂xᶠᶜᶠ, ∂ⱼ_τ₃ⱼ,
            model.closure, model.diffusivity_fields, model.clock, model_fields, model.buoyancy))
        
        # Eddy torques?
        ωuFLUX = -∂x(dfm((ω₂ - ω_dfm) * (u - u_dfm)))
        ωwFLUX = -∂z(dfm((ω₂ - ω_dfm) * (w - w_dfm)))
        
        $simulation.output_writers[:secondary_circulation_evolution] = JLD2OutputWriter(
            model,
            (; ω_dfm, tiltx′, tiltz′, ATW, τ = -τ₃₁ - τ₁₃, ωuFLUX, ωwFLUX),
            filename = "$output_folder/secondary_circulation.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        @info "Created secondary circulation writer"
        # Down-front momentum
        # Friction
        τ₂ = -dfm(KernelFunctionOperation{Center, Center, Center}(ℑyᵃᶜᵃ, grid, ∂ⱼ_τ₂ⱼ,
            model.closure, model.diffusivity_fields, model.clock, model_fields, model.buoyancy))
        
        # Eddy torques?
        vuFLUX = -∂x(dfm((v - v_dfm) * (u - u_dfm)))
        vwFLUX = -∂z(dfm((v - v_dfm) * (w - w_dfm)))
        
        $simulation.output_writers[:down_front_evolution] = JLD2OutputWriter(
            model,
            (; τ₂, vuFLUX, vwFLUX),
            filename = "$output_folder/down_front.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        @info "Created down-front momentum writer"
        # Across-front momentum
        # Friction
        τ₁ = -dfm(KernelFunctionOperation{Center, Center, Center}(ℑxᶜᵃᵃ, grid, ∂ⱼ_τ₁ⱼ,
            model.closure, model.diffusivity_fields, model.clock, model_fields, model.buoyancy))
        
        # Eddy torques?
        uuFLUX = -∂x(dfm((u - u_dfm) * (u - u_dfm)))
        uwFLUX = -∂z(dfm((u - u_dfm) * (w - w_dfm)))
        
        pHY′ = dfm(model.pressures.pHY′)
        pNHS = dfm(model.pressures.pNHS)
        
        $simulation.output_writers[:across_front_evolution] = JLD2OutputWriter(
            model,
            (; τ₁, uuFLUX, uwFLUX, φH₁=-∂x(pHY′), φNH₁=-∂x(pNHS)),
            filename = "$output_folder/across_front.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        @info "Created across-front momentum writer"
        # Vertical momentum
        # Friction
        τ₃ = -dfm(KernelFunctionOperation{Center, Center, Center}(ℑzᵃᵃᶜ, grid, ∂ⱼ_τ₃ⱼ,
            model.closure, model.diffusivity_fields, model.clock, model_fields, model.buoyancy))
        
        # Eddy torques?
        wuFLUX = -∂x(dfm((w - w_dfm) * (u - u_dfm)))
        wwFLUX = -∂z(dfm((w - w_dfm) * (w - w_dfm)))
        
        pHY′ = dfm(model.pressures.pHY′)
        pNHS = dfm(model.pressures.pNHS)
        
        $simulation.output_writers[:vertical_evolution] = JLD2OutputWriter(
            model,
            (; τ₃, wuFLUX, wwFLUX, φH₃=-∂z(pHY′), φNH₃=-∂z(pNHS)),
            filename = "$output_folder/vertical.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        @info "Created vertical momentum writer"
        
        # Buoyancy
        buFLUX = -∂x(dfm((b - b_dfm) * (u - u_dfm)))
        bwFLUX = -∂z(dfm((b - b_dfm) * (w - w_dfm)))
        val_tracer_index = Val(1)
        B = -dfm(KernelFunctionOperation{Center, Center, Center}(∇_dot_qᶜ, grid,
            model.closure, model.diffusivity_fields, val_tracer_index, b, model.clock, model_fields, model.buoyancy))
        
        $simulation.output_writers[:buoyancy_evolution] = JLD2OutputWriter(
            model,
            (; B, buFLUX, bwFLUX),
            filename = "$output_folder/buoyancy.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        @info "Created buoyancy writer"
    end
end

#= 
return quote
        # Define the lateral and geostrophic shear production and buoyancy flux
        
        LSP = Field(Average((u - u_dfm) * (v - v_dfm); dims=2)) * ∂x(v_dfm)
        GSP = Field(Average((w - w_dfm) * (v - v_dfm); dims=2)) * ∂z(v_dfm)
        
        uFLUX = Average((w - w_dfm) * (u - u_dfm); dims=2)
        vFLUX = Average((w - w_dfm) * (v - v_dfm); dims=2)
        wFLUX = Average((w - w_dfm) * (w - w_dfm); dims=2)
        bFLUX = Average((w - w_dfm) * (b - b_dfm); dims=2)
        
        $simulation.output_writers[:turbulent_flux] = JLD2OutputWriter(
            model,
            (; bFLUX, uFLUX, vFLUX, wFLUX),
            filename = "$output_folder/turbulent_flux.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        
        $simulation.output_writers[:shearproduction] = JLD2OutputWriter(
            model,
            (; LSP, GSP),
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
        
        pHY′ = Average(model.pressures.pHY′; dims=2)
        pNHS = Average(model.pressures.pNHS; dims=2)
        $simulation.output_writers[:pressure] = JLD2OutputWriter(
            model,
            (; pHY′, pNHS),
            filename = "$output_folder/pressure.jld2",
            schedule = TimeInterval(1/(sp.f*write_freq)),
            overwrite_existing = true
        )
        
        $simulation.output_writers[:full_state] = JLD2OutputWriter(
            model,
            (; u, v, w, b),
            filename = "$output_folder/full_state.jld2",
            schedule = TimeInterval(1/(sp.f)),
            overwrite_existing = true
        )
    end
=#