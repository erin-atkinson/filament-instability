using Oceananigans
using JLD2
using CairoMakie
using Statistics
using ImageFiltering: imfilter, Kernel.gaussian

@inline function qᶜᶜᶜ(vᶜᶠᶜ, bᶜᶜᶜ, Δzᵃᵃᶜ, Δx, f)
    # Integrate
    
    ∂xbᶜᶜᶜ = (circshift(bᶜᶜᶜ, (-1, 0)) - circshift(bᶜᶜᶜ, (1, 0))) ./ (2Δx)
    ∂xvᶜᶠᶜ = (circshift(vᶜᶠᶜ, (-1, 0)) - circshift(vᶜᶠᶜ, (1, 0))) ./ (2Δx)

    ∂zbᶜᶜᶜ = (bᶜᶜᶜ - circshift(bᶜᶜᶜ, (0, 1))) ./ Δzᵃᵃᶜ
    ∂zbᶜᶜᶜ[:, 1] .= ∂zbᶜᶜᶜ[:, 2]

    ∂zvᶜᶠᶜ = (vᶜᶠᶜ - circshift(vᶜᶠᶜ, (0, 1))) ./ Δzᵃᵃᶜ
    ∂zvᶜᶠᶜ[:, 1] .= ∂zvᶜᶠᶜ[:, 2]
    return (∂xvᶜᶠᶜ .+ f) .* ∂zbᶜᶜᶜ .- ∂zvᶜᶠᶜ .* ∂xbᶜᶜᶜ
end


@inline function neg_pv_volume(foldername; σ=0)
    #foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    
    paramfilename = "parameters.jld2"
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen("$foldername/$filename")
    
    frames = keys(file["timeseries/t"])[101:2000]
    grid = file["serialized/grid"]
    zᵃᵃᶜ = znodes(grid, Center())
    zᵃᵃᶠ = znodes(grid, Face())
    xᶠᵃᵃ = xnodes(grid, Face())
    
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    ΔVᶜᵃᶜ = Δx .* repeat(Δzᵃᵃᶜ, length(xᶠᵃᵃ), 1)
    ts = [file["timeseries/t/$frame"] for frame in frames] .- 1
    Vq = map(frames) do frame
        v = file["timeseries/v_dfm/$frame"][:, 1, :]
        b = file["timeseries/b_dfm/$frame"][:, 1, :]
        (v, b) = map((v, b)) do field
            imfilter(field, gaussian((σ, 0)), "circular")
        end
        sum(ΔVᶜᵃᶜ[qᶜᶜᶜ(v, b, Δzᵃᵃᶜ, Δx, sp.f) .< 0])
    end
    close(file)
    return (; ts, Vq)
end

@inline function sc_by_pv_sign(foldername; σ=0)
    #foldername = "../scratch/filament-instability/$runname"
    filename = "down_front_mean.jld2"
    
    paramfilename = "parameters.jld2"
    sp = jldopen("$foldername/$paramfilename") do file
        file["parameters/simulation"]
    end
    
    file = jldopen("$foldername/$filename")
    
    frames = keys(file["timeseries/t"])[101:end-1]
    grid = file["serialized/grid"]
    zᵃᵃᶜ = znodes(grid, Center())
    zᵃᵃᶠ = znodes(grid, Face())
    xᶠᵃᵃ = xnodes(grid, Face())
    
    Δx = xᶠᵃᵃ[2] - xᶠᵃᵃ[1]
    Δzᵃᵃᶜ = reshape(diff(zᵃᵃᶠ), 1, length(zᵃᵃᶜ))
    ΔVᶜᵃᶜ = Δx .* repeat(Δzᵃᵃᶜ, length(xᶠᵃᵃ), 1)
    ts = [file["timeseries/t/$frame"] for frame in frames] .- 1
    U = map(frames) do frame
        # Interpolate
        u = file["timeseries/u_dfm/$frame"][:, 1, :]
        u = (circshift(u, (1, 0)) + circshift(u, (-1, 0))) / 2
        w = file["timeseries/w_dfm/$frame"][:, 1, :]
        w = (w[:, 2:end] + w[:, 1:end-1]) / 2
        v = file["timeseries/v_dfm/$frame"][:, 1, :]
        b = file["timeseries/b_dfm/$frame"][:, 1, :]
        (v, b) = map((v, b)) do field
            imfilter(field, gaussian((σ, 0)), "circular")
        end
        
        (u, w) = map((u, w)) do field
            imfilter(field, gaussian((σ, 0)), "circular")
        end
        U = u .^2 .+ w .^2
        pinds, ninds = let q = qᶜᶜᶜ(v, b, Δzᵃᵃᶜ, Δx, sp.f); (q .>= 0, q .< 0) end
        
        (sum(U[pinds] .* ΔVᶜᵃᶜ[pinds]) / sum(ΔVᶜᵃᶜ[pinds]), sum(U[ninds] .* ΔVᶜᵃᶜ[ninds]) / sum(ΔVᶜᵃᶜ[ninds]))
    end
    Upos = map(x->x[1], U)
    Uneg = map(x->x[2], U)
    close(file)
    return (; ts, Upos, Uneg)
end


@inline function plot_Vq!(layout_cell, runnames; σ=0)
    plot_datas = map(rn->neg_pv_volume(rn; σ), runnames)
    
    axis_kwargs = (;
        xlabel="t / 2π",
        ylabel=L"V_{q < 0}"
    )
    
    ax = Axis(layout_cell; axis_kwargs...)
    αs = reverse(range(0.3, 1.0, length(runnames)))
    map(zip(plot_datas, αs)) do (plot_data, α)
        lines!(ax, plot_data.ts ./ (2π), plot_data.Vq; color=(:red, α))
    end
end

@inline function plot_U!(layout_cell, runnames; σ=0)
    plot_datas = map(rn->sc_by_pv_sign(rn; σ), runnames)
    
    axis_kwargs = (;
        xlabel="t / 2π",
        ylabel=L"U_{q < 0}"
    )
    
    ax = Axis(layout_cell; axis_kwargs...)
    αs = reverse(range(0.3, 1.0, length(runnames)))
    map(zip(plot_datas, αs)) do (plot_data, α)
        lines!(ax, plot_data.ts, plot_data.Upos; color=(:red, α), linestyle=:dash)
        lines!(ax, plot_data.ts, plot_data.Uneg; color=(:red, α))
    end
end

@inline function frontogenesis_pv(runnames, legendtitle, runlabels; σ=0, size=(1000, 500))
    fig = Figure(; size)
    lns = plot_Vq!(fig[1, 1], runnames; σ)
    #lns = plot_U!(fig[1, 2], runnames)
    Legend(fig[1, 2], lns, runlabels, legendtitle)
    colgap!(fig.layout, 15)
    fig
end